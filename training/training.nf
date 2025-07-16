// Define a list of required parameters
def requiredParams = ['input_bams', 'annotations', 'output']

// Check each required parameter
requiredParams.each { param ->
    if (!params.containsKey(param) || params."$param" == '') {
        // If a required parameter is not set or is an empty string, exit with an error message
        log.error "Missing required parameter: --$param"
        params.help = true
    }
}

// Default parameter values
def defaultParams = [
    nb_threads: 20,
    nb_reads_process: 2000,
    subsampling_step: 10,
    max_read_length: 200000,
    target_mod_type: 'N+b?',
    data_type: 'Megalodon',
    prop_train: 0.6,
    prop_valid: 0.2,
    prop_test: 0.2,
    noise_control: null,
    with_augmentation: false,
    with_fork_distinction: false,
    padding_method: 'clean',
    nb_layers: [4, 5, 6],
    nb_filters: [8, 16],
    size_kernel: [10],
    with_training: false,
    with_hp_tuning: false,
    with_gpus: false,
    tmpdir: null
]

// Apply default values
defaultParams.each { key, value ->
    params[key] = params.getOrDefault(key, value)
}

workflow.onComplete {
    log.info "${workflow.runName} complete"
}

//////////// PARAMETERS ////////////
// Help display
help = """
Pipeline parameters:
  --input_bams            Comma-separated BAM files to process
  --annotations           Glob path or list of RDS annotation files
  --output                Output directory

Processes parameters:
  --subsampling_step X    Make the subsampling at 1/(10*X)
  --prop_train X          Proportion of data for training
  --prop_valid X          Proportion of data for validation
  --prop_test X           Proportion of data for test
  --noise_control FILE    Path to noise control file with read_id
  --with_augmentation     Perform augmentation of the training set (true/false)
  --with_fork_distinction Reverse signal to distinguish leading/lagging (true/false)
  --padding_method X      Configure masking types for padding (clean/rough/default)
  --max_read_length X     Longest read to conserve before splitting it (bp)
  --data_type X           Data type (default: 'Megalodon', or 'DNAscent')
  --target_mod_type X     Modification type to target (default: N+b?)
  --nb_layers X           Number of downsampling layers in the U-Net network (default: 4,5,6)
  --nb_filters X          Number of filters per convolutional layer (default: 8,16)
  --size_kernel X         Kernel size (default: 10)
  --with_training         Perform model training (default: false)
  --with_hp_tuning        Perform model training (default: false)
  --with_gpus X           Use GPU or only CPU (true/false)
  --tmpdir X              Temporary directory for training only. Default to output directory
"""
params.help = false
if (params.size()==1 || params.help) {
    println help
    exit 1
}

// Convert array inputs if necessary
def nb_layers = params.nb_layers // Needs to modify parameter name donwstream too
if (nb_layers instanceof String) {
    nb_layers = nb_layers.tokenize(',').collect { it.trim().toInteger() }
}

// Calculate max layer for pooling
nb_pooling = nb_layers.max()

// Parameters extraction
if( params.tmpdir == null ) {
    tmp_dir = params.output
}else{
    tmp_dir = params.tmpdir
}

out_dir = java.nio.file.Paths.get(params.output).toString()

def input_bams = (params.input_bams instanceof String)
    ? params.input_bams.split(',').collect { 
        def it_trimmed = it.trim()
        it_trimmed.startsWith('/') ? it_trimmed : "${workflow.launchDir}/${it_trimmed}"
    }
    : params.input_bams.collect { 
        def path = it instanceof String ? it : it.toString()
        path.startsWith('/') ? path : "${workflow.launchDir}/${path}"
    }

def annotations_list = (params.annotations instanceof String) 
    ? params.annotations.split(',').collect { it.trim() }
    : params.annotations

def annotations_abs = annotations_list.collect { path ->
    path.startsWith('/') ? path : "${workflow.launchDir}/${path}"
}.join(',')

//////////// PIPELINE ////////////
process ProcessAnnotation {
    output:
        file "annotation.read_id.txt" into list_read_id
        file "annotation.rds" into annotations_selected
    script:
        def noise_opt = params.noise_control ? "--noise_control ${file(params.noise_control).toAbsolutePath()}" : ""
        """
        prepare_annotation_data -a "${annotations_abs}" -o annotation.rds ${noise_opt}
        """
}

Channel
    .value(input_bams)
    .set { bam_list }

process MergeBam {
    input:
        file list_id from list_read_id
        val bam_paths from bam_list

    output:
        path "merged.mod_mappings.sorted.bam" into merged_bam
        path "merged.mod_mappings.sorted.bam.bai"

    script:
        """
        samtools cat ${bam_paths.join(' ')} | samtools view -bN $list_id - | samtools sort -o merged.mod_mappings.sorted.bam -
        samtools index merged.mod_mappings.sorted.bam
        """
}

process SplitBam {
    input:
        file bam from merged_bam
    output:
        path 'mod_splitted_*' into bam_chunks
    script:
        """
        NUM_DIGITS=6
        NB_CPUS=${params.nb_threads}
        
        CPUS_FOR_VIEW=\$((NB_CPUS * 80 / 100))
        CPUS_FOR_FILTER=\$((NB_CPUS * 20 / 100))
        export CPUS_FOR_FILTER
        
        samtools view -H ${bam} > header
        samtools view -@ \$CPUS_FOR_VIEW ${bam} | split - mod_splitted_ -l $params.nb_reads_process -a \$NUM_DIGITS -d --filter='cat header - | samtools view -@ \$CPUS_FOR_FILTER -b - > \$FILE.bam'
        """
}

process AnnotateSignal {
    input:
        path x from bam_chunks.flatten()
        file annotation_table from annotations_selected
    output:
        set val(id), file("*.prepared.rds") optional true into data_chunks_prepared
    script:
        id = (x =~ /(\w+)\_(\w+)\_(\d+)/)[0][3]

        """
        prepare_training_data -b $x -a $annotation_table -o ./chunk_${id}.prepared.rds --target_mod_type ${params.target_mod_type} --data_type ${params.data_type}
        """
}

process SmoothSignal {
    input:
        set val(id), file(x) from data_chunks_prepared
    output:
        set val(id), file("*.smoothed.rds") into data_smooth_prepared
    script:
        """
        smooth_training_data -i $x -m by_bin -w 100 -o ./chunk_${id}.smoothed.rds
        """
}
// TODO maybe publishDir pdf.
// Could further clean up fork junctions and read ends from misclicks

// Prepare data for training
process PrepareData {
    publishDir path: out_dir, mode: 'copy'
    input:
        file all_smoothed from data_smooth_prepared.map { it[1] }.collect()
    output:
        tuple(file("train_signals_balanced_0xempty.npy"), 
            file("train_labels_balanced_0xempty.npy"),
            file("train_signals_balanced_3xempty.npy"), 
            file("train_labels_balanced_3xempty.npy"),
            file("train_signals_nonbalanced.npy"), 
            file("train_labels_nonbalanced.npy"), 
            file("valid_signals.npy"), 
            file("valid_labels.npy"),
        ) into training_arrays
        file '*'
    script:
        def smoothed_files = all_smoothed.join(",")
        def with_augmentation = params.with_augmentation ? '--with_augmentation' : ''
        def with_fork_distinction = params.with_fork_distinction ? '--with_fork_distinction' : ''

        """
        export RETICULATE_PYTHON=\$(which python)
        save_training_data \
            -i $smoothed_files \
            -o ./ \
            --prop_train $params.prop_train \
            --prop_valid $params.prop_valid \
            --prop_test $params.prop_test \
            --nb_pooling $nb_pooling \
            --padding_method $params.padding_method \
            --max_read_length $params.max_read_length \
            ${with_augmentation} ${with_fork_distinction}
        """
}

// Perform training (Optional)
if(params.with_training) {
    process Training {
        label (params.with_gpus ? 'with_gpus': 'with_cpus')
        env {
            if (!params.with_gpus) {
                CUDA_VISIBLE_DEVICES = '-1'
            }
        }
        publishDir path: out_dir, mode: 'move'
        input:
            tuple(file(train_signals_balanced_0xempty), 
                file(train_labels_balanced_0xempty),
                file(train_signals_balanced_3xempty), 
                file(train_labels_balanced_3xempty),
                file(train_signals_nonbalanced), 
                file(train_labels_nonbalanced), 
                file(valid_signals), 
                file(valid_labels)
            ) from training_arrays
        output:
            file '*'
        script:
            """
            CONDA_PREFIX=\$(dirname \$(dirname \$(which python)))
            export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH:-}:\$CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/
            training \
                $train_signals_balanced_0xempty \
                $train_labels_balanced_0xempty \
                $train_signals_balanced_3xempty \
                $train_labels_balanced_3xempty \
                $train_signals_nonbalanced \
                $train_labels_nonbalanced \
                $valid_signals \
                $valid_labels \
                $params.padding_method \
                $params.with_hp_tuning \
                ./ \
                ${tmp_dir} \
                $params.nb_threads
            """       
    }
}