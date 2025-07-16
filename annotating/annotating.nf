// Define a list of required parameters
def requiredParams = ['bam', 'output_dir', 'sample']

// Check each required parameter
requiredParams.each { param ->
    if (!params.containsKey(param) || params."$param" == '') {
        log.error "Missing required parameter: --$param"
        params.help = true
    }
}

// Default parameter values
def defaultParams = [
    // SelectReadsWithSignal defaults
    res_step: 100,
    thre_level: 0.1,
    num_points: 10,
    size_thre: 1,
    threads_select: 2,

    // PrepareAnnotation defaults
    reviewers: 'default',
    nb_random: 5000,
    nb_with_signal: 1000,
    min_read_length: 10000,
    batch_size: 50,
    path_discard: null,
    with_independent_review: true,
    target_mod_type: 'N+b?',
    data_type: 'Megalodon',
    threads_prepare: 8
]

// Apply default values
defaultParams.each { key, value -> params[key] = params.getOrDefault(key, value) }

workflow.onComplete {
    def duration = workflow.duration.toMinutes()
    def status = workflow.success ? "SUCCESS" : "FAILED"
    def status_emo1 = workflow.success ? "✅" : "⁉️"
    def status_emo2 = workflow.success ? "🎉" : "😬"

    log.info """
    ────────────────────────────────────────────────────
    ${status_emo1}  Workflow Execution Complete! ${status_emo2}
    
    🔹 Workflow Name : ${workflow.runName}
    🔹 Status        : ${status}
    🔹 Duration      : ${duration} minutes
    🔹 Output Dir    : ${params.output_dir}
    
    📢 If the workflow completed successfully, results are now available in:
       → ${params.output_dir} # Datasets need to be retrieve on computer with a browser for a manual annotation
       
       conda activate forkml_annotating

       Rscript app.R # Run shiny application
       # Open http://127.0.0.1:9999/ in browser (tested with Chrome)

    
    ❗ If there were errors, check the logs using:
       → nextflow log ${workflow.runName}
    ────────────────────────────────────────────────────
    """
}

//////////// PARAMETERS ////////////
// Help display
help = """
Pipeline parameters:
  --bam                     Input BAM file
  --output_dir              Directory where results will be stored
  --sample                  Sample name

Process: SelectReadsWithSignal
  --res_step X              Resolution step for signal filtering (default: 100)
  --thre_level X            Threshold level for signal detection (default: 0.1)
  --num_points X            Number of points to evaluate (default: 10)
  --size_thre X             Size threshold for filtering (default: 1)
  --threads_select X        Number of threads for computation (default: 2)

Process: PrepareAnnotation
  --reviewers X             List of reviewers (comma-separated, default: Alan,Benoit)
  --nb_random X             Number of randomly selected reads (default: 0)
  --nb_with_signal X        Number of signal-based selected reads (default: 1000)
  --min_read_length X       Minimum read length filter (default: 10000)
  --batch_size X            Number of reads per batch (default: 50)
  --threads_prepare X       Number of threads for processing (default: 8)
  --path_discard            Path to file with read_id to discard
  --with_independent_review Include independent review step (default: true)
  --data_type X             Data type (default: 'Megalodon', or 'DNAscent')
  --target_mod_type X       Modification type to target (default: N+b?)
"""

params.help = false
if (params.size() == 1 || params.help) {
    println help
    exit 1
}

// Print all params at the start of the script
log.info "Starting workflow with parameters:"
params.each { param, value ->
    log.info "${param}: ${value}"
}

process SelectReadsWithSignal {
    output:
        path "listids_thre${params.thre_level}_points_${params.num_points}.${params.sample}.txt" into filtered_read_id_files
    script:
        def path_bam = file(params.bam).toAbsolutePath().toString()
        """
        select_reads_with_signal \
            --in_bam ${path_bam} \
            --out_folder ./ \
            --base_name ${params.sample} \
            --res_step ${params.res_step} \
            --thre_level ${params.thre_level} \
            --num_points ${params.num_points} \
            --size_thre ${params.size_thre} \
            --threads ${params.threads_select} \
            |& tee
        """
}

process PrepareAnnotation {
    publishDir "${params.output_dir}/forkml_annotation", mode: 'copy', pattern: "dataset.subset.*"
    input:
        path read_ids from filtered_read_id_files
    output:
        path "dataset.subset.rds" into annotation_output
        path "dataset.subset.read_id.txt" into read_id_output
        path "dataset.subset.read_id.source.txt" into read_id_source_output
        path "dataset.subset.sorted.bam" into bam_output
        path "dataset.subset.sorted.bam.bai" into bam_index_output
        path "data_*" into reviewers_data

    script:
        def path_bam = file(params.bam).toAbsolutePath().toString()
        def with_independent_review = params.with_independent_review ? '--with_independent_review' : ''
        """
        prepare_annotation \
            -o ${params.output_dir} \
            -n ${params.sample} \
            -b ${path_bam} \
            -r ${params.reviewers} \
            -i $read_ids \
            --path_discard ${params.path_discard} \
            --nb_random ${params.nb_random} \
            --nb_with_signal ${params.nb_with_signal} \
            --min_read_length ${params.min_read_length} \
            --batch_size ${params.batch_size} \
            --nb_threads ${params.threads_prepare} \
            --target_mod_type ${params.target_mod_type} \
            --data_type ${params.data_type} \
            ${with_independent_review}
        mv ${params.output_dir}/${params.sample}_all_reviewers/* ./
        """
}

process SetupForkMLAnnotation {
    publishDir "${params.output_dir}", mode: 'copy'
    input:
        path dataset from reviewers_data.flatten()

    output:
        path "forkml_annotation/for_${dataset.name.replaceFirst('data_', '')}" into forkml_annotation_dir

    script:
        """
        # Extract reviewer name (remove "data_" prefix)
        base_output=\$(basename "$dataset" | sed 's/data_//')

        # Create annotation directory for this reviewer
        mkdir -p "forkml_annotation/for_\$base_output"

        # Copy template files
        cp -r ${workflow.projectDir}/template/* "forkml_annotation/for_\$base_output/"

        # Modify samples.csv
        sed -i "s/<name>/\${base_output}/g" "forkml_annotation/for_\$base_output/samples.csv"
        sed -i "s/<time>/30/g" "forkml_annotation/for_\$base_output/samples.csv"

        # Move reviewer data folder
        mv "$dataset" "forkml_annotation/for_\$base_output/"

        """
}
