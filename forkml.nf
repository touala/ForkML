// Define a list of required parameters
def requiredParams = ['bam', 'model', 'output_dir', 'sample']

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
    nb_reads_process: 2000,
    nb_process: -1,
    min_length: 10000,
    resolution: 'low',
    legacy: false,
    duration_pulse: 4,
    frequency_pulse: 30,
    size_margin: 3000,
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
  --bam                     Input mod_mapping.bam
  --model                   Model name (see in ./models/)
  --output_dir              Directory where to put the results
  --sample                  Sample name

Processes parameters:
  --nb_process              Number of processes in total
  --nb_reads_process        Number of reads per process
  --min_length X            Minimum length of the reads
  --resolution X            Signal subsampling resolution. Either "high" (1pt/20bp) or "low" (1pt/100bp)
  --legacy                  Run in legacy mode
  --duration_pulse X        Duration of the BrdU pulse, in minutes
  --frequency_pulse X       Time between pulses (minutes)
  --size_margin X           Size of the margin on both sides of the detected segments
"""
params.help = false
if (params.size()==1 || params.help) {
    println help
    exit 1
}

// Set main parameters
// Resolution
if (params.resolution=='low') {
    subsampling_step = 100
} else if (params.resolution=='high') {
    subsampling_step = 20
} else {
    println('ERROR : --resolution argument needs to be either "low" or "high".')
    exit 1
}
// Outputs
tracks_out_file  = 'tracks_'+params.sample+'.tsv.gz'
signals_out_file = 'signals_'+params.sample+'.tsv.gz'
analysis_output = params.output_dir+'/analysis_'+params.sample+'/'

// Print all params at the start of the script
log.info "Starting workflow with parameters:"
params.each { param, value ->
    log.info "${param}: ${value}"
}


//////////// PIPELINE ////////////

process SplitBam {
    output:
        path 'mod_splitted_*' into bam_chunks
    script:
        def path_bam = file(params.bam).toAbsolutePath().toString()
        """
        NUM_DIGITS=6
        NB_CPUS=40
        
        CPUS_FOR_VIEW=\$((NB_CPUS * 80 / 100))
        CPUS_FOR_FILTER=\$((NB_CPUS * 20 / 100))
        export CPUS_FOR_FILTER
        
        samtools view -H ${path_bam} > header
        samtools view -@ \$CPUS_FOR_VIEW ${path_bam} | split - mod_splitted_ -l $params.nb_reads_process -a \$NUM_DIGITS -d --filter='cat header - | samtools view -@ \$CPUS_FOR_FILTER -b - > \$FILE.bam'
        """
}

process ParseBam {
    //publishDir path: out_dir, mode: 'copy'
    input:
        path x from bam_chunks.flatten()
    output:
        set file("*.fa.gz"), file("*.fa_ratio.gz"), file("*.fa_cigar.gz") into reads
    script:
        """
        preprocess_modbam \
            --in_bam $x \
            --out_folder ./ \
            --name_files $params.sample
        """
}

// Transform Fork-seq data from .fa and .fa_ratio files to a R table as a .tsv.gz file
process PreprocessFasta {
    input:
        set file(reads_fa), file(reads_ra), file(reads_ca) from reads
    output:
        set val(id), file("reads_*.rds") optional true into mapping_info
    script:
        chunk = (reads_fa.name =~ /(\d+)\.\w+.\w+$/)[0][1]
        id = reads_fa.baseName + '.' + chunk
        legacy = params.legacy ? '--legacy' : ''
        """
        get_reads \
            -f $reads_fa \
            -r $reads_ra \
            -c $reads_ca \
            -o reads_${id}.rds \
            --min_length $params.min_length \
            --subsampling_step $subsampling_step \
            --sample $params.sample $legacy 
        """
}
mapping_info.into { mapping_info_for_segmentation; mapping_info_for_detection }

// Segment reads by direction of replication
process Segmentation {
    input:
        set val(id), file(reads) from mapping_info_for_segmentation
    output:
        set val(id), file("signals_seg_${id}.tsv.gz"), file("tracks_seg_${id}.tsv.gz") optional true into segmentation
    script:
        """
        export RETICULATE_PYTHON=\$(which python)
        segmentation \
            -i $reads \
            -t tracks_seg_${id}.tsv.gz \
            -s signals_seg_${id}.tsv.gz \
            --subsampling $subsampling_step \
            --model_name $params.model
        sleep 30
        """
}

segmentation
    .join(mapping_info_for_detection, by: 0)
    .map { entry ->
        def id      = entry[0]
        def reads   = entry[3]
        def signals = entry[1]
        def tracks  = entry[2]
        return tuple(id, reads, signals, tracks)
    }
    .set { detection_inputs }

// Detect replication tracks in Fork-seq data
process Detection {
    input:
        tuple val(id), file(reads), file(signals), file(tracks) from detection_inputs
    output:
        set val(id), file("tracks_fit_${id}.tsv.gz") optional true into tracks
        set val(id), file("signals_fit_${id}.tsv.gz") optional true into signals
    script:
        legacy = params.legacy ? '--legacy' : ''
        """
        detect \
            -a $signals \
            -b $tracks \
            -c $reads \
            -y signals_fit_${id}.tsv.gz \
            -z tracks_fit_${id}.tsv.gz \
            --duration_pulse $params.duration_pulse \
            --frequency_pulse $params.frequency_pulse \
            --subsampling_step $subsampling_step \
            --size_margin $params.size_margin \
            $legacy
        """
}

// Merge all tables of tracks together
process MergeTracks {
    publishDir path: params.output_dir, mode: 'copy'
    input:
        file all_tracks from tracks.map { it[1] }.collect()
    output:
        file "$tracks_out_file" into tracks_merged
    script:
        """
        zcat $all_tracks | awk 'NR > 1 && /^read_id/ { next } 1' | gzip > $tracks_out_file
        """
}

// Merge all list of signals together
process MergeSignals {
    publishDir path: params.output_dir, mode: 'copy'
    input:
        file all_signals from signals.map { it[1] }.collect()
    output:
        file "$signals_out_file" into signals_merged
    script:
        """
        zcat $all_signals | awk 'NR > 1 && /^read_id/ { next } 1' | gzip > $signals_out_file
        """
}

process ComputeScore {
    publishDir path: params.output_dir, mode: 'copy'
    output:
        file "scores_${params.sample}.txt" into scores_file
    script:
        def path_bam = file(params.bam).toAbsolutePath().toString()
        """
        compute_signal_score ${path_bam} scores_${params.sample}.txt --threads 2
        """
}

process InputStats {
    publishDir path: params.output_dir, mode: 'copy'
    output:
        file "inputStats_${params.sample}.txt" into stats_file
    script:
        def path_bam = file(params.bam).toAbsolutePath().toString()
        """
        samtools stats ${path_bam} | grep ^SN | cut -f 2- > inputStats_${params.sample}.txt
        """
}

// Produce basic analyses of ForkML results.
process AnalysisTracks {
    publishDir path: analysis_output, mode: 'copy'
    input:
        file tracks from tracks_merged
        file signals from signals_merged
        file scores from scores_file
        file stats from stats_file
    output:
        file '*' optional true
    script:
        legacy = params.legacy ? '--legacy' : ''
        """
        basic_analysis --sample ${params.sample} --input_dir . --output_dir ./ --pulse_duration ${params.frequency_pulse} --bin_size 20000 $legacy --plot
        """
}

