# ForkML: DNA Replication Fork Detection Pipeline

This repository provides a main Nextflow-based pipeline for detecting DNA replication forks from modified BAM files (Nanopore sequencing R9 or R10). This main pipeline performs detection, filtering, and model inference. The complete ForkML toolkit contains:
- "Detection" pipeline with `nextflow run forkml.nf`, see below for full description.
- "Annotating" pipeline with `nextflow run annotating/annotating.nf` and associated Shiny web application for manual annotation of DNA replication tracks. See [README](https://github.com/touala/ForkML/blob/main/annotating/readme.md) for full description.
- "Training" pipeline with `nextflow run training/training.nf` for training models on DNA replication fork signal from double pulse experiments. See [README](https://github.com/touala/ForkML/blob/main/training/readme.md) for full description.

In `figures/`, you can find the script making all of the figures with data from our paper (see [here](https://github.com/touala/ForkML/blob/main/figures/readme.md)).

## 🚀 Quick Start

Make sure you install and activate the required conda environment (< 1h):

```bash
conda env create -f forkml.yml -n forkml 
conda activate forkml

# Get Nanopore R10.4.1 model, ~3 GB
wget -O ./models/default_R10.h5 https://zenodo.org/records/15706385/files/default_R10.h5
# Or get Nanopore R9.4.1 model, ~3 GB
wget -O ./models/default_R9.h5 https://zenodo.org/records/15706385/files/default_R9.h5
```

This includes all necessary R packages (`tidyverse`, `rsamtools`, etc.) and Python tools (TensorFlow, Pandas, etc.).

### Example Run

```bash
nextflow run forkml.nf -c params_local.config \
  --bam /path/to/mod_mappings.sorted.bam \
  --output_dir output_dir \
  --sample sample_name \
  --model default_R10.h5 \
  -resume
```

## 🧾 Parameters

| Parameter     | Description                                        |
|---------------|----------------------------------------------------|
| `--bam`       | Input BAM file (from Megalodon or DNAscent).        |
| `--model`     | Path to the pre-trained ForkML model (`.h5`).      |
| `--sample`    | Sample name (used for output filenames).           |
| `--output_dir`| Directory where outputs will be stored.            |

## 🛠 Configuration

Edit the file `params_local.config` to set defaults or cluster-specific parameters like CPU, memory, or environment module loading (e.g., `params.config`).

## 📦 Output

The pipeline outputs:

- ForkML prediction `.rds` files, forks, speeds, and events (in `analysis_<sample>/`):

| File Name                                      | Description                                                                              |
| ---------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `ForkML.basic_results_summary.<sample>.txt`    | Plain-text summary report with overall yield, detection, and fork speed stats.           |
| `ForkML.forks.<sample>.tsv`                    | Table of detected replication forks per read, including genomic positions and direction. |
| `ForkML.speeds.<sample>.tsv`                   | Per-fork speed estimates, with read and fork metadata.                                   |
| `ForkML.events.<sample>.tsv`                   | Detected initiation and termination events for each read.                                |
| `ForkML.speeds_distribution.<sample>.pdf`      | PDF plots showing fork speed distributions by read length bin.                           |
| `ForkML.typical_signals_examples.<sample>.pdf` | Example plots of raw signals and detected forks/events on representative reads.          |

- ForkML input dataset information, processed signal, and raw prediction files:

| File Name                              | Description                                                                                  |
| -------------------------------------- | -------------------------------------------------------------------------------------------- |
| `inputStats_<sample>.txt`              | Sequencing summary: total reads, mapped reads, yield, average read length, etc.              |
| `scores_<sample>.txt`                  | Signal quality scores for filtering reads (used for QC and legacy mode filtering).            |
| `signals_<sample>.tsv.gz`              | Compressed table of all binned signals for the sample (per-read, per-position).              |
| `tracks_<sample>.tsv.gz`               | Compressed table of detected tracks for the sample; raw detection results before summarizing. |

- Logs and trace files (if `-with-report`, `-with-trace`, etc. are used):

| File Name                              | Description                                                                                  |
| -------------------------------------- | -------------------------------------------------------------------------------------------- |
| `trace.txt`                            | Execution log for workflow steps (useful for debugging or reproducibility).                  |
| `timeline.html`                        | Workflow timeline or HTML report (if available; shows step execution and timing).            |
| `report.html`                          | Main analysis report in HTML format (open in browser).                                       |

## 📦 Process Descriptions from `forkml.nf`

### `SplitBam`
**Purpose**: Splits the input BAM file into smaller chunks to allow parallel processing.
- Extracts BAM header
- Divides reads into batches using `samtools view` and GNU `split`
- Balances CPU usage between viewing and filtering

### `ParseBam`
**Purpose**: Converts BAM chunks into fasta and ratio representations used for signal extraction.
- Outputs `.fa.gz`, `.fa_ratio.gz`, `.fa_cigar.gz` using `preprocess_modbam`

### `PreprocessFasta`
**Purpose**: Converts fasta-based Fork-seq data into structured `.rds` format.
- Filters reads by length
- Subsamples based on `subsampling_step`
- Output: single `.rds` file with cleaned and aligned reads

### `Segmentation`
**Purpose**: Uses machine learning to segment reads by replication direction.
- Inputs `.rds` files
- Outputs tracks and signal segments in `.tsv.gz`
- Relies on a defined trained model

### `Detection`
**Purpose**: Refines and fits tracks to fork models using detected segments.
- Combines segmentation and raw signal data
- Outputs updated `tracks_fit_*.tsv.gz` and `signals_fit_*.tsv.gz`

### `MergeTracks`
**Purpose**: Merges track tables from all read chunks into a single track file.
- Removes duplicate headers
- Produces a compressed merged `.tsv.gz` file

### `MergeSignals`
**Purpose**: Merges signal files from all chunks into one.
- Similar to `MergeTracks` but for signal data

### `AnalysisTracks`
**Purpose**: Analyzes the fork detection results to produce summary stats.
- Uses `forkml_speed_analysis`
- Outputs QC reports and possibly plots depending on analysis tool behavior

## Generating input dataset

### With `megalodon` for R9
Needs `Guppy` (v4.4.1) installed. `BrdU_megalodon.cfg` and `BrdU_megalodon.json` can be obtained from https://doi.org/10.5281/zenodo.15706384.
```bash
megalodon.forkml (){
  input_data=$1
  reference_genome=$2
  output_dir=$3
  nb_threads=$4 # minimum 8

  conda env create -f megalodon.yml -n megalodon
  conda activate megalodon

  mkdir -p $output_dir
  megalodon --guppy-server-path <path_to_guppy_install>/guppy/4.4.1/bin/guppy_basecall_server $input_data --guppy-config BrdU_megalodon.cfg --outputs basecalls mod_mappings --output-dir $output_dir --reference $reference_genome --processes $((nb_threads - 7)) --overwrite --device cuda:0 --guppy-timeout 900 --disable-mod-calibration --mod-min-prob 0
  
  samtools sort -@ $nb_threads -m 4G -o $output_dir/mod_mappings.sorted.bam $output_dir/mod_mappings.bam
  samtools index -@ $nb_threads $output_dir/mod_mappings.sorted.bam
  if [ -e "$output_dir/mod_mappings.sorted.bam" ]; then
    rm $output_dir/mod_mappings.bam
  fi

  pigz -p $nb_threads $output_dir/basecalls.fastq
  conda deactivate
}
```

### With `DNAscent` for R10
Needs `dorado` (v0.7.3) and `singularity` (v3.8.6)  installed, as well as `DNAscent` image (v4.0.3 with `singularity pull DNAscent.sif library://mboemo/dnascent/dnascent:4.0.3`).
```bash
dnascent.forkml() {
  # Usage check
  if [[ $# -ne 5 ]]; then
    echo "Usage:"
    echo "  dnascent.forkml <Pod5_dir> <Outdir> <Reference.fa> <naming_prefix> <cuda_gpu_number>"
    echo "Example:"
    echo "  dnascent.forkml /path/pod5 /path/outdir /path/ref.fa <sample_name> 0"
    return 1
  fi

  # Arguments
  local Pod5="$1"
  local Outdir="$2/bamout"
  local Index="$2/DS4_index"
  local DS4bam="$2/DS4_bam"
  local Ref="$3"
  local naming="$4"
  local gpu_num="$5"

  # Edit these as needed
  local main_path="/path/to/DNAscent4"
  local Dorado="${main_path}/dorado-0.7.3-linux-x64/bin/dorado"
  local Model="${main_path}/dorado_mdl/dna_r10.4.1_e8.2_400bps_fast@v5.0.0"
  local DS_sif="${main_path}/DNAscent.sif"

  mkdir -p "$Outdir" "$Index" "$DS4bam"

  echo "Step 1: Basecalling with Dorado"
  "$Dorado" basecaller -r -x cuda:"$gpu_num" "$Model" "$Pod5" --reference "$Ref" > "$Outdir/$naming.bam"

  echo "Step 2: Indexing with DNAscent"
  singularity run --bind "$Pod5":"$Pod5" --bind "$Index":"$Index" --nv "$DS_sif" \
    index -f "$Pod5/" -o "$Index/$naming.index"

  echo "Step 3: Detection with DNAscent"
  singularity run --bind "$Pod5":"$Pod5" --bind "$Index":"$Index" --bind "$Ref":"$Ref" --bind "$Outdir":"$Outdir" --bind "$DS4bam":"$DS4bam" --nv "$DS_sif" \
    detect -b "$Outdir/$naming.bam" -r "$Ref" -i "$Index/$naming.index" -o "$DS4bam/${naming}_DS4.bam" --GPU "$gpu_num" -t 20
}

```

