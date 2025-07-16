
# ForkML Annotating Pipeline

This Nextflow pipeline prepares datasets for manual annotation of DNA replication tracks using the ForkML toolkit.

---

## 🧩 Overview
This pipeline extracts signal-rich reads from a BAM file and packages them into datasets readable with an annotation app. Reviewers can then manually annotate these reads using a Shiny-based interface.

---

## 🔧 Usage
```bash
nextflow run annotating.nf \
  -c annotating/params_local.config \
  -resume \
  -with-report output_dir/report.html \
  -with-trace output_dir/trace.txt \
  -with-timeline output_dir/timeline.html \
  --bam /path/to/input.bam \
  --output_dir output_dir \
  --sample sample_name \
  --res_step 100 \
  --thre_level 0.1 \
  --num_points 10 \
  --size_thre 1 \
  --threads_select 2 \
  --threads_prepare 8 \
  --reviewers "Reviewer1,Reviewer2" \
  --nb_random 0 \
  --nb_with_signal 1000 \
  --min_read_length 10000 \
  --batch_size 50 \
  --with_independent_review true \
  --target_mod_type "N+b?"
```

---

## 🧪 Parameters

### Required
- `--bam` : Path to the input BAM file.
- `--output_dir` : Output directory for annotated data.
- `--sample` : Sample name used as basename.

### SelectReadsWithSignal
- `--res_step` : Signal resolution (e.g., 100 for 1pt/100bp).
- `--thre_level` : Threshold for signal presence.
- `--num_points` : Minimum number of points with signal.
- `--size_thre` : Minimum segment size (kb).
- `--threads_select` : Number of threads for signal filtering.

### PrepareAnnotation
- `--reviewers` : Comma-separated names of reviewers.
- `--nb_random` : Number of random reads to include.
- `--nb_with_signal` : Number of signal-rich reads to include.
- `--min_read_length` : Minimum read length.
- `--batch_size` : Number of reads per annotation batch.
- `--threads_prepare` : Number of threads for preparation.
- `--with_independent_review` : Whether to include both reviewers independently.
- `--target_mod_type` : Modification string used to filter reads.

---

## 🧬 Environment
Create and activate the environment with:
```bash
conda env create -f forkml_annotating.yml -p forkml_annotating
conda activate forkml_annotating
```

---

## 📦 Outputs
The pipeline generates:
- Filtered BAM + index.
- Annotation-ready `.rds` datasets.
- Shiny app directory for each reviewer in `output_dir/forkml_annotation/for_<reviewer>/`.

Each review folder contains:
- `samples.csv` : Sample metadata file.
- `reads_toannotate/` : Data to annotate.
- `app.R` : Shiny app for annotation.

Launch app locally:
```bash
Rscript app.R
# Then open http://127.0.0.1:9999/ in your browser.
```

---

## 📈 Reports
Outputs also include:
- `report.html` : Summary of execution.
- `trace.txt` : Per-process runtime info.
- `timeline.html` : Execution timeline.

---

## 📁 Example
```bash
nextflow run annotating.nf \
  -c params_local.config \
  --bam example.bam \
  --output_dir results \
  --sample EX123 \
  --reviewers "Alice,Bob" \
  --nb_with_signal 500 \
  -resume
```

---

## 📫 Contact
Issues or questions? Open an issue on the [GitHub repository](https://github.com/touala/forkml).
