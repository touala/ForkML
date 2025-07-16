# ForkML Training Pipeline

A **Nextflow** pipeline to preprocess and train models on DNA replication fork signal data using double pulse annotation datasets.

---

## 📦 Requirements

- [Nextflow ≥ 22.10.7](https://www.nextflow.io/)
- Conda (e.g., via [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
- Bash, GNU coreutils, and a Linux environment

### 🧪 Environment Setup

```bash
conda env create -p forkml_training -f forkml_training.yml
conda activate forkml_training

bash post_install.sh # (Optional) Infrastructure dependent
```

---

## 🚀 Usage

### Run the pipeline

```bash
nextflow run path/to/training.nf \
  -c path/to/training_local.config \
  -with-report output_dir/report.html \
  -with-trace output_dir/trace.txt \
  -with-timeline output_dir/timeline.html \
  --input_bams "/path/to/annotated.sorted.bam" \
  --annotations "/path/to/annotation/sample_name.batch_*.rds" \
  --output output_dir \
  --nb_reads_process 200 \
  --max_read_length 200000 \
  --with_augmentation \
  --padding_method clean \
  --data_type DNAscent \
  --with_training \
  --tmpdir "/localtmp/forkml" \
  -resume
```

---

## ⚙️ Parameters

| Parameter             | Description                                                | Default        |
|----------------------|------------------------------------------------------------|----------------|
| `--input_bams`       | Comma-separated BAM files to process                       | **Required**   |
| `--annotations`      | Glob path or list of RDS annotation files                  | **Required**   |
| `--output`           | Output directory                                            | **Required**   |
| `--nb_reads_process` | Number of reads per BAM chunk                              | `2000`         |
| `--max_read_length`  | Maximum length to retain a read (bp)                       | `200000`       |
| `--padding_method`   | Type of padding to apply (`clean`, `rough`, `default`)     | `clean`        |
| `--with_augmentation`| Enable synthetic data augmentation                         | `false`        |
| `--with_training`    | Train a model at the end of preprocessing                  | `false`        |
| `--tmpdir`           | Directory for temporary files                              | `$output_dir`  |

---

## 📁 Outputs

- Trained ForkML model `.h5`
- Preprocessed `*.rds` and `*.bam` files
- Balanced and unbalanced `.npy` arrays for training and validation
- Training logs, model files (if `--with_training`)
- Nextflow execution reports:
  - `trace.txt` – task timeline and resource usage
  - `timeline.html` – interactive task breakdown
  - `report.html` – process summary

---

## 💡 Tips

- Use `-resume` to skip completed steps and reuse cached results

---

## 📚 Citation

Please cite **[ForkML](https://github.com/touala/ForkML)** if you use this pipeline in your research.

---
