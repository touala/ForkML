# ForkML: Make figures

This directory provides a script for making most of the figures.

## 🚀 Quick Start

```bash
conda env create -f forkml_figure.yml -n forkml_figure
conda activate forkml_figure

wget -O datasets_figures.tar.gz https://zenodo.org/records/15706385/files/datasets_figures.tar.gz
tar -xzf datasets_figures.tar.gz
Rscript --vanilla forkml_make_figures.R ./datasets_figures # figures made in local directory

```
