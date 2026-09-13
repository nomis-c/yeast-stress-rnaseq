# yeast-stress-rnaseq

## Introduction

This repository contains a Snakemake workflow to study transcriptional response of *Saccharomyces cerevisiae* to heat stress.
The workflow has been tested on Oracle Linux 9 using the Life Science Compute Cluster at the University of Vienna.

## Pipeline

The workflow starts from raw sequencing reads and ends with a list of genes
that change their activity under heat stress, together with the biological
processes those genes belong to. It consists of 4 key tools.

| Tool | Step |
|------|------|
| `fastp` | cleans the raw reads (adapters, low-quality sequence) |
| `Salmon` | measures how active each gene is, by matching reads against the yeast reference transcriptome |
| `DESeq2` | compares heat-stressed to control samples and finds the genes that changed significantly |
| `clusterProfiler` | asks which biological functions are over-represented among those genes |

The pipeline is not tied to this dataset: to run it on other samples, edit
`resources/samples.tsv` (which samples to download) and `config/config.yaml`
(which groups to compare).

## Usage

Users should first install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), then clone the repository and activate the environment:

```
git clone https://github.com/nomis-c/yeast-stress-rnaseq.git
cd yeast-stress-rnaseq
mamba env create -f workflow/envs/env.yaml
mamba activate ysr
```

Run the analysis locally:
```
snakemake -c 1
```

Run the analysis on HPC:

```
snakemake -c 1 --profile config/slurm
```

## Results

The results match known observations about heat-stressed yeast. Cells activate
heat-shock genes like chaperones, which protect other proteins from being
damaged. At the same time, they switch off genes needed for growth and for
building new ribosomes.

For a detailed report, please read [report.md](report.md).
