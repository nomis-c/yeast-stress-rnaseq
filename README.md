# yeast-stress-rnaseq

## Introduction

This repository contains a Snakemake workflow to study transcriptional response of *Saccharomyces cerevisiae* to heat stress.
The workflow has been tested on Oracle Linux 9 using the Life Science Compute Cluster at the University of Vienna.

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

For a detailed report, please read [report.md](report.md).
