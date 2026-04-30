# TAM Screening Analysis

Reusable analysis code for TAM/PAM screening data.

This repository keeps scripts, notebook templates, configuration, and documentation under version control. Raw sequencing data, downloaded archives, extracted FASTQ files, pickle count files, and full result folders should stay outside Git.

## Current Fixed Workflow

The fixed workflow compares each experiment against a control N5 library and runs six analysis configurations:

- `raw`: direct `Y / X`
- `corrected`: `Y` divided by `median(Y / X)`
- thresholds: `1.5sigma`, `5sigma`, and hard `FC=0.15`
- both depleted and enriched PAMs are reported
- each plot contains scatter, FC distribution, selected motif, and TopN motif

No low-abundance filter is applied in the fixed workflow.

## Quick Start

Create an environment:

```bash
conda env create -f environment.yml
conda activate tam-screening-analysis
```

Prepare a sample sheet:

```bash
cp configs/sample_sheet_template.csv sample_sheet.csv
```

Extract 5N counts:

```bash
python scripts/extract_5n_counts.py   --sample-sheet sample_sheet.csv   --out-dir processed_data
```

Run fixed six-configuration analysis:

```bash
python scripts/run_fixed6_analysis.py   --counts-dir processed_data   --control N5   --samples 118p2_200 118p3_200 347p0_200 347p1_200 CasX   --out-dir results/fixed6
```

Compare candidate N5 libraries:

```bash
python scripts/compare_n5_libraries.py   --n5 N5_20260312=/path/to/N5_counts.p   --n5 N5_20260401=/path/to/N5_counts.p   --n5 N5_20260429=/path/to/N5_counts.p   --out-dir results/n5_qc
```

## What Should Not Be Committed

Do not commit:

- raw FASTQ files
- downloaded vendor archives
- OSS URLs containing signatures or access keys
- large generated image folders
- full per-batch result directories
- local machine-specific paths

Commit:

- scripts
- templates
- configuration examples
- small example tables if needed
- documentation of analysis choices

## Recommended GitHub Flow

```bash
git init
git add .
git commit -m "Initialize TAM screening analysis workflow"
git branch -M main
git remote add origin git@github.com:<your-user-or-org>/TAM-screening-analysis.git
git push -u origin main
```
