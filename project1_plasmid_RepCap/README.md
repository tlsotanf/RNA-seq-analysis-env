# Project 1: Plasmid RepCap Expression Analysis

## Overview
RNA-seq analysis of RepCap plasmid expression in AAV production.

## Experimental Design
- **Cell lines**: 3 (cell_line_1, cell_line_2, cell_line_3)
- **Time points**: 0h, 24h, 72h
- **Replicates**: 2-3 biological replicates

## Pipeline
- Nextflow workflow
- Salmon quantification
- DEG analysis with t-test and FDR correction

## Results
- Differential expression: 0h vs 24h, 0h vs 72h
- Significant upregulation at 24h and 72h in all cell lines

## Files
- `main.nf`: Nextflow pipeline
- `nextflow.config`: Pipeline configuration
- `reference/`: Plasmid sequences (FASTA, GTF)
- `results/`: MultiQC reports
- `docs/`: Analysis documentation

## Usage
```bash
nextflow run main.nf -c nextflow.config
```
