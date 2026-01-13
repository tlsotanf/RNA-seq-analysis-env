# RNA-seq Analysis Environment

## Overview
This repository documents a **reproducible RNA-seq analysis environment** built on  
**WSL2 (Ubuntu)**, designed to bridge **Linux-based RNA-seq pipelines** and  
**R-based downstream analysis**.

The goal of this project is to establish a clean and reproducible setup for  
RNA-seq analysis on a Windows machine while maintaining compatibility with  
widely used **Linux-only bioinformatics tools**.

---

## Motivation
Many core RNA-seq tools (e.g. **STAR**, **samtools**, **featureCounts**, **Nextflow**)  
are developed primarily for Linux environments and are difficult or unstable to  
run natively on Windows.

At the same time, **R-based downstream analyses**  
(DEG analysis, enrichment, visualization) are often more convenient in a  
Windows-based **RStudio workflow**.

This project adopts the following strategy:

### Analysis Strategy
**Linux (WSL2 / Ubuntu)**  
→ Alignment, quantification, workflow orchestration  

**Windows / RStudio**  
→ Differential expression analysis, functional enrichment, visualization  

This separation improves:
- Reproducibility
- Tool compatibility
- Developer productivity

---

## Environment Architecture
```text
Windows 11
│
├── RStudio (Windows)
│   ├── DESeq2 / edgeR
│   ├── clusterProfiler
│   └── Visualization & reporting
│
└── WSL2 (Ubuntu 24.04)
    ├── Conda environments
    ├── RNA-seq tools (STAR, samtools, etc.)
    └── Workflow engine (Nextflow)

## Repository Scope
This repository focuses on:
- Documenting environment setup
- Managing configuration and scripts
- Ensuring reproducibility across machines

It does **not** store:
- Large binary files
- Installer packages
- Raw sequencing data

---

## Version Control & Security Practices
- GitHub **fine-grained, repository-scoped Personal Access Tokens**
- Explicit permission control (**Contents: Read & Write**)
- `.gitignore` configured to exclude:
  - Installer files
  - Local notes
  - Temporary artifacts

  ## Author
**Geonwoo**  
Background in bioengineering and RNA-seq data analysis  
Interested in reproducible bioinformatics pipelines and scalable analysis workflows

---

## License
This repository is intended for **educational and portfolio purposes**.  
Reuse or modification is allowed with proper attribution.