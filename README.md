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

## Pipeline Design Statement
This RNA-seq pipeline is designed as a modular and reproducible workflow
orchestrated by Nextflow, with explicit separation between
data processing (alignment and quantification) and downstream statistical analysis.

The pipeline prioritizes robustness, transparency, and reproducibility,
while remaining flexible for future extensions such as alternative quantification
methods or additional quality control steps.

## Design Rationale

This pipeline was designed by referencing the overall structure and best practices
of the **nf-core/rnaseq** workflow, which serves as a widely adopted community standard
for reproducible RNA-seq analysis.

Rather than reusing the nf-core pipeline directly, this project focuses on
independently implementing a simplified and transparent RNA-seq workflow
to deepen understanding of workflow design, tool integration, and data flow
management using Nextflow.

Key design decisions include:

- **Workflow orchestration with Nextflow**  
  to ensure scalability, parallel execution, and reproducibility across systems.

- **Conda-based dependency management**  
  to provide isolated, reproducible environments for each analysis step without
  manual tool installation.

- **Alignment-based quantification (STAR + featureCounts)**  
  as the default strategy, prioritizing interpretability, robustness, and
  compatibility with widely used differential expression tools such as
  DESeq2 and edgeR.

- **Clear separation between upstream processing and downstream analysis**,  
  where computationally intensive steps are executed in a Linux environment (WSL2),
  while statistical analysis and visualization are performed in RStudio on Windows.

The nf-core/rnaseq project is acknowledged as a conceptual reference for workflow
organization and reproducibility principles, while the implementation, tool
selection, and parameterization in this repository are intentionally customized
for learning and portfolio purposes.

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