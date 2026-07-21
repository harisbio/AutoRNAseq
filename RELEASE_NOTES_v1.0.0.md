# AutoRNAseq v1.0.0

AutoRNAseq v1.0.0 is the first public release of the project.

AutoRNAseq is a Dockerized, browser-based RNA-seq analysis platform for local, reproducible transcriptomics workflows. It is designed to help users move from raw FASTQ files to downstream biological interpretation without needing to write code or maintain a hosted server.

## Highlights

- Local Docker-based deployment for private, reproducible analysis
- Browser interface running at `http://localhost:3838`
- Guided workflow for non-coders and researchers
- Quality control, trimming, alignment, quantification, and DESeq2 differential expression
- Downstream GO, KEGG, and GSEA enrichment analysis
- STRINGdb-supported protein-protein interaction network analysis
- Machine-learning-assisted biomarker ranking
- Downloadable plots, tables, and report outputs
- Built-in help text and workflow guidance across modules

## Main modules

This release includes the following major modules:

1. Data Setup
2. Quality Control
3. Filtering and Trimming
4. Alignment
5. Quantification
6. Differential Expression
7. GO Enrichment
8. KEGG Pathway Analysis
9. GSEA
10. PPI Network Analysis
11. Biomarker Discovery with Machine Learning
12. Automated Pipeline

## Installation

1. Install Docker Desktop.
2. Download the release assets.
3. Run the included startup script for your operating system.
4. Open `http://localhost:3838` in your browser.

## Release assets

Recommended files to attach to this release:

- `docker-compose.yml`
- Windows startup script
- macOS/Linux startup script
- quick-start guide
- example metadata or test dataset, if included

## Notes

- This release is intended for local execution, not hosted server deployment.
- The recommended workflow is to keep user data on the local machine and run analyses inside Docker.
- Users should confirm that their input metadata, FASTQ files, reference genome, and annotation files are consistent before starting the automated pipeline.

## Citation

If you use AutoRNAseq in a publication, please cite the software and the underlying tools used in the analysis workflow.

