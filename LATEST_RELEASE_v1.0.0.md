# AutoRNAseq Latest Release v1.0.0

This file marks the current public release of AutoRNAseq.

## Release version

`v1.0.0`

## Release summary

AutoRNAseq v1.0.0 is the first public release of the project. It provides a Dockerized, browser-based RNA-seq analysis workflow for local use, along with direct R package usage for advanced users who want to run the application from R.

## Included capabilities

- Data setup and study metadata handling
- Quality control
- Filtering and trimming
- Alignment
- Quantification
- Differential expression analysis
- GO enrichment
- KEGG pathway analysis
- GSEA
- PPI network analysis
- Biomarker discovery with machine learning
- Automated pipeline execution

## Release files

The release package should include:

- `docker-compose.yml`
- `Dockerfile`
- `start-autornaseq-windows.bat`
- `start-autornaseq-mac-linux.sh`
- `README_R.md`
- `RELEASE_NOTES_v1.0.0.md`
- the source code package
- the GitHub Pages landing page in `docs/`

## User install paths

### Docker

1. Install Docker Desktop.
2. Run the startup script or `docker compose up -d --build shinyapp`.
3. Open `http://localhost:3838`.

### R

1. Install the package and its dependencies.
2. Run `library(AutoRNAseq)`.
3. Launch the app with `runAnalyser()`.

## Links

- Landing page: https://harisbio.github.io/AutoRNAseq/
- Source repository: https://github.com/harisbio/AutoRNAseq
- Release notes: `RELEASE_NOTES_v1.0.0.md`

