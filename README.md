# AutoRNAseq

[Open the landing page](https://harisbio.github.io/AutoRNAseq/)

AutoRNAseq is a Dockerized, browser-based RNA-seq analysis platform for end-to-end transcriptomics workflows and machine-learning-assisted biomarker discovery. It is designed for researchers who want a guided local application instead of a hosted web service, so sequencing data stays on their own machine or institutional system.

The app runs locally in a browser at `http://localhost:3838` after Docker starts. This makes it practical for non-coders, keeps analysis private, and gives every user the same reproducible software environment.

## Installation at a glance

Choose one of the two supported paths:

1. Docker, recommended for most users
2. R package installation, for advanced users who already manage R dependencies

## 1. Docker install

### Prerequisites

- Docker Desktop
- Internet access for the first build
- Enough disk space for your FASTQ files, reference files, and outputs

### Start AutoRNAseq

From the project folder:

```powershell
docker compose up -d --build shinyapp
```

Or use the provided launchers:

```text
start-autornaseq-windows.bat
start-autornaseq-mac-linux.sh
```

Then open:

```text
http://localhost:3838
```

### Stop AutoRNAseq

```powershell
docker compose down
```

## 2. R install

If you want to run AutoRNAseq directly inside R:

```r
install.packages("remotes")
remotes::install_local("D:/dr_murad/antigravity_MS_project/SeqExpressionAnalyser-master/AutoRNAseq")
library(AutoRNAseq)
runAnalyser()
```

If dependencies are missing on Windows, install `Rtools` first and then use `BiocManager` for the Bioconductor packages.

```r
install.packages("BiocManager")
BiocManager::install(c(
  "Rqc","QuasR","Rsubread","DESeq2","clusterProfiler","enrichplot",
  "DOSE","org.Hs.eg.db","org.Mm.eg.db","pathview","fgsea","STRINGdb","rtracklayer"
))
```

## What AutoRNAseq does

AutoRNAseq brings the common RNA-seq workflow into one Shiny interface:

1. Import study metadata and raw FASTQ inputs.
2. Run quality control with Rqc, FastQC, and MultiQC-based workflows.
3. Trim or filter reads with QuasR, fastp, or Trimmomatic.
4. Align reads with supported RNA-seq alignment tools.
5. Generate gene-level count matrices with featureCounts/Rsubread workflows.
6. Run DESeq2-based differential expression analysis.
7. Interpret results with GO, KEGG, and GSEA modules.
8. Explore protein-protein interactions with STRINGdb.
9. Rank candidate biomarkers with machine-learning models.
10. Produce downloadable tables, plots, and report-ready outputs.

## What is included in the release

- `docker-compose.yml`
- `Dockerfile`
- `start-autornaseq-windows.bat`
- `start-autornaseq-mac-linux.sh`
- `README_R.md`
- source package code
- landing page under `docs/`

## Why this distribution model

RNA-seq analysis is resource-heavy and often awkward to deploy for non-technical users. AutoRNAseq uses a local Docker workflow instead of a hosted server so it can:

- avoid central hosting costs;
- keep raw sequencing data private;
- reduce installation friction;
- provide a reproducible runtime;
- work on laptops, workstations, and institutional machines.

## Quick start

### Docker users

1. Install and start Docker Desktop.
2. Run `docker compose up -d --build shinyapp`.
3. Open `http://localhost:3838`.
4. Upload your data and run the workflow.

### R users

1. Install R and the required dependencies.
2. Install AutoRNAseq from the local folder or GitHub.
3. Run `library(AutoRNAseq)`.
4. Launch the app with `runAnalyser()`.

### Run your analysis

Upload your FASTQ files and metadata, then move through the sidebar modules in order:

- Data Setup
- Quality Control
- Filtering and Trimming
- Alignment
- Quantification
- Differential Expression
- GO / KEGG / GSEA
- PPI Network Analysis
- Biomarker Discovery

## Public download and release flow

The easiest public distribution pattern is:

1. Publish the source code on GitHub.
2. Create GitHub Releases for versioned downloads.
3. Add a GitHub Pages landing page for non-coders.
4. Keep the app running locally through Docker.

The landing page can link to the latest release and give a short step-by-step quick start. A ready-to-publish landing page lives in [`docs/`](docs/).

## Installation options

### Recommended: Docker Desktop

This is the best option for most users. It bundles the app with the required R, Bioconductor, and command-line dependencies in one consistent runtime.

### Advanced: R package usage

If you already manage R dependencies yourself, you can launch the app directly:

```r
library(AutoRNAseq)
runAnalyser()
```

## Main features

- Guided browser interface for RNA-seq analysis
- Local-first deployment with no hosted server requirement
- Quality control, trimming, alignment, and quantification modules
- Differential expression analysis with DESeq2
- Functional enrichment with GO, KEGG, and GSEA
- STRINGdb-supported PPI network analysis
- Supervised biomarker ranking with caret-based machine learning
- Downloadable plots, tables, and report outputs
- Built-in help text and module guidance

## Example workflow

1. Put FASTQ files in the mounted user data folder or provide SRA accessions.
2. Upload metadata and choose the condition/sample columns.
3. Run quality control and trimming if needed.
4. Align reads and generate BAM files.
5. Quantify counts with featureCounts.
6. Run DESeq2 and inspect the volcano plot.
7. Interpret the DEG list with GO, KEGG, GSEA, and PPI modules.
8. Rank biomarkers with the ML module.

## Example data

The package includes example FASTQ files and sample metadata under `inst/extdata/` for demonstration and interface testing.

## Repository structure

- `R/` - Shiny modules and package functions
- `inst/extdata/` - help text, guidebook content, example inputs
- `inst/www/` - UI assets, stylesheet, workflow image, and branding
- `inst/report_templates/` - R Markdown report template
- `docs/` - public GitHub Pages landing page
- `manuscript/` - manuscript draft and submission support files

## License

AutoRNAseq is licensed under the GNU General Public License v3.0 or later. See [LICENSE.md](LICENSE.md) for the full license text.

## Citation

Under publication yet.

## Contact

**Muhammad Haris**  
Department of Bioinformatics and Biotechnology  
Government College University Faisalabad, Pakistan  
Email: [mharis.202101862@gcuf.edu.pk](mailto:mharis.202101862@gcuf.edu.pk)  
Portfolio: [harisbio.github.io](https://harisbio.github.io/)

## Links

- [Project landing page](https://harisbio.github.io/AutoRNAseq/)
- [Latest release](https://github.com/harisbio/AutoRNAseq/releases/tag/v1.0.0)
- [Author portfolio](https://harisbio.github.io/)
