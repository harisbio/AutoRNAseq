# Using AutoRNAseq in R

AutoRNAseq can be used directly as an R package if you already have R installed and are comfortable managing package dependencies.

## Install from a local folder

```r
install.packages("remotes")
remotes::install_local("D:/dr_murad/antigravity_MS_project/SeqExpressionAnalyser-master/AutoRNAseq")
library(AutoRNAseq)
runAnalyser()
```

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("AutoRNAseq/AutoRNAseq")
library(AutoRNAseq)
runAnalyser()
```

## Windows notes

On Windows, install `Rtools` if R asks for it, and use `BiocManager` to install the required Bioconductor packages if they are missing.

```r
install.packages("BiocManager")
BiocManager::install(c(
  "Rqc","QuasR","Rsubread","DESeq2","clusterProfiler","enrichplot",
  "DOSE","org.Hs.eg.db","org.Mm.eg.db","pathview","fgsea","STRINGdb","rtracklayer"
))
```

## Docker notes

If you prefer not to manage dependencies yourself, use the Docker route instead:

```powershell
docker compose up -d --build shinyapp
```

Then open:

```text
http://localhost:3838
```
