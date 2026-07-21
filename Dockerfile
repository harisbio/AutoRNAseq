FROM rocker/shiny:4.6.0

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxt-dev \
    libglpk-dev \
    libgmp-dev \
    libcairo2-dev \
    libpq-dev \
    samtools \
    fastqc \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('remotes','BiocManager'), repos='https://cloud.r-project.org')"

COPY DESCRIPTION /tmp/AutoRNAseq/DESCRIPTION

RUN R -e "options(repos = c(CRAN='https://cloud.r-project.org')); bioc <- c('Rqc','QuasR','Rsubread','DESeq2','EnhancedVolcano','apeglm','pheatmap','rintrojs','clusterProfiler','enrichplot','DOSE','AnnotationDbi','org.Hs.eg.db','org.Mm.eg.db','pathview','KEGGREST','fgsea','msigdbr','STRINGdb','igraph','visNetwork','ggraph','tidygraph','caret','randomForest','e1071','glmnet','xgboost','pROC','plotly','rmarkdown','BiocParallel','rtracklayer','Biostrings'); BiocManager::install(bioc, ask = FALSE, update = FALSE)"

COPY . /srv/AutoRNAseq

RUN R -e "remotes::install_local('/srv/AutoRNAseq', dependencies = FALSE, upgrade = 'never')"

EXPOSE 3838

CMD ["R", "-e", "AutoRNAseq::runAnalyser(host='0.0.0.0', port=3838)"]
