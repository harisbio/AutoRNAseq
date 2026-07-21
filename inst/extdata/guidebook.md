<div class="autornaseq-doc">
  <section class="doc-hero">
    <span class="doc-chip">Operational Manual</span>
    <h1>AutoRNAseq Guidebook</h1>
    <p>This guide walks you through the recommended analysis sequence from raw reads to biomarker prioritization. If you are new to RNA-Seq, follow the modules in order for the cleanest experience.</p>
  </section>

  <section class="doc-block">
    <h2>Start Here</h2>
    <ol>
      <li>Go to <strong>Data Setup</strong> and upload your FASTQ files plus the sample metadata CSV.</li>
      <li>Review <strong>Quality Control</strong> outputs before running trimming or alignment.</li>
      <li>Continue through the modules one by one, saving outputs as you go.</li>
    </ol>
  </section>

  <section class="doc-block">
    <h2>Workflow Blueprint</h2>
    <div class="doc-grid">
      <div><h3>1. Data Setup</h3><p>Load FASTQ/FASTQ.GZ files, upload metadata CSV, validate sample mapping, and define the grouping column.</p></div>
      <div><h3>2. Quality Control</h3><p>Profile read quality with <strong>Rqc</strong> and decide whether trimming is needed.</p></div>
      <div><h3>3. Filtering &amp; Trimming</h3><p>Apply read cleanup with <strong>QuasR</strong> and related preprocessing tools.</p></div>
      <div><h3>4. Alignment</h3><p>Map reads with <strong>Rsubread</strong> and inspect alignment behavior.</p></div>
      <div><h3>5. Quantification</h3><p>Generate gene-level counts from aligned reads and annotation files.</p></div>
      <div><h3>6. Differential Expression</h3><p>Use <strong>DESeq2</strong> and <strong>apeglm</strong> for robust inference and shrinkage.</p></div>
      <div><h3>7. Functional Analysis</h3><p>Run GO/KEGG/GSEA using <strong>clusterProfiler</strong>, <strong>pathview</strong>, <strong>fgsea</strong>, and <strong>msigdbr</strong>.</p></div>
      <div><h3>8. PPI &amp; ML</h3><p>Build interaction networks and rank candidate biomarkers with ML methods.</p></div>
    </div>
  </section>

  <section class="doc-block">
    <h2>Data Requirements</h2>
    <ul>
      <li>Input reads: FASTQ or FASTQ.GZ.</li>
      <li>Metadata: CSV with one row per sample and matching sample identifiers.</li>
      <li>Large projects may take time and memory, especially during alignment and quantification.</li>
      <li>For public deployment, keep uploads within the configured server limit and recommend compressed input when possible.</li>
    </ul>
  </section>

  <section class="doc-block">
    <h2>Execution Standards</h2>
    <ul>
      <li>Use consistent sample names across metadata and file names.</li>
      <li>Match annotation files to the reference genome build.</li>
      <li>Use adjusted p-values and effect sizes when interpreting differential expression.</li>
      <li>Interpret ML biomarkers together with DEG, pathway, and network evidence.</li>
    </ul>
  </section>
</div>
