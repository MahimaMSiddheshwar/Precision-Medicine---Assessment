# SARS-CoV-2 RNA-Seq Analysis

## Project Title

**Transcriptomic Analysis of SARS-CoV-2–Infected Human Respiratory Cells Using RNA-Seq**

---

## Overview

This project implements an end-to-end **RNA-Seq analysis pipeline** to study host transcriptomic responses in **human respiratory cells infected with SARS-CoV-2**. The workflow spans raw data acquisition to downstream **differential gene expression (DEG)** and **Gene Ontology (GO) enrichment analysis**, comparing control and infected samples at **24H and 72H time points**.

---

## Objectives

* Process raw RNA-Seq data from SRA to gene-level count matrices
* Identify differentially expressed genes (DEGs)
* Perform GO Biological Process enrichment analysis
* Compare time-dependent (24H vs 72H) and infection-specific (Control vs Infected) responses

---

## Dataset

* **Source:** NCBI Sequence Read Archive (SRA)
* **Samples:** Control vs SARS-CoV-2 infected respiratory cells
* **Time Points:** 24H and 72H
* **Reference Genome:** UCSC hg38 (genome + GTF)

---

## Workflow Summary

### 1. Data Acquisition & Quality Control

* Downloaded SRA files using **sra-toolkit**
* Converted SRA to FASTQ using **fasterq-dump**
* Performed initial QC using **FastQC** and aggregated reports with **MultiQC**

### 2. Read Preprocessing

* Adapter and quality trimming using **Trim Galore**
* Post-trimming QC validation using FastQC

### 3. Reference Genome Preparation

* Downloaded UCSC hg38 reference genome and annotation files
* Generated STAR genome index

### 4. Read Alignment

* Aligned trimmed reads to hg38 using **STAR**
* Generated coordinate-sorted BAM files
* Indexed BAM files using **samtools**
* Filtered **primary alignments only** to ensure accurate read quantification

### 5. Gene Expression Quantification

* Performed gene-level quantification using **featureCounts**
* Generated count matrices for all BAM files and primary-alignment BAM files

### 6. Differential Expression Analysis (R)

* Normalized count data using **edgeR**
* Modeled expression using **limma** (GLM framework)
* Differential comparisons performed:

  * 24H vs 72H
  * Control vs SARS-CoV-2 infected
* Visualized results using **EnhancedVolcano**

### 7. GO Enrichment Analysis

* Filtered DEGs using |logFC| > 1 and p-value < 0.05
* Converted gene symbols to Entrez IDs
* Performed GO Biological Process enrichment using **clusterProfiler**
* Visualized enriched pathways using barplots

---

## Key Results

* **Control vs SARS-CoV-2 infected** comparison revealed extensive transcriptional reprogramming
* Significant upregulation of **snoRNAs (SNORD116 family)** observed
* Enrichment of immune response, ribonucleoprotein biogenesis, and cell-cycle–related pathways
* Time-dependent changes (24H vs 72H) were present but comparatively subtle

---

## Biological Interpretation

The upregulation of **SNORD family genes** aligns with prior observations in lung diseases such as COPD, suggesting a potential role in **viral replication, host immune modulation, and inflammatory signaling**. These findings indicate that dysregulated snoRNA expression may serve as **novel biomarkers or therapeutic targets** in SARS-CoV-2–associated respiratory pathology.

---

## Outputs

* Volcano plots (PNG)
* GO enrichment result tables (CSV)
* FeatureCounts gene expression matrices
* Alignment statistics and BAM index files

---

## Tech Stack & Skills

**RNA-Seq | NGS | Linux | HPC | SLURM | Bash | R | STAR | FastQC | MultiQC | Trim Galore | samtools | featureCounts | edgeR | limma | clusterProfiler | EnhancedVolcano**

### Tool Links

* STAR: [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
* FastQC: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* MultiQC: [https://multiqc.info](https://multiqc.info)
* Trim Galore: [https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* Subread / featureCounts: [https://subread.sourceforge.net](https://subread.sourceforge.net)
* edgeR: [https://bioconductor.org/packages/edgeR](https://bioconductor.org/packages/edgeR)
* limma: [https://bioconductor.org/packages/limma](https://bioconductor.org/packages/limma)
* clusterProfiler: [https://bioconductor.org/packages/clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
* EnhancedVolcano: [https://bioconductor.org/packages/EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano)

---

