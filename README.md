## **Transcriptomic Analysis of SARS-CoV-2–Infected Human Respiratory Cells Using RNA-Seq**

This project implements an end-to-end **RNA-Seq analysis pipeline** to investigate host transcriptomic responses in **human respiratory cells infected with SARS-CoV-2**. The workflow spans raw data acquisition through downstream **differential gene expression (DEG)** and **Gene Ontology (GO) enrichment analysis**, comparing control and infected samples across **24H and 72H time points**.

---

## Objectives

* Process raw RNA-Seq data from SRA into gene-level count matrices
* Identify differentially expressed genes (DEGs)
* Perform Gene Ontology (GO) Biological Process enrichment analysis
* Compare time-dependent (24H vs 72H) and infection-specific (Control vs Infected) transcriptional responses

---

## Dataset

* **Source:** NCBI Sequence Read Archive (SRA)
* **Samples:** Control and SARS-CoV-2–infected human respiratory cells
* **Time Points:** 24H and 72H
* **Reference Genome:** UCSC hg38 (genome and GTF annotation)

---

## Workflow Summary

### 1. Data Acquisition & Quality Control

* Downloaded RNA-Seq data using **sra-toolkit**
* Converted SRA files to FASTQ format using **fasterq-dump**
* Performed initial quality control using **FastQC** and aggregated reports with **MultiQC**

### 2. Read Preprocessing

* Adapter and quality trimming using **Trim Galore**
* Post-trimming quality assessment using FastQC

### 3. Reference Genome Preparation

* Downloaded UCSC hg38 reference genome and gene annotation files
* Generated STAR genome indices

### 4. Read Alignment

* Aligned trimmed reads to hg38 using **STAR**
* Generated coordinate-sorted BAM files
* Indexed BAM files using **samtools**
* Filtered **primary alignments only** to ensure accurate and unbiased read quantification

### 5. Gene Expression Quantification

* Quantified gene-level expression using **featureCounts**
* Generated count matrices for both all BAM files and primary-alignment-only BAM files

### 6. Differential Expression Analysis (R)

* Normalized count data using **edgeR**
* Modeled expression using **limma** within a generalized linear modeling framework
* Differential comparisons performed:

  * **24H vs 72H**
  * **Control vs SARS-CoV-2 infected**
* Visualized DEGs using **EnhancedVolcano**

### 7. GO Enrichment Analysis

* Filtered DEGs using |logFC| > 1 and p-value < 0.05
* Converted gene symbols to Entrez IDs
* Performed GO Biological Process enrichment using **clusterProfiler**
* Visualized enriched pathways using bar plots

---

## Key Results

* The **Control vs SARS-CoV-2 infected** comparison revealed extensive transcriptional reprogramming
* Significant upregulation of **snoRNAs (SNORD116 family)** was observed
* Enrichment of immune response, ribonucleoprotein biogenesis, and cell-cycle–related pathways
* Time-dependent changes (24H vs 72H) were present but comparatively subtle

---

## Outputs

* `volcano_24h_vs_72h.png` – Volcano plot comparing gene expression between 24H and 72H
* `volcano_control_vs_infected.png` – Volcano plot comparing control and SARS-CoV-2 infected samples
* `GO_enrichment_24H_vs_72H.csv` – GO enrichment results for the 24H vs 72H comparison
* `GO_enrichment_Control_vs_Infected.csv` – GO enrichment results for control vs SARS-CoV-2 comparison
* GO enrichment bar plots visualizing significantly enriched biological processes

---

## Biological Interpretation

**Volcano Plot: 24H vs 72H**
The x-axis represents log2 fold change, and the y-axis represents −log10(p-value). Fewer significant DEGs were observed compared to the infection-based comparison, indicating more subtle transcriptional changes over time. Notably, upregulated genes at 72H included members of the **SNORD family** and other small nucleolar RNAs such as **SNORA5C** and **SCARNA8**. These findings suggest that temporal regulatory changes occur but are less pronounced than those induced by viral infection.

**Volcano Plot: Control vs SARS-CoV-2 Infected**
This comparison revealed a large number of significant DEGs, indicating widespread transcriptional disruption following infection. Upregulated genes included **SNORD116 family members** and **CENPF**, a gene involved in chromosome segregation during mitosis. **CXCL8**, a chemokine associated with immune response, was notably downregulated. These results are consistent with known viral infection mechanisms involving immune modulation and altered cell-cycle regulation.

**GO Enrichment Analysis**
The 24H vs 72H comparison yielded only one significantly enriched GO term (*ribonucleoprotein complex biogenesis*), consistent with limited transcriptional divergence over short time intervals. In contrast, the Control vs SARS-CoV-2 comparison revealed **165 enriched GO terms**, reflecting the strong biological impact of viral infection on host cellular processes.

Previous studies (Liu et al.) have reported upregulation of snoRNAs such as **SNORD116**, **SNORD66**, and **SNORD18** in chronic obstructive pulmonary disease (COPD), where they may contribute to inflammation and tissue remodeling. Similarly, this study identified significant upregulation of **SNORD116** in SARS-CoV-2–infected cells, suggesting a potential role in viral replication or host immune response. These findings indicate that dysregulated snoRNA expression may serve as **novel biomarkers or therapeutic targets** for COVID-19–associated respiratory diseases.

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
