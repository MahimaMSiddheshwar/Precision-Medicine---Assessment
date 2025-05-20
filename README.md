# Mahima-MS_PM_-Assignment-02


Scope
This assignment involves a comprehensive RNA-Seq data analysis pipeline to study gene expression in human respiratory cells infected with SARS-CoV-2. The workflow includes data preprocessing, quality control, alignment, quantification, differential expression analysis, and GO enrichment analysis to identify significant genes and biological pathways affected by the infection.

Assignment Description
This project uses publicly available RNA-Seq data to perform differential expression analysis and gene ontology (GO) enrichment analysis. The analysis compares control and infected samples across two time points (24H and 72H), highlighting changes in gene expression and potential biological processes impacted by SARS-CoV-2 infection.

Objective
The objective of this assignment is to process RNA-Seq data to identify differentially expressed genes (DEGs) and enriched biological pathways in SARS-CoV-2 infected vs. control samples at two different time points. This analysis can reveal insights into host response mechanisms and cellular pathways affected by viral infection.

Programming Language
Bash (for data processing, alignment, and quantification)
R (for differential expression analysis and GO enrichment)

Input Files
SRR Files: RNA-Seq data from SRA (e.g., SRR22269872, SRR22269873, etc.)
Reference Genome: Homo_sapiens_UCSC_hg38.tar.gz
Annotation File: genes.gtf from UCSC (hg38)
gene_counts File: DEG/GO Analysis

Required Files/Dependencies
FASTQ Files: Generated from SRR files using SRA Toolkit (fasterq-dump)
Reference Genome and GTF Files: Downloaded from UCSC


Required Libraries

Modules:
sra-toolkit: For downloading and converting SRA files to FASTQ.
fastqc: For quality control of raw and trimmed reads.
STAR: For aligning reads to the reference genome.
samtools: For indexing and filtering BAM files.
subread: For gene quantification with FeatureCounts.
multiqc: For aggregating quality control reports.

Python Packages: MultiQC (installed via pip)

R Packages:
edgeR, limma: For DEG analysis.
clusterProfiler, org.Hs.eg.db: For GO enrichment analysis.
EnhancedVolcano: For visualizing DEGs.


# Execution Step
This workflow for RNAseq pipeline involves main steps, including alignment and quantification of transcriptome profiles for control and SARS-CoV-2 infected respiratory cells. Here is the step-by-step work flow for the following:

## Data Acquisition and Quality Control
1. Downloading the required samples using SRA toolkit.
2. Converted SRA files to FASTQ format.
3. Performed the quality control using FastQC and MultiQC.

## Data Preprocessing
Trimmed the reads using Trim Galore with quality and length thresholds.
and then performed quality control on trimmed reads using FastQC.

## Reference Genome Preparation
1. Downloaded the human reference genome (UCSC hg38).
2. Created a STAR index for the reference genome.

## Read Alignment
1. Aligned the trimmed reads to the reference genome using STAR.
2. Sorted and indexed the resulting BAM files.

## Gene Expression Quantification
1. Used featureCounts to quantify gene expression levels.
2. Repeated the process with primary alignments only.


## Codes for the following steps are as follows:
### Note: Path to the File Directory is Hardcoded please update the File before running the script

## Change directory to slate
cd /N/slate/msiddhe

## Create a new directory for Assignment-2
mkdir Assignment-2
cd Assignment-2


## Load the SRA toolkit module to download the files
module load sra-toolkit

## Download the required samples
prefetch SRR22269883 SRR22269882 SRR22269881 SRR22269880 SRR22269879 SRR22269878 SRR22269877 SRR22269876 SRR22269875 SRR22269874 SRR22269873 SRR22269872

## Convert the downloaded SRA files to fastq format
fasterq-dump SRR22269883 SRR22269882 SRR22269881 SRR22269880 SRR22269879 SRR22269878 SRR22269877 SRR22269876 SRR22269875 SRR22269874 SRR22269873 SRR22269872

## Create a new directory for Raw Data
mkdir Raw_Data
mv /N/slate/msiddhe/Assignment-2/*.fastq /N/slate/msiddhe/Assignment-2/Raw_Data/
mv /N/slate/msiddhe/Assignment-2/* /N/slate/msiddhe/Assignment-2/Raw_Data/

## Load FastQC module
module load fastqc

## Create a directory for FastQC reports
mkdir fastqc_reports

## Run FastQC on all fastq files
fastqc /N/slate/msiddhe/Assignment-2/Raw_Data/*.fastq -o /N/slate/msiddhe/Assignment-2/fastqc_reports

## Navigate to your quartz to install multiqc
$ cd /geode2/home/u060/msiddhe/Quartz

## Make a directory called installations
$ mkdir installations
$ cd installations

## Install multiqc
$ module load python
$ pip install multiqc

## Add multiqc to your path 
$ pip show multiqc
$ nano ~/.bashrc
$ export PATH="/geode2/home/u060/msiddhe/Quartz/.local/bin:$PATH"
$ source ~/.bashrc

## Check if it is added to your path
$ multiqc --version
cd /N/slate/msiddhe/Assignment-2

## Run multiqc on raw data
$ multiqc /N/slate/msiddhe/Assignment-2/fastqc_reports
$ ls

## Trimming Step
$ conda install -c bioconda trim-galore

$ mkdir  -p /N/slate/msiddhe/Assignment-2/trimmed_reads
$ cd /N/slate/msiddhe/Assignment-2/Raw_Data
$ trim_galore --quality 20 --length 30 -o /N/slate/msiddhe/Assignment-2/trimmed_reads /N/slate/msiddhe/Assignment-2/Raw_Data/*.fastq

$ mkdir -p /N/slate/msiddhe/Assignment-2/trimmed_fastqc_reports
$ cd /N/slate/msiddhe/Assignment-2/trimmed_reads
$ fastqc /N/slate/msiddhe/Assignment-2/trimmed_reads/*.fq -o /N/slate/msiddhe/Assignment-2/trimmed_fastqc_reports


## Download referece genome and genes.gtf file
$ cd ..
$ wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
$ tar -xzf Homo_sapiens_UCSC_hg38.tar.gz


## Star Indexing Step
$ mkdir /N/slate/msiddhe/Assignment-2/star_index
$ nano star_index.slurm
$ sbatch star_index.slurm
$ squeue -u msiddhe

----------------#############-------------
#!/bin/bash

#SBATCH --mail-user=msiddhe@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=2
#SBATCH --partition=gpu
#SBATCH --mem=150gb
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=indexing
#SBATCH --output=indexing_output
#SBATCH --error=indexing_error
#SBATCH -A c01064

## Load the STAR module
module load star/2.7.11a

## Run STAR genomeGenerate command
STAR --runThreadN 11 \
     --runMode genomeGenerate \
     --genomeDir /N/slate/msiddhe/Assignment-2/star_index \
     --genomeFastaFiles /N/slate/msiddhe/Assignment-2/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --sjdbGTFfile /N/slate/msiddhe/Assignment-2/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf
----------------#############-------------


## READ Alignment
$ cd /N/slate/msiddhe/Assignment-2/aligned_reads
$ nano aligned_reads.slurm
$ sbatch aligned_reads.slurm

----------------#############-------------
#!/bin/bash

#SBATCH --mail-user=msiddhe@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=11
#SBATCH --partition=general
#SBATCH --mem=150gb
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=align_reads
#SBATCH --output=align_reads_output
#SBATCH --error=align_reads_error
#SBATCH -A c01064

## Load the STAR module
module load star/2.7.11a

## Define input and output directories
input_dir="/N/slate/msiddhe/Assignment-2/trimmed_reads"
output_dir="/N/slate/msiddhe/Assignment-2/aligned_reads"
genome_dir="/N/slate/msiddhe/Assignment-2/star_index"

## Make the output directory if it doesn't exist
mkdir -p "$output_dir"

## Loop over each FASTQ file in the input directory
for read_file in "$input_dir"/*.fq; do
    # Extract the sample name from the file name
    sample_name=$(basename "$read_file" | cut -d '_' -f 1)

    # Run STAR alignment for the current FASTQ file
    STAR --runThreadN 11 \
         --genomeDir "$genome_dir" \
         --readFilesIn "$read_file" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "$output_dir/aligned_reads_${sample_name}_"
done
----------------#############------------------------------------------------------------------------

## BAM Files Preparation
$ cd /N/slate/msiddhe/Assignment-2
$ conda activate multiqc_env
$ nano index_bam_files.sh
$ chmod +x index_bam_files.sh
$ ./index_bam_files.sh

----------------#############-------------
#!/bin/bash

# Load samtools module if necessary (depends on your system)
module load samtools

# Define the directory containing BAM files
bam_dir="/N/slate/msiddhe/Assignment-2/aligned_reads"

# Loop over each BAM file in the directory and index it
for bam_file in "$bam_dir"/*.bam; do
    echo "Indexing $bam_file..."
    samtools index "$bam_file"
done

echo "Indexing complete for all BAM files."
----------------#############---------------------------------------------------------------


## Sorting Primary Bam Files (To ensure only primary bam files are aligned in order to ensures that read counts are accurate, unique, and interpretable, avoiding potential biases and complications from multiple mappings)
$ nano move_and_flagstat.sh
$ chmod +x move_and_flagstat.sh
$ ./move_and_flagstat.sh

----------------#############-------------	
#!/bin/bash

# Create a directory to store the primary BAM files
mkdir -p /N/slate/msiddhe/Assignment-2/aligned_reads/primary_bam_files

# Loop through each BAM file and filter primary alignments
for bam_file in /N/slate/msiddhe/Assignment-2/aligned_reads/*.bam; do
    base_name=$(basename "$bam_file" .bam)
    output_bam="/N/slate/msiddhe/Assignment-2/aligned_reads/primary_bam_files/${base_name}_primary.bam"
    samtools view -F 256 -o "$output_bam" "$bam_file"
    samtools flagstat "$output_bam" > "/N/slate/msiddhe/Assignment-2/aligned_reads/primary_bam_files/${base_name}_flagstat.txt"
done

echo "Process complete for all Primary_BAM files."
----------------#############---------------------------------------------------------------------------------------------------------------


## Feature Count
$ nano run_featurecounts.sh
$ chmod +x run_featurecounts.sh
$ ./run_featurecounts.sh

----------------#############-------------
#!/bin/bash

## Load the Subread module (FeatureCounts is part of Subread)
module load subread

## Define the directory containing BAM files and the output directory
bam_dir="/N/slate/msiddhe/Assignment-2/aligned_reads"
output_dir="/N/slate/msiddhe/Assignment-2/featurecounts_output"
mkdir -p "$output_dir"

## Path to the GTF file
gtf_file="/N/slate/msiddhe/Assignment-2/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"

## Run FeatureCounts
featureCounts -T 8 -t exon -g gene_id -a "$gtf_file" -o "$output_dir/gene_counts.txt" "$bam_dir"/*.bam

echo "Gene expression quantification completed. Results are in $output_dir/gene_counts.txt"
----------------#############-----------------------------------------------------------------------------------------

## New Feature Counts with Primary Bam Files (To ensure only primary bam files are aligned in order to ensures that read counts are accurate, unique, and interpretable, avoiding potential biases and complications from multiple mappings)
$ mkdir -p /N/slate/msiddhe/Assignment-2/NewFeatureCounts
$ chmod +x run_featurecounts.sh
$ ./run_featurecounts.sh

----------------#############-------------
#!/bin/bash

## Load the Subread module (FeatureCounts is part of Subread)
module load subread

## Define the directory containing primary BAM files and the output directory
primary_bam_dir="/N/slate/msiddhe/Assignment-2/aligned_reads/primary_bam_files"
output_dir="/N/slate/msiddhe/Assignment-2/NewFeatureCounts"
mkdir -p "$output_dir"

## Path to the GTF file
gtf_file="/N/slate/msiddhe/Assignment-2/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"

## Run FeatureCounts
featureCounts -T 8 -t exon -g gene_id -a "$gtf_file" -o "$output_dir/gene_counts.txt" "$primary_bam_dir"/*.bam

echo "Gene expression quantification completed. Results are in $output_dir/gene_counts.txt"
----------------#############---------------------------------------------------------------------------------------------------------------------

## Transfer Gene Count file to local machine
$ scp msiddhe@quartz.uits.iu.edu:/N/slate/msiddhe/Assignment-2/featurecounts_output/gene_counts.txt / C:/Users/mmsid/Downloads/
$ scp msiddhe@quartz.uits.iu.edu:/N/slate/msiddhe/Assignment-2/featurecounts_output/gene_counts.txt.summary / C:/Users/mmsid/Downloads/
$ scp msiddhe@quartz.uits.iu.edu:/N/slate/msiddhe/Assignment-2/NewFeatureCounts/gene_counts.txt C:/Users/mmsid/Downloads/
$ scp msiddhe@quartz.uits.iu.edu:/N/slate/msiddhe/Assignment-2/NewFeatureCounts/gene_counts.txt.summary C:/Users/mmsid/Downloads/  


# DEG/GO Enrichment Analysis

### Note: Path to the File Directory is Hardcoded please update the File before running the script


## Install BiocManager if you haven't installed it already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Install edgeR and other packages required
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")



## Set the working directory
setwd("C:/Users/mmsid/OneDrive/Documents/1. Mahima Mahabaleshwar_IUPUI Documents/1. IUPUI Semester Folder/3rd Semester IUI Fall 2024/2. Precesion Medicine/Assignment/Assignment-2/New file 2")

## Load necessary libraries
library(edgeR)            # For differential expression analysis
library(limma)            # Linear models for RNA-Seq data
library(org.Hs.eg.db)     # Annotation package for GO term enrichment (human data)
library(clusterProfiler)  # For GO enrichment analysis
library(EnhancedVolcano)  # For Volcano plots


## Checking if file exists and load gene counts data
if (!file.exists("gene_counts.txt")) stop("File 'gene_counts.txt' not found.")
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1)


## Columns with SRR identifiers (count data) and cleaning the column
counts_filtered <- counts[, grepl("^X\\.N\\.slate\\.msiddhe\\.Assignment\\.2\\.aligned_reads\\.primary_bam_files\\.aligned_reads_SRR", colnames(counts))]
colnames(counts_filtered) <- sub("_Aligned\\.sortedByCoord\\.out_primary\\.bam$", "", colnames(counts_filtered))

## Define sample groups
group <- factor(c(rep("Control_24H", 3), rep("Control_72H", 3), rep("Infected_24H", 3), rep("Infected_72H", 3)))


## Create DGEList object and normalize for library sizes
if (length(group) != ncol(counts_filtered)) stop("Mismatch between group length and number of count columns.")
y <- DGEList(counts = counts_filtered, group = group)
y <- calcNormFactors(y)

## Design matrix and fit the model
design <- model.matrix(~0 + group)  # Avoid intercept term
colnames(design) <- levels(group)
fit <- glmQLFit(y, design)


## Differential Expression Analysis

## 1. 24H vs 72H comparison
contrast_24h_vs_72h <- makeContrasts((Control_72H + Infected_72H) / 2 - (Control_24H + Infected_24H) / 2, levels = design)
result_24h_vs_72h <- glmQLFTest(fit, contrast = contrast_24h_vs_72h)
deg_24h_vs_72h <- topTags(result_24h_vs_72h, n = Inf)$table

## 2. Control vs Infected comparison
contrast_control_vs_infected <- makeContrasts((Infected_24H + Infected_72H) / 2 - (Control_24H + Control_72H) / 2, levels = design)
result_control_vs_infected <- glmQLFTest(fit, contrast = contrast_control_vs_infected)
deg_control_vs_infected <- topTags(result_control_vs_infected, n = Inf)$table


## Generating the Volcano Plots

## 1. 24H vs 72H Volcano Plot
png("volcano_24h_vs_72h.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(deg_24h_vs_72h,
                lab = rownames(deg_24h_vs_72h),
                x = 'logFC',
                y = 'PValue',
                title = 'Volcano Plot: 24H v/s 72H',
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 3.0,
                pointSize = 1.5,
                col = c("grey", "lightgreen", "blue", "lightcoral"))
dev.off()


## 2. Control vs Infected Volcano Plot
png("volcano_control_vs_infected.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(deg_control_vs_infected,
                lab = rownames(deg_control_vs_infected),
                x = 'logFC',
                y = 'PValue',
                title = 'Volcano Plot: Mock Control v/s SARS-CoV-2',
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 3.0,
                pointSize = 1.5,
                col = c("grey", "lightgreen", "blue", "lightcoral"))
dev.off()



## GO Enrichment Analysis

## Filtered DEGs with adjusted p-value < 0.05 and |logFC| > 1
deg_24h_vs_72h_filtered <- deg_24h_vs_72h[deg_24h_vs_72h$PValue < 0.05 & abs(deg_24h_vs_72h$logFC) > 1, ]
deg_control_vs_infected_filtered <- deg_control_vs_infected[deg_control_vs_infected$PValue < 0.05 & abs(deg_control_vs_infected$logFC) > 1, ]

## Converting gene symbols to Entrez IDs
genes_24h_vs_72h <- bitr(rownames(deg_24h_vs_72h_filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes_control_vs_infected <- bitr(rownames(deg_control_vs_infected_filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


## Performing GO Enrichment for each comparison

## 1. 24H vs 72H GO Enrichment
go_enrich_24h_vs_72h <- enrichGO(
  gene = genes_24h_vs_72h$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

## Ploting GO Enrichment Barplot for 24H vs 72H
barplot(go_enrich_24h_vs_72h, showCategory = 10, title = "GO Enrichment: 24H v/s 72H", font.size = 10)


## 2. Control vs Infected GO Enrichment
go_enrich_control_vs_infected <- enrichGO(
  gene = genes_control_vs_infected$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

## Ploting GO Enrichment Barplot for Control vs Infected
barplot(go_enrich_control_vs_infected, showCategory = 10, title = "GO Enrichment: Control v/s SARS-CoV-2", font.size = 10)


## Top 10 enriched GO terms
head(go_enrich_24h_vs_72h, 10)
head(go_enrich_control_vs_infected, 10)

## Save GO enrichment results to CSV files
write.csv(as.data.frame(go_enrich_24h_vs_72h), "GO_enrichment_24H_vs_72H.csv", row.names = FALSE)
write.csv(as.data.frame(go_enrich_control_vs_infected), "GO_enrichment_Control_vs_Infected.csv", row.names = FALSE)


## Output Files:
volcano_24h_vs_72h.png: Volcano plot comparing gene expression between 24H and 72H.
volcano_control_vs_infected.png: Volcano plot comparing gene expression between control and SARS-CoV-2 infected cells.
GO_enrichment_24H_vs_72H.csv: GO enrichment results for 24H vs 72H comparison.
GO_enrichment_Control_vs_Infected.csv: GO enrichment results for control vs SARS-CoV-2 comparison.
GO_enrichment_24H_vs_72H.png: Barplot showing enriched GO terms for the 24H vs 72H comparison.
GO_enrichment_Control_vs_SARS-CoV-2.png: Barplot showing enriched GO terms for the control vs SARS-CoV-2 comparison.


## Interpretation for each output files:

Volcano Plot: 24H vs 72H
The x-axis represents log2 fold change, and the y-axis represents -log10P. Genes with significant changes between time points appear as red dots. There are fewer significant DEGs compared to the control vs infected comparison, indicating more subtle changes over time. Notable upregulated genes at 72H include members of the SNORD family and other small nucleolar RNAs (SNORA5C, SCARNA8).
In conclusion the volcano plot shows that there are fewer significant DEGs between time points compared to infection status. This suggests that while some regulatory changes occur over time, they may not be as pronounced as those caused by viral infection itself.

Volcano Plot: Control (mock) vs SARS-CoV-2 infected cells
The x-axis represents the log2 fold change (log2FC) in gene expression between conditions, while the y-axis represents the negative log10 of the p-value (-log10P). Genes farther to the right are upregulated in SARS-CoV-2 infected cells, while those farther to the left are downregulated. Genes with both a high fold change and low p-value appear as red dots (significant DEGs).
Some notable genes include: SNORD116 family members (upregulated).CENPF (upregulated), which is involved in chromosome segregation during mitosis. CXCL8 (downregulated), a chemokine involved in immune response. The plot shows a large number of significant DEGs, indicating that SARS-CoV-2 infection induces widespread changes in gene expression.
In conclusion this plot highlights key genes that are differentially expressed between control and SARS-CoV-2 infected cells. Many of these genes are involved in cell division and immune response, consistent with viral infection mechanisms.

The GO Enrichment: 24H vs 72H plot shows only one significant GO term ("ribonucleoprotein complex biogenesis"), which aligns with the expectation that miRNAs may not induce widespread changes over short time periods (like 24H vs 72H). This significantly explain why only one enriched term was observed.

In contrast, the GO Enrichment: Control vs SARS-CoV-2 plot shows more enriched GO terms (165 terms), which is expected because viral infections like SARS-CoV-2 tend to cause more pronounced gene expression changes, even at the level of miRNA regulation.


Liu et al., reported that snoRNAs such as SNORD116, SNORD66, and SNORD18 have been found to be upregulated in COPD (chronic obstructive pulmonary disease) patients, potentially contributing to chronic inflammation and tissue remodeling. Similarly from my findings in COVID-19, snoRNAs like SNORD116 are significantly upregulated, suggesting that they may be involved in viral replication or host immune responses that may infect the human respiratory cells severly. These findings indicate that dysregulated snoRNA expression could serve as novel biomarkers for lung diseases and offer promising therapeutic targets for managing conditions  related to COVID-19.


---------------------------------------------------------------################################################--------------------------------------------------------------------------



