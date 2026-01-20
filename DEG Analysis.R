
# Install BiocManager if you haven't installed it already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install edgeR and other packages required
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")



# Set the working directory
setwd("C:/Users/mmsid/OneDrive/Documents/1. Mahima Mahabaleshwar_IUPUI Documents/1. IUPUI Semester Folder/3rd Semester IUI Fall 2024/2. Precesion Medicine/Assignment/Assignment-2/New file 2")

# Load necessary libraries
library(edgeR)            # For differential expression analysis
library(limma)            # Linear models for RNA-Seq data
library(org.Hs.eg.db)     # Annotation package for GO term enrichment (human data)
library(clusterProfiler)  # For GO enrichment analysis
library(EnhancedVolcano)  # For Volcano plots


# Checking if file exists and load gene counts data
if (!file.exists("gene_counts.txt")) stop("File 'gene_counts.txt' not found.")
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1)


# Columns with SRR identifiers (count data) and cleaning the column
counts_filtered <- counts[, grepl("^X\\.N\\.slate\\.msiddhe\\.Assignment\\.2\\.aligned_reads\\.primary_bam_files\\.aligned_reads_SRR", colnames(counts))]
colnames(counts_filtered) <- sub("_Aligned\\.sortedByCoord\\.out_primary\\.bam$", "", colnames(counts_filtered))

# Define sample groups
group <- factor(c(rep("Control_24H", 3), rep("Control_72H", 3), rep("Infected_24H", 3), rep("Infected_72H", 3)))


# Create DGEList object and normalize for library sizes
if (length(group) != ncol(counts_filtered)) stop("Mismatch between group length and number of count columns.")
y <- DGEList(counts = counts_filtered, group = group)
y <- calcNormFactors(y)

# Design matrix and fit the model
design <- model.matrix(~0 + group)  # Avoid intercept term
colnames(design) <- levels(group)
fit <- glmQLFit(y, design)


# Differential Expression Analysis

# 1. 24H vs 72H comparison
contrast_24h_vs_72h <- makeContrasts((Control_72H + Infected_72H) / 2 - (Control_24H + Infected_24H) / 2, levels = design)
result_24h_vs_72h <- glmQLFTest(fit, contrast = contrast_24h_vs_72h)
deg_24h_vs_72h <- topTags(result_24h_vs_72h, n = Inf)$table

# 2. Control vs Infected comparison
contrast_control_vs_infected <- makeContrasts((Infected_24H + Infected_72H) / 2 - (Control_24H + Control_72H) / 2, levels = design)
result_control_vs_infected <- glmQLFTest(fit, contrast = contrast_control_vs_infected)
deg_control_vs_infected <- topTags(result_control_vs_infected, n = Inf)$table


# Generating the Volcano Plots

# 1. 24H vs 72H Volcano Plot
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


# 2. Control vs Infected Volcano Plot
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



# GO Enrichment Analysis

# Filtered DEGs with adjusted p-value < 0.05 and |logFC| > 1
deg_24h_vs_72h_filtered <- deg_24h_vs_72h[deg_24h_vs_72h$PValue < 0.05 & abs(deg_24h_vs_72h$logFC) > 1, ]
deg_control_vs_infected_filtered <- deg_control_vs_infected[deg_control_vs_infected$PValue < 0.05 & abs(deg_control_vs_infected$logFC) > 1, ]

# Converting gene symbols to Entrez IDs
genes_24h_vs_72h <- bitr(rownames(deg_24h_vs_72h_filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes_control_vs_infected <- bitr(rownames(deg_control_vs_infected_filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Performing GO Enrichment for each comparison

# 1. 24H vs 72H GO Enrichment
go_enrich_24h_vs_72h <- enrichGO(
  gene = genes_24h_vs_72h$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Ploting GO Enrichment Barplot for 24H vs 72H
barplot(go_enrich_24h_vs_72h, showCategory = 10, title = "GO Enrichment: 24H v/s 72H", font.size = 10)


# 2. Control vs Infected GO Enrichment
go_enrich_control_vs_infected <- enrichGO(
  gene = genes_control_vs_infected$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Ploting GO Enrichment Barplot for Control vs Infected
barplot(go_enrich_control_vs_infected, showCategory = 10, title = "GO Enrichment: Control v/s SARS-CoV-2", font.size = 10)


# Top 10 enriched GO terms
head(go_enrich_24h_vs_72h, 10)
head(go_enrich_control_vs_infected, 10)

# Save GO enrichment results to CSV files
write.csv(as.data.frame(go_enrich_24h_vs_72h), "GO_enrichment_24H_vs_72H.csv", row.names = FALSE)
write.csv(as.data.frame(go_enrich_control_vs_infected), "GO_enrichment_Control_vs_Infected.csv", row.names = FALSE)
