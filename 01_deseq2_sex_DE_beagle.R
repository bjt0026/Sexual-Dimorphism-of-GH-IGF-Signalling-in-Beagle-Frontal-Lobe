# DESeq2 pipeline for sex differential expression in beagle frontal cortex
# Course: Functional Genomics (BIOL 6850)
# Author: Beverly Thomas
# Input: gene_count_matrix.csv (gene-level output)
# Output: DESeq2 results tables + GH/IGF figures (600 DPI PNG)
# This script:
#   1) runs DESeq2 (male vs female),
#   2) generates QC plots (sample distance heatmap, MA plot, PCA plot),
#   3) extracts GH/IGF pathway genes,
#   4) produces hypothesis-driven heatmaps and boxplots,
#   5) writes a ranked .rnk file and normalized expression table
#      formatted for downstream tools (e.g., GSEA, Cytoscape).

library(DESeq2)
library(matrixStats)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)

### 1. Set working directory first (in RStudio: Session -> Set Working Directory)
# setwd("~/scratch")

##########  Input data ##############

### Input the count data
# Expectation: gene_count_matrix.csv with:
#   - rows = genes,
#   - columns = samples,
#   - first column = gene IDs,
#   - raw integer counts
countdata <- as.matrix(read.csv("gene_count_matrix.csv",
                                row.names = 1,
                                check.names = FALSE))
dim(countdata)
head(countdata)

### Build meta data (phenotype data) inline: Male vs Female
# Define sample-level information. Here, we only model sex.
sample <- colnames(countdata)

sex <- c(
  "female","female","female", # SAMN31934616, SAMN31934617, SAMN31934618
  "male","male","male"        # SAMN31934759, SAMN31934760, SAMN31934761
)

coldata <- data.frame(
  sample = sample,
  sex = factor(sex, levels = c("female","male")),  # set reference level order
  row.names = sample
)

dim(coldata)
head(coldata)

# Sanity checks: column order in count matrix must match row order in coldata.
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

## Create DESeq dataset and define model
# Design: counts ~ sex (female vs male).
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ sex
)

##### Prefiltering
# Remove very low-count genes to speed up DESeq2 and reduce multiple testing.
dds <- dds[ rowSums(counts(dds)) > 20, ]
dds

##### Differential expression analysis: male vs female
# Relevel so that female is baseline and log2FC > 0 = higher in males.
dds$sex <- relevel(dds$sex, ref = "female") # female = reference
dds <- DESeq(dds)

# Wald test results for male vs female.
res <- results(dds, contrast = c("sex","male","female"))
res

### Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]
resOrdered

### DEG tallies at different stringency levels
deg_05 <- res[ which(res$padj < 0.05), ]
deg_10 <- res[ which(res$padj < 0.10), ]
deg_20 <- res[ which(res$padj < 0.20), ]

cat("DEGs at FDR 0.05 (padj < 0.05):", nrow(deg_05), "\n")
cat("DEGs at FDR 0.10 (padj < 0.10):", nrow(deg_10), "\n")
cat("DEGs at FDR 0.20 (padj < 0.20):", nrow(deg_20), "\n")

# Exploratory subset: nominal p < 0.05 AND |log2FC| >= 1.
cand_p05_lfc1 <- res[ which(res$pvalue < 0.05 &
                              !is.na(res$pvalue) &
                              abs(res$log2FoldChange) >= 1), ]
cat("Genes with p < 0.05 and |log2FC| >= 1 (exploratory):",
    nrow(cand_p05_lfc1), "\n")

### MA-plot (on screen only)
# Quick visual QC: log2FC vs mean expression.
plotMA(res, main = "DESeq2: male vs female", ylim = c(-8,8))

## Plot counts for top DEG
# Counts per sample for lowest FDR gene.
plotCounts(dds, gene = which.min(res$padj), intgroup = "sex")

## Write ordered results
# Full DESeq2 result table, ordered by FDR.
write.csv(as.data.frame(resOrdered),
          file = "DGESeq_results_male_vs_female.csv")

## Transformations
# rlog / VST for downstream visualization (heatmaps, PCA, distances).
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

### Optional overview heatmaps (top DE, top variable) ----
# Top 20 DE genes by FDR.
top20_ids <- rownames(resOrdered)[which(!is.na(resOrdered$padj))[1:20]]
mat_top20 <- assay(vsd)[top20_ids, ]
# Center each gene across samples.
mat_top20 <- mat_top20 - rowMeans(mat_top20)
anno_col <- as.data.frame(colData(dds)[, "sex", drop = FALSE])

pheatmap(
  mat_top20,
  annotation_col = anno_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = " "
)

# Top 20 most variable genes across all samples.
topVarGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 20)
mat_var <- assay(vsd)[ topVarGenes, ]
mat_var <- mat_var - rowMeans(mat_var)
anno_var <- as.data.frame(colData(vsd)[, "sex", drop = FALSE])

pheatmap(
  mat_var,
  annotation_col = anno_var,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = " "
)

# GH–IGF / insulin / downstream signalling gene set
# Based on canonical GH–IGF and JAK–STAT / PI3K–AKT–mTOR / MAPK pathway
# schematics and reviews, for example:
# Le, T. K. C., Dao, X. D., Nguyen, D. V., Luu, D. H., Bui, T. M. H., Le, T. H., Nguyen, H. T., Le, T. N., Hosaka, T., & Nguyen, T. T. T. (2023). Insulin signaling and its application. Frontiers in endocrinology, 14, 1226655. https://doi.org/10.3389/fendo.2023.1226655

gh_symbols <- c(
  "GH1","GH2","GHR","GHRHR","GHRH",
  "IGF1","IGF2","IGF1R","IGF2R",
  "IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IGFALS",
  "JAK2","STAT5A","STAT5B","STAT3",
  "IRS1","IRS2","SHC1","GRB2","SOS1",
  "PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2",
  "AKT1","AKT2","AKT3",
  "TSC1","TSC2","MTOR","RPTOR","RICTOR",
  "FOXO1","FOXO3","FOXO4","GSK3B","BAD",
  "HRAS","KRAS","NRAS","RAF1","BRAF",
  "MAP2K1","MAP2K2","MAPK1","MAPK3","ELK1","JUN","FOS",
  "SOCS1","SOCS2","SOCS3","CISH","PTPN1","PTPN2",
  "MAS1","SLC22A1","SLC22A2","SLC22A3","ACAT2","RPS6KA2"
)

### Map DESeq2 results to GH/IGF pathway genes -------------
# Add gene_id and symbol columns to DESeq2 table.
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$symbol <- sub(".*\\|", "", res_df$gene_id)

# Subset to GH/IGF pathway list.
gh_res <- res_df[ res_df$symbol %in% gh_symbols, ]
gh_res <- gh_res[order(gh_res$padj), ]

# Save full hypothesis table
write.csv(gh_res,
          file = "GH_IGF_hypothesis_genes_results.csv",
          row.names = FALSE)

### Exploratory GH/IGF subsets ---------------
# FDR cutoffs for GH/IGF hypothesis set.
g050 <- subset(gh_res, !is.na(padj) & padj < 0.05)
gh10 <- subset(gh_res, !is.na(padj) & padj < 0.10)
gh20 <- subset(gh_res, !is.na(padj) & padj < 0.20)

write.csv(gh10, file = "GH_IGF_genes_FDR0.05.csv", row.names = FALSE)
write.csv(gh10, file = "GH_IGF_genes_FDR0.10.csv", row.names = FALSE)
write.csv(gh20, file = "GH_IGF_genes_FDR0.20.csv", row.names = FALSE)

### Sample-to-sample distances + MA plot + PCA -----------
# QC Panel:
#   Fig_QC_SampleDistances: Euclidean distances between samples (rlog scale).
#   Fig_QC_MAplot: MA plot 
#   Fig_QC_PCA_sex: PCA of rlog-transformed counts colored by sex.

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sex)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Distance heatmap PNG (QC 1)
png("Fig_QC_SampleDistances.png", width = 7, height = 6, units = "in", res = 600)
sampleDistPlot <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors,
                           main = "Sample-to-sample distances")
print(sampleDistPlot)
dev.off()

# MA plot PNG (QC 2)
png("Fig_QC_MAplot.png", width = 7, height = 6, units = "in", res = 600)
plotMA(res, main = "DESeq2: male vs female", ylim = c(-8,8))
dev.off()

# PCA PNG (QC 3)
png("Fig_QC_PCA_sex.png", width = 7, height = 6, units = "in", res = 600)
print(plotPCA(rld, intgroup = "sex"))
dev.off()

### =========================================================
### Hypothesis-targeted FIGURES
###Fig 1: MA, PCA, sample-to-sample distances
### Fig2: Key GH/IGF genes (boxplots)
### Fig3A: GH/IGF FDR < 0.05
### Fig3B: GH/IGF FDR < 0.20
### =========================================================

## Fig3 ---------------------
# Heatmap of GH/IGF genes with FDR < 0.05 (most stringent set).

gh05 <- subset(gh_res, !is.na(padj) & padj < 0.05)
gh05_ids <- gh05$gene_id[gh05$gene_id %in% rownames(vsd)]

if (length(gh05_ids) > 1) {
  mat_gh05 <- assay(vsd)[gh05_ids, , drop = FALSE]
  mat_gh05 <- mat_gh05 - rowMeans(mat_gh05)
  row_syms05 <- gh05$symbol[match(gh05_ids, gh05$gene_id)]
  rownames(mat_gh05) <- make.unique(row_syms05)
  
  anno_col4 <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
  sex_cols <- list(sex = c("female" = "salmon", "male" = "steelblue"))
  
  png("Fig1_GH_IGF_FDR0.05_heatmap.png",
      width = 7, height = 4.5, units = "in", res = 600)
  pheatmap(
    mat_gh05,
    annotation_col    = anno_col4,
    annotation_colors = sex_cols,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 8,
    main              = " "
  )
  dev.off()
} else {
  message("No GH/IGF genes with padj < 0.05.")
}

## Fig3B ---------------------
# Heatmap of GH/IGF genes with FDR < 0.20 (exploratory set).

gh20_ids <- gh20$gene_id[gh20$gene_id %in% rownames(vsd)]

if (length(gh20_ids) > 1) {
  mat_gh20 <- assay(vsd)[gh20_ids, , drop = FALSE]
  mat_gh20 <- mat_gh20 - rowMeans(mat_gh20)
  row_syms20 <- gh20$symbol[match(gh20_ids, gh20$gene_id)]
  rownames(mat_gh20) <- make.unique(row_syms20)
  
  anno_col3 <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
  # Consistent sex colors: female = salmon, male = steelblue.
  sex_cols <- list(sex = c("female" = "salmon", "male" = "steelblue"))
  
  png("Fig2B_GH_IGF_FDR0.20_exploratory_heatmap.png",
      width = 7, height = 6, units = "in", res = 600)
  pheatmap(
    mat_gh20,
    annotation_col = anno_col3,
    annotation_colors = sex_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = " "
  )
  dev.off()
}

## ---------------------
# Heatmap of GH/IGF genes with FDR < 0.10 (middle-set).
####JUST FOR FUN!####

gh10_ids <- gh10$gene_id[gh10$gene_id %in% rownames(vsd)]

if (length(gh10_ids) > 1) {
  mat_gh10 <- assay(vsd)[gh10_ids, , drop = FALSE]
  mat_gh10 <- mat_gh10 - rowMeans(mat_gh10)
  row_syms10 <- gh10$symbol[match(gh10_ids, gh10$gene_id)]
  rownames(mat_gh10) <- make.unique(row_syms10)
  
  # annotation and colors for sex
  anno_col2 <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
  sex_cols <- list(sex = c("female" = "salmon", "male" = "steelblue"))
  
  png("Fig2A_GH_IGF_FDR0.10_heatmap.png",
      width = 7, height = 4.5, units = "in", res = 600)
  pheatmap(
    mat_gh10,
    annotation_col = anno_col2,
    annotation_colors = sex_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = " "
  )
  dev.off()
}


## Fig2: Box/strip plots for key GH/IGF genes ----
# Normalized-count boxplots for selected GH/IGF genes (per sex).

key_genes2 <- c("FOS", "IGF1", "IGF2", "GHR", "IGF1R", "STAT5B", "SOCS2")

norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$symbol <- sub(".*\\|", "", norm_counts_df$gene_id)

plot_genes_df2 <- norm_counts_df[norm_counts_df$symbol %in% key_genes2, ]

if (nrow(plot_genes_df2) > 0) {
  long_df2 <- melt(
    plot_genes_df2,
    id.vars = c("gene_id", "symbol"),
    variable.name = "sample",
    value.name = "normalized_count"
  )
  
  long_df2$sex <- coldata[match(long_df2$sample, rownames(coldata)), "sex"]
  
  p2 <- ggplot(long_df2, aes(x = sex, y = normalized_count, color = sex)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.08, size = 2.5) +
    facet_wrap(~symbol, scales = "free_y") +
    theme_bw(base_size = 12) +
    labs(
      title = "Key GH/IGF genes including FOS",
      x = "Sex",
      y = "Normalized counts"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  png("Fig3_Key_GH_IGF_genes.png", width = 10, height = 6, units = "in", res = 600)
  print(p2)
  dev.off()
}

### ---- Save key DESeq2 outputs to working directory -------
# Export key DESeq2 tables for downstream analyses and sharing.

write.csv(as.data.frame(resOrdered),
          file = "DESeq2_results_male_vs_female_ordered.csv")

write.csv(as.data.frame(deg_05),
          file = "DEGs_FDR_0.05.csv")
write.csv(as.data.frame(deg_10),
          file = "DEGs_FDR_0.10.csv")
write.csv(as.data.frame(deg_20),
          file = "DEGs_FDR_0.20.csv")

write.csv(as.data.frame(cand_p05_lfc1),
          file = "DEGs_p0.05_log2FC1_exploratory.csv")

top20_table <- resOrdered[top20_ids, ]
write.csv(as.data.frame(top20_table),
          file = "Top20_DEGs_for_heatmap.csv")

# already wrote gh_res, gh10, gh20 above

### 1) Ranked list with symbols only  ----
# Build a rank metric for GSEA-style analyses:
#   rank = sign(log2FC) * -log10(pvalue)
# and output as a two-column .rnk file with gene symbols only.

res_df <- as.data.frame(resOrdered)
res_df$gene_id <- rownames(res_df)

# SYMBOL ONLY from rowname 
res_df$symbol <- sub(".*\\|", "", res_df$gene_id)

res_rank <- within(res_df,
                   rank <- sign(log2FoldChange) * -log10(pvalue))

DGErank_withName <- subset(res_rank,
                           !is.na(symbol) & symbol != "" & !is.na(rank),
                           select = c(symbol, rank))
# Assumes row names are formatted like geneID|SYMBOL; if row names are already symbols,
# this still returns the full row name.
## Use symbol as Name, matching OG column name
colnames(DGErank_withName)[1] <- "Name"

write.table(as.data.frame(DGErank_withName),
            file = "DGErankName.rnk",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

### 2) Normalized expression table, symbols only, OG filename ----
# Create log2-normalized expression matrix (log2(counts+1)) with:
#   - first column = Name (symbol),
#   - second column = gene_id (symbol),
#   - remaining columns = samples.
# This matches expected input format for many network / pathway tools.

nt <- normTransform(dds)  # log2(x + 1)
NormTransExp <- as.data.frame(assay(nt))

# SYMBOL ONLY from rowname
symbol <- sub(".*\\|", "", rownames(NormTransExp))

# Both Name and gene_id = symbol (no prefixes)
NormTransExp_Anno <- cbind(
  Name    = symbol,
  gene_id = symbol,
  NormTransExp
)

# Drop rows without symbol
NormTransExp_Anno_withName <- subset(NormTransExp_Anno,
                                     !is.na(Name) & Name != "")

write.table(as.data.frame(NormTransExp_Anno_withName),
            file = "NormTransExp_Anno_Names.txt",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")
