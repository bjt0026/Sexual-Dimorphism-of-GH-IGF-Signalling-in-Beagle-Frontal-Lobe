# DESeq2 pipeline for sex differential expression in beagle frontal cortex
# Course: Functional Genomics (BIOL 6850)
# Author: Beverly Thomas
# Input: gene_count_matrix.csv (prepDE.py gene-level output)
# Output: DESeq2 results tables + GH/IGF figures (PDF)

library(DESeq2)
library(matrixStats)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)

### 1. Set working directory first (in RStudio: Session -> Set Working Directory)
# setwd("~/scratch")

########## 1.3 Input data ##############

### Input the count data (from PrepDE.py output)
countdata <- as.matrix(read.csv("gene_count_matrix.csv",
                                row.names = 1,
                                check.names = FALSE))
dim(countdata)
head(countdata)

### Build meta data (phenotype data) inline: Male vs Female
sample <- colnames(countdata)

sex <- c(
  "female","female","female", # SAMN31934616, SAMN31934617, SAMN31934618
  "male","male","male" # SAMN31934759, SAMN31934760, SAMN31934761
)

coldata <- data.frame(
  sample = sample,
  sex = factor(sex, levels = c("female","male")),
  row.names = sample
)

dim(coldata)
head(coldata)

# Check that sample IDs match and are in the same order
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

## Create DESeq dataset and define model
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ sex
)

##### Prefiltering
dds <- dds[ rowSums(counts(dds)) > 20, ]
dds

##### Differential expression analysis: male vs female
dds$sex <- relevel(dds$sex, ref = "female") # female = reference
dds <- DESeq(dds)

res <- results(dds, contrast = c("sex","male","female"))
res

### Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]
resOrdered

### DEG tallies at different stringency levels
deg_05 <- res[ which(res$padj < 0.05), ]
deg_10 <- res[ which(res$padj < 0.10), ]
deg_20 <- res[ which(res$padj < 0.20), ]

cat("DEGs at FDR 0.05 (padj < 0.05):", nrow(deg_05), "\\n")
cat("DEGs at FDR 0.10 (padj < 0.10):", nrow(deg_10), "\\n")
cat("DEGs at FDR 0.20 (padj < 0.20):", nrow(deg_20), "\\n")

cand_p05_lfc1 <- res[ which(res$pvalue < 0.05 &
                              !is.na(res$pvalue) &
                              abs(res$log2FoldChange) >= 1), ]
cat("Genes with p < 0.05 and |log2FC| >= 1 (exploratory):",
    nrow(cand_p05_lfc1), "\\n")

### MA-plot (on screen only)
plotMA(res, main = "DESeq2: male vs female", ylim = c(-8,8))

## Plot counts for top DEG
plotCounts(dds, gene = which.min(res$padj), intgroup = "sex")

## Write ordered results
write.csv(as.data.frame(resOrdered),
          file = "DGESeq_results_male_vs_female.csv")

## Transformations
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

### Optional overview heatmaps (top DE, top variable) ----
top20_ids <- rownames(resOrdered)[which(!is.na(resOrdered$padj))[1:20]]
mat_top20 <- assay(vsd)[top20_ids, ]
mat_top20 <- mat_top20 - rowMeans(mat_top20)
anno_col <- as.data.frame(colData(dds)[, "sex", drop = FALSE])

pheatmap(
  mat_top20,
  annotation_col = anno_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 20 DE genes (male vs female)"
)

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
  main = "Top 20 most variable genes"
)

# GH–IGF / insulin / downstream signalling gene set
# Based on canonical GH–IGF and JAK–STAT / PI3K–AKT–mTOR / MAPK pathway
# schematics and reviews, for example:
# - GH–IGF axis and receptor signalling:
# PMID: 23443822 (Le Roith & Yakar, Endocrinology, IGF system overview)[web:309]
# PMID: 25757452 (Waters & Brooks, Growth hormone receptor signalling)[web:309]
# - JAK–STAT and SOCS negative feedback:
# PMID: 19112499 (Starr & Hilton, JAK–STAT signalling and SOCS proteins)[web:323]
# - PI3K–AKT–mTOR and FOXO:
# PMID: 17011494 (Manning & Cantley, AKT/PKB signalling)[web:314]
# PMID: 12150915 (Brunet et al., FOXO transcription factors)[web:314]
# - RAS–RAF–MEK–ERK MAPK cascade:
# PMID: 17496910 (Wellbrock et al., RAF proteins and ERK signalling)[web:323]

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

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$symbol <- sub(".*\\\\|", "", res_df$gene_id)

gh_res <- res_df[ res_df$symbol %in% gh_symbols, ]
gh_res <- gh_res[order(gh_res$padj), ]

# Save full hypothesis table
write.csv(gh_res,
          file = "GH_IGF_hypothesis_genes_results.csv",
          row.names = FALSE)

### Stringent and exploratory GH/IGF subsets ---------------

gh10 <- subset(gh_res, !is.na(padj) & padj < 0.10)
gh20 <- subset(gh_res, !is.na(padj) & padj < 0.20)

write.csv(gh10, file = "GH_IGF_genes_FDR0.10.csv", row.names = FALSE)
write.csv(gh20, file = "GH_IGF_genes_FDR0.20.csv", row.names = FALSE)

### Sample-to-sample distances + PCA (on screen) -----------

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sex)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-sample distances")

plotPCA(rld, intgroup = "sex")

### =========================================================
### Hypothesis-targeted FIGURES (saved to disk)
### Fig1: PCA by sex
### Fig2A: GH/IGF FDR < 0.10
### Fig2B: GH/IGF FDR < 0.20 (exploratory)
### Fig3: Key GH/IGF genes including FOS (boxplots)
### =========================================================

## Fig1: PCA by sex ----------------------------------------

pdf("Fig1_PCA_sex.pdf", width = 7, height = 6)
print(plotPCA(rld, intgroup = "sex"))
dev.off()

## Fig2A: GH/IGF genes with FDR < 0.10 ---------------------

gh10_ids <- gh10$gene_id[gh10$gene_id %in% rownames(vsd)]

if (length(gh10_ids) > 1) {
  mat_gh10 <- assay(vsd)[gh10_ids, , drop = FALSE]
  mat_gh10 <- mat_gh10 - rowMeans(mat_gh10)
  row_syms10 <- gh10$symbol[match(gh10_ids, gh10$gene_id)]
  rownames(mat_gh10) <- make.unique(row_syms10)
  
  anno_col2 <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
  
  pdf("Fig2A_GH_IGF_FDR0.10_heatmap.pdf", width = 7, height = 4.5)
  pheatmap(
    mat_gh10,
    annotation_col = anno_col2,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = "GH/IGF genes with FDR < 0.10"
  )
  dev.off()
}

## Fig2B: Exploratory GH/IGF genes with FDR < 0.20 ---------

gh20_ids <- gh20$gene_id[gh20$gene_id %in% rownames(vsd)]

if (length(gh20_ids) > 1) {
  mat_gh20 <- assay(vsd)[gh20_ids, , drop = FALSE]
  mat_gh20 <- mat_gh20 - rowMeans(mat_gh20)
  row_syms20 <- gh20$symbol[match(gh20_ids, gh20$gene_id)]
  rownames(mat_gh20) <- make.unique(row_syms20)
  
  anno_col3 <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
  
  pdf("Fig2B_GH_IGF_FDR0.20_exploratory_heatmap.pdf",
      width = 7, height = 6)
  pheatmap(
    mat_gh20,
    annotation_col = anno_col3,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = "Exploratory GH/IGF genes (FDR < 0.20)"
  )
  dev.off()
}

## Fig3: Box/strip plots for key GH/IGF genes incl. FOS ----

key_genes2 <- c("FOS", "IGF1", "IGF2", "GHR", "IGF1R", "STAT5B", "SOCS2")

norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$symbol <- sub(".*\\\\|", "", norm_counts_df$gene_id)

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
  
  ggsave("Fig3_Key_GH_IGF_genes.pdf",
         plot = p2, width = 10, height = 6)
}

### ---- Save key DESeq2 outputs to working directory -------

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

top20_table <- resOrdered[ top20_ids, ]
write.csv(as.data.frame(top20_table),
          file = "Top20_DEGs_for_heatmap.csv")

# already wrote gh_res, gh10, gh20 above

norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          file = "Normalized_counts_all_genes.csv")
