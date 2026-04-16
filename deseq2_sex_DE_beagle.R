# DESeq2 pipeline for sex differential expression in beagle frontal cortex
# Course: Functional Genomics (BIOL 6850)
# Author: Beverly Thomas
# Input: gene_count_matrix.csv (gene-level output)
# Output: DESeq2 results tables + QC Figures + GH/IGF figures (600 DPI PNG)
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

# Drop mitochondrial rna-NC_* features from the ordered results
is_mito_rna <- grepl("^rna-NC_", rownames(resOrdered))
cat("Number of rna-NC_ features removed:", sum(is_mito_rna), "\n")

resOrdered <- resOrdered[!is_mito_rna, ]

### DEG tallies at different stringency levels
deg_05 <- res[ which(res$padj < 0.05), ]
deg_10 <- res[ which(res$padj < 0.10), ]
deg_20 <- res[ which(res$padj < 0.20), ]

cat("DEGs at FDR 0.05 (padj < 0.05):", nrow(deg_05), "\n")
cat("DEGs at FDR 0.10 (padj < 0.10):", nrow(deg_10), "\n")
cat("DEGs at FDR 0.20 (padj < 0.20):", nrow(deg_20), "\n")

# Nominal p < 0.05 AND |log2FC| >= 1.
cand_p05_lfc1 <- res[ which(res$pvalue < 0.05 &
                              !is.na(res$pvalue) &
                              abs(res$log2FoldChange) >= 1), ]
cat("Genes with p < 0.05 and |log2FC| >= 1:",
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

### Optional GENERAL DGE overview heatmaps (top DE, top variable) ----


# Top 20 DE genes by FDR
top20_ids <- rownames(resOrdered)[which(!is.na(resOrdered$padj))[1:20]]
mat_top20 <- assay(vsd)[top20_ids, ]
mat_top20 <- mat_top20 - rowMeans(mat_top20)

row_syms_top20 <- sub(".*\\|", "", top20_ids)
rownames(mat_top20) <- make.unique(row_syms_top20)

anno_sex <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
sex_cols <- list(sex = c("female" = "salmon", "male" = "steelblue"))

png("Fig_QC_Top20_DE_heatmap.png", width = 7, height = 4.5, units = "in", res = 600)
pheatmap(
  mat_top20,                  # <- FIXED: use mat_top20 here
  annotation_col    = anno_sex,
  annotation_colors = sex_cols,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 6,
  main              = " "
)
dev.off()

# Top 20 most variable genes across all samples
topVarGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 20)
mat_var <- assay(vsd)[topVarGenes, ]
mat_var <- mat_var - rowMeans(mat_var)

row_syms_var <- sub(".*\\|", "", rownames(mat_var))
rownames(mat_var) <- make.unique(row_syms_var)

png("Fig_QC_Top20_variable_genes_heatmap.png", width = 7, height = 4.5, units = "in", res = 600)
pheatmap(
  mat_var,
  annotation_col    = anno_sex,
  annotation_colors = sex_cols,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 6,
  main              = " "
)
dev.off()

#A few notes:
##### Hypothesis-centered Figures
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
# The rest of the references are in the other script titled "IGF–insulin_pathway_fgsea_analysis"

# GH–IGF / insulin / downstream signalling gene set
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
  "MAS1","SLC22A1","SLC22A2","SLC22A3","ACAT2","RPS6KA2", "JAK1","JAK2","JAK3","MAPK1","MAPK3","MTOR","STAT3","TYK2",
  "TYK2",        
  "STAT1","STAT2",
  "STAT4","STAT5A","STAT5B","STAT6",   # STAT transcription factors
  "IL6","IL6R","IL10RA","IL10RB",      # key IL cytokines/receptors
  "IL2","IL2RA","IL2RB","IL2RG",       # IL‑2 receptor complex
  "OSMR","LIFR","IL7R",                # gp130/IL‑6 family receptors
  "SOCS1","SOCS2","SOCS3", "IGF1","IGF1R",
  "IRS1","IRS2",
  # Ras–MAPK arm
  "SHC1","GRB2","SOS1",
  "HRAS","RAF1","MAP2K1","MAPK3","MAPK8",
  "ELK1","FOS","JUN","SRF",
  # PI3K–AKT survival arm
  "PIK3CA","PIK3R1","PDPK1","AKT1","BAD",
  # modulation / GAP/PTP
  "PTPN11","RASA1", "IGF1","IGF2","IGF1R","IGF2R","INSR",
  "BAD","GRB2","HRAS","IGF1","IGF1R","IRS1","MAPK3",
  "PIK3CA","PIK3R1","RAF1","SHC1","SOS1","YWHAH",
  "GRB2","HRAS","INS","INSR","IRS1","PIK3CA","PIK3R1",
  "PTPN11","SHC1","SLC2A14","SOS1",
  "CISH","GH1","GHR","GRB2","HNF1A","INS","INSR","IRS1","JAK2",
  "PIK3CA","PIK3R1","PLCG1","PRKCA","PRKCB","PTPN6","RPS6KA1",
  "SHC1","SLC2A14","SOS1","SRF", "CAV1","CISH","DAB1","HGS","LEPROT","NEUROD1","NF2",
  "PPP2CA","PPP2R1A","PTPN2",
  "SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS6","SOCS7",
  "CALM1","CAMK2A","CCL5","CENPJ","CRLF3",
  "CSH1","CSH2","CSHL1",
  "CYP1B1","EP300","ERBB4",
  "F2","F2R",
  "GH1","GH2","GHR",
  "HES1","HES5",
  "IFNL1","IL10",
  # PI3K catalytic subunits
  "PIK3CA","PIK3CB","PIK3CD","PIK3CG",
  # PI3K regulatory subunits
  "PIK3R1","PIK3R2","PIK3R3","PIK3R5",
  # AKT isoforms
  "AKT1","AKT2","AKT3",
  # PDK1 / mTOR nodes
  "PDPK1","MTOR","RPS6KB1","RPS6KB2","EIF4EBP1",
  # canonical downstream substrates / integrators
  "GSK3B","FOXO1","FOXO3","TSC1","TSC2","RHEB",
  # RTK adaptor linking to PI3K
  "GAB1","FRS2",
  "AKT1","AKT2","AKT3",
  "PHLPP1","PHLPP2","PTEN",
  "THEM4","TRIB3", "INPP4B","TSC1","TSC2",
  "AGO1","AGO2","AGO3","AGO4",
  "AKT1","AKT1S1","AKT2","AKT3",
  "BAD","BTC","CASP9",
  "CD19","CD28","CD80","CD86",
  "CDKN1A","CDKN1B","CHUK","CREB1",
  "EGF","EGFR","ERBB2","ERBB3","ERBB4","EREG",
  "FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9",
  "FGF10","FGF16","FGF17","FGF18","FGF19","FGF20","FGF22","FGF23",
  "FGFR1","FGFR2","FGFR3","FGFR4",
  "FOXO1","FOXO3","FOXO4","FRS2",
  "HRAS","KRAS","NRAS",
  "RAF1","BRAF",
  "MAP2K1","MAP2K2",
  "MAPK1","MAPK3",
  "FOS","JUN","ELK1","SRF","MYC"
  
)

# Print total unique GH/IGF hypothesis genes once
cat("Total unique GH/IGF hypothesis genes:", length(unique(gh_symbols)), "\n")

### Map DESeq2 results to GH/IGF pathway genes -------------
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$symbol <- sub(".*\\|", "", res_df$gene_id)


gh_res <- res_df[ res_df$symbol %in% gh_symbols, ]
gh_res <- gh_res[order(gh_res$padj), ]

# Save full hypothesis table
write.csv(gh_res,
          file = "GH_IGF_hypothesis_genes_results.csv",
          row.names = FALSE)


# How many of those GH/IGF genes appear in  DESeq2 results
gh_in_data <- intersect(unique(gh_symbols), res_df$symbol)
cat("GH/IGF genes present in DESeq2 results:",
    length(gh_in_data), "\n")

### Exploratory GH/IGF subsets 

gh050 <- subset(gh_res, !is.na(padj) & padj < 0.05)
gh10 <- subset(gh_res, !is.na(padj) & padj < 0.10)
gh20 <- subset(gh_res, !is.na(padj) & padj < 0.20)

write.csv(gh050, file = "GH_IGF_genes_FDR0.05.csv", row.names = FALSE)
write.csv(gh10, file = "GH_IGF_genes_FDR0.10.csv", row.names = FALSE)
write.csv(gh20, file = "GH_IGF_genes_FDR0.20.csv", row.names = FALSE)

### Sample-to-sample distances + MA plot + PCA 

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sex)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("Fig_QC_SampleDistances.png", width = 7, height = 6, units = "in", res = 600)
sampleDistPlot <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors,
                           main = "Sample-to-sample distances")
print(sampleDistPlot)
dev.off()

png("Fig_QC_MAplot.png", width = 7, height = 6, units = "in", res = 600)
plotMA(res, main = "DESeq2: male vs female", ylim = c(-8,8))
dev.off()

png("Fig_QC_PCA_sex.png", width = 7, height = 6, units = "in", res = 600)
print(plotPCA(rld, intgroup = "sex"))
dev.off()

### Hypothesis-targeted figures
### Fig 1: QC
### Fig 2: Key GH/IGF genes
### Fig 3A: GH/IGF FDR < 0.05
### Fig 3B: GH/IGF FDR < 0.10
### GH/IGF FDR < 0.20


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
    fontsize_row = 5,
    main = " "
  )
  dev.off()
}

# Heatmap of GH/IGF genes with FDR < 0.10 

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
    fontsize_row = 6,
    main = " "
  )
  dev.off()
}


## Fig2: Box/strip plots for key GH/IGF gene

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
      title = " ",
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

### Save key DESeq2 outputs to working directory 
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

### Manuscript-ready DEG summary table ------------------------

# Helper to count up- and down-regulated genes at a given padj cutoff
deg_summary_fun <- function(res_obj, padj_cutoff) {
  sig <- res_obj[which(!is.na(res_obj$padj) & res_obj$padj < padj_cutoff), ]
  up   <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
  data.frame(
    padj_threshold       = padj_cutoff,
    n_significant_genes  = nrow(sig),
    n_upregulated_genes  = up,
    n_downregulated_genes = down
  )
}

# Build summary for padj values
deg_summary_01 <- deg_summary_fun(res, 0.10)
deg_summary_005 <- deg_summary_fun(res, 0.05)
deg_summary_020 <- deg_summary_fun(res, 0.20)

manuscript_ready_deg_table <- rbind(deg_summary_01, deg_summary_005, deg_summary_020)
row.names(manuscript_ready_deg_table) <- NULL

manuscript_ready_deg_table
write.csv(
  manuscript_ready_deg_table,
  file = "ManuscriptReady_DEG_summary_padj.csv",
  row.names = FALSE
)


#### Manuscript-ready DEG summary *restricted to GH/IGF hypothesis genes*

# Subset the full DESeq2 result to GH/IGF symbols only
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$symbol  <- sub(".*\\|", "", res_df$gene_id)

res_gh <- res_df[res_df$symbol %in% gh_symbols, ]

# Helper to count up- and down-regulated genes at a given padj cutoff
deg_summary_fun_gh <- function(res_obj, padj_cutoff) {
  sig <- res_obj[which(!is.na(res_obj$padj) & res_obj$padj < padj_cutoff), ]
  up   <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
  data.frame(
    padj_threshold        = padj_cutoff,
    n_significant_genes   = nrow(sig),
    n_upregulated_genes   = up,
    n_downregulated_genes = down
  )
}

# Build summary for padj values
deg_gh_005 <- deg_summary_fun_gh(res_gh, 0.05)
deg_gh_010 <- deg_summary_fun_gh(res_gh, 0.10)
deg_gh_020 <- deg_summary_fun_gh(res_gh, 0.20)

gh_manuscript_deg_table <- rbind(deg_gh_005, deg_gh_010, deg_gh_020)
row.names(gh_manuscript_deg_table) <- NULL
gh_manuscript_deg_table

write.csv(
  gh_manuscript_deg_table,
  file = "ManuscriptReady_GH_IGF_DEG_summary_padj.csv",
  row.names = FALSE
)

# Total unique GH/IGF genes in your hypothesis list
cat("Total unique GH/IGF hypothesis genes:",
    length(unique(gh_symbols)), "\n")
