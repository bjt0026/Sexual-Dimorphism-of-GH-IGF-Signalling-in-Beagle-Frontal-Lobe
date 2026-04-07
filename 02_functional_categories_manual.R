# Functional categorization of DEGs using manual gene sets
# Course: Functional Genomics (BIOL 6850)
# Author: Beverly Thomas
# Input: DESeq2_results_male_vs_female_ordered.csv (from DESeq2 script)
# Output:
#   - DEG_functional_categories_manual.csv
#   - Functional_category_counts_manual.csv
#   - Functional_categories_manual_barplot.pdf
# Description:
#   Filters DEGs at padj < 0.05, assigns each to a broad functional
#   category using curated gene lists and visualizes counts per category.
### Categories:
###   1) synaptic/neuronal
###   2) metabolic/mitochondrial
###   3) immune/inflammatory & glial
###   4) hormonal/signaling (incl. GH–IGF)
###   5) other/structural

### =========================================================

library(ggplot2)

### 1. Read DESeq2 results --------------------------------------------

res_df <- read.csv("DESeq2_results_male_vs_female_ordered.csv",
                   row.names = 1,
                   check.names = FALSE)

res_df$gene_id <- rownames(res_df)
res_df$symbol  <- sub(".*\\|", "", res_df$gene_id)

dim(res_df)
head(res_df)

### 2. Define DEG threshold -------------------------------------------

deg <- subset(res_df, !is.na(padj) & padj < 0.05)
deg <- deg[order(deg$padj), ]

cat("Number of DEGs (padj < 0.05):", nrow(deg), "\n")

### 3. Manual broad functional category gene sets ---------------------
### These are the literature-informed seed lists and can be edited.

# synaptic / neuronal signalling
synaptic_genes <- c(
  "SNAP25","SYN1","SYN2","SYP","SYT1","SYT2","VAMP2","STX1A","STXBP1",
  "RAB3A","DLG4","SHANK1","SHANK2","SHANK3","NRXN1","NRXN2","NLGN1",
  "NLGN2","NLGN3","GRIN1","GRIN2A","GRIN2B","GRIA1","GRIA2","GABRA1",
  "GABRB2","GABRG2","CAMK2A","CAMK2B","SLC17A7","SLC32A1","ELAVL3",
  "RBFOX3","MAP2","TUBB3","DCX","NEFL","NEFM","NEFH","BDNF","ARC","FOS","JUN"
)

# metabolic / mitochondrial / oxidative stress
metabolic_genes <- c(
  "NDUFA1","NDUFA2","NDUFA9","NDUFS1","NDUFS2","NDUFS3","NDUFV1","SDHA",
  "SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","COX6A1","ATP5F1A","ATP5F1B",
  "ATP5MC1","ATP5MC2","ATP5PF","PDHA1","PDHB","CS","IDH3A","IDH3B","MDH2",
  "SLC25A4","SLC25A5","SLC25A3","TFAM","TOMM20","TIMM23","CPT1A","CPT2",
  "ACADM","ACADL","HSPD1","HSPE1","SOD1","SOD2","GPX1","CAT","PRDX3","PRDX5"
)

# immune / inflammatory / glial
immune_genes <- c(
  "AIF1","ITGAM","CSF1R","CX3CR1","TREM2","TYROBP","C1QA","C1QB","C1QC",
  "C3","C4A","C4B","TLR2","TLR4","IL1B","IL6","TNF","CCL2","CCL3","CCL4",
  "CXCL10","NFKB1","RELA","PTGS2","NOS2","GFAP","VIM","S100B","ALDH1L1",
  "AQP4","OLIG1","OLIG2","MBP","MOG","PLP1","MOBP","HLA-DRA","HLA-DRB1",
  "HLA-A","HLA-B","HLA-C","STAT1","IRF1","IRF7"
)

# hormonal / signalling / endocrine (includes GH-IGF-insulin-steroid pathways)
hormonal_genes <- c(
  "GH1","GH2","GHR","GHRHR","GHRH","IGF1","IGF2","IGF1R","IGF2R",
  "IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IGFALS",
  "INS","INSR","IRS1","IRS2","JAK2","STAT3","STAT5A","STAT5B",
  "PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","AKT1","AKT2","AKT3",
  "TSC1","TSC2","MTOR","RPTOR","RICTOR","MAPK1","MAPK3","MAP2K1","MAP2K2",
  "RAF1","BRAF","ELK1","SOCS1","SOCS2","SOCS3","CISH","PTPN1","PTPN2",
  "ESR1","ESR2","AR","PGR","NR3C1","NR3C2","THRA","THRB","CRHR1","CRHR2",
  "GNRHR","KISS1R","LHCGR","FSHR"
)

# structural / cytoskeleton 
structural_genes <- c(
  "ACTB","ACTG1","TUBA1A","TUBA1B","TUBB","TUBB2A","TUBB2B","TUBB3",
  "MAP1A","MAP1B","MAP2","NEFL","NEFM","NEFH","VCL","PXN","TLN1","DST",
  "FLNA","FLNB","SPTAN1","SPTBN1","SPTBN2","VIM","GFAP","COL1A1","COL1A2",
  "COL4A1","COL4A2","LAMA1","LAMB1","LAMC1","FN1","TNC","NCAM1","L1CAM",
  "ITGB1","ITGA6","CDH2","CDH11","CTNNA1","CTNNB1"
)

### 4. Assign each DEG to one main category ---------------------------
### NOTE : Priority order matters; later assignments overwrite earlier ones.

deg$category <- "other/unclassified"

deg$category[deg$symbol %in% structural_genes] <- "Cytoskeletal"
deg$category[deg$symbol %in% synaptic_genes]   <- "Synaptic / neuronal"
deg$category[deg$symbol %in% metabolic_genes]  <- "Metabolic"
deg$category[deg$symbol %in% immune_genes]     <- "Immune / inflammatory"
deg$category[deg$symbol %in% hormonal_genes]   <- "Hormonal / GH-IGF signaling"

### 5. Summarize category counts --------------------------------------

cat_tab <- table(deg$category)

### OPTIONAL ### -------------------------------
# remove other/unclassified completely from graph if it overpowers the y axis 
cat_tab <- cat_tab[names(cat_tab) != "other/unclassified"]

# sort after removing it
cat_tab <- sort(cat_tab, decreasing = TRUE)
cat_tab

cat_df <- as.data.frame(cat_tab)
colnames(cat_df) <- c("category", "count")

write.csv(deg,
          file = "DEG_functional_categories_manual.csv",
          row.names = FALSE)

write.csv(cat_df,
          file = "Functional_category_counts_manual.csv",
          row.names = FALSE)
### 6. Make barplot ---------------------------------------------------

p_cat <- ggplot(cat_df, aes(x = reorder(category, -count), y = count)) +
  geom_col(fill = "purple") +
  theme_bw(base_size = 12) +
  labs(
    title = " ",
    x = "Functional category",
    y = "Number of DEGs"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave("Functional_categories_manual_barplot.pdf",
       plot = p_cat, width = 8, height = 5)

### 7. Print genes by category for easy review ------------------------

split(deg$symbol, deg$category)
