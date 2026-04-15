### =========================================================
#Course: Functional Genomics (BIOL 6850)
# Author: Beverly Thomas
#### Functional categories using PANTHER GO BP overrepresentation table
# Purpose:
#   Start from a PANTHER GO biological process overrepresentation .csv,
#   map enriched GO BP terms into broad functional categories,
#   and visualize counts per category with a barplot.
#
# Input files expected in working directory:
#   1) PANTHER_GO_BP_overrep.xml
#   2) panther_term_to_category.csv
#
# panther_term_to_category.csv must contain:
#   PantherTermName,BroadCategory
#
# Example:
#   PantherTermName,BroadCategory
#   glycolipid translocation,Metabolic
#   regulation of localization,Synaptic / neuronal
#
# Output files:
#   PANTHER_GO_BP_terms_with_categories.csv
#   Functional_category_counts_panther_terms.csv
#   Functional_categories_panther_terms_barplot.pdf
#
# Output files:
#   - DEG_functional_categories_panther_fgsea.csv
#   - Functional_category_counts_panther_fgsea.csv
#   - Functional_categories_panther_fgsea_barplot.pdf
#
### =========================================================
### 1. Read PANTHER GO BP overrepresentation file
### =========================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
library(xml2)
library(dplyr)})

# Read UTF-16 XML as text, then parse
xml_lines <- readLines("PANTHER_GO_BP_overrep.xml", encoding = "UTF-16", warn = FALSE)

# Collapse back to one XML string
xml_txt <- paste(xml_lines, collapse = "\n")

# Parse from text instead of file path
doc <- read_xml(xml_txt)

# Pull all result nodes
res_nodes <- xml_find_all(doc, ".//result")

# Build gene x term table

term_gene_df <- bind_rows(lapply(res_nodes, function(res) {
  x <- as_list(res)
  
  term_label <- x$term$label[[1]]
  fdr_text   <- x$input_list$fdr[[1]]
  pm_text    <- x$input_list$plus_minus[[1]]
  
  mapped_ids <- x$input_list$mapped_id_list
  if (is.null(mapped_ids)) return(NULL)
  
  mapped_ids <- unlist(mapped_ids, use.names = FALSE)
  if (length(mapped_ids) == 0) return(NULL)
  
  data.frame(
    GeneSymbol = mapped_ids,
    PantherTermName = term_label,
    plus_minus = pm_text,
    fdr = as.numeric(fdr_text),
    stringsAsFactors = FALSE
  )
}))
# Keep only significantly enriched (+) terms
term_gene_sig <- term_gene_df %>%
  filter(!is.na(fdr) & fdr < 0.20 & plus_minus == "+")

write.csv(
  term_gene_sig,
  "PANTHER_gene_to_term_map_no_category.csv",
  row.names = FALSE
)

print(dim(term_gene_df))
print(dim(term_gene_sig))
print(head(term_gene_sig))

# Save a clean file for manual category labeling
unique_terms <- data.frame(
  PantherTermName = sort(unique(term_gene_sig$PantherTermName)),
  stringsAsFactors = FALSE
)

write.csv(
  unique_terms,
  "panther_unique_terms_to_label.csv",
  row.names = FALSE
)

print(dim(unique_terms))
print(head(unique_terms, 20))

cat("Done.\n")
cat("Rows in full gene-term table:", nrow(term_gene_df), "\n")
cat("Rows in significant positive gene-term table:", nrow(term_gene_sig), "\n")
head(term_gene_sig)


### =========================================================
### 2. Automatically assign broad functional categories
###    from PANTHER GO term names using curated keyword bins
### =========================================================

term_gene_sig$term_lower <- tolower(term_gene_sig$PantherTermName)
term_gene_sig$BroadCategory <- NA_character_

## Cytoskeletal
cyto_pat <- "actin|microtubule|cytoskeleton|cell-cell adhesion|cell adhesion|junction|extracellular matrix|extracellular structure|epithelium|projection morphogenesis|projection organization|plasma membrane bounded cell projection"
term_gene_sig$BroadCategory[
  grepl(cyto_pat, term_gene_sig$term_lower)
] <- "Cytoskeletal"

## Synaptic 
neuro_pat <- "axon|dendrit|synap|neur|hippocampus|cortex|forebrain|telencephalon|pallium|subpallium|olfactory bulb|nervous system|memory|limbic"
term_gene_sig$BroadCategory[
  is.na(term_gene_sig$BroadCategory) &
    grepl(neuro_pat, term_gene_sig$term_lower)
] <- "Synaptic "

## Metabolic
met_pat <- "metabolic process|catabolic process|biosynthetic process|homeostasis|amino acid transport|amino acid transmembrane transport|glutamate transmembrane transport|carboxylic acid transmembrane transport|d-xylulose 5-phosphate|guanine metabolic process|nitrogen compound transport|organic cation transport|myo-inositol transport|nad transmembrane transport|nad transport|sulfur amino acid transport|sulfur compound transport|proteoglycan biosynthetic process"
term_gene_sig$BroadCategory[
  is.na(term_gene_sig$BroadCategory) &
    grepl(met_pat, term_gene_sig$term_lower)
] <- "Metabolic"

##  Hormonal 

horm_pat <- "insulin|igf|growth factor|egf|tgf|hormone|endocrine||serine/threonine kinase signaling|growth|differentiation|development|osteoblast|hematopoietic progenitor"
term_gene_sig$BroadCategory[
  is.na(term_gene_sig$BroadCategory) &
    grepl(horm_pat, term_gene_sig$term_lower)
] <- "Hormonal"


## Immune
immune_pat <- "immune|inflammat|viral|virus|interferon|chemokine|leukocyte|hematopoietic|defense response, oxidat"
term_gene_sig$BroadCategory[
  is.na(term_gene_sig$BroadCategory) &
    grepl(immune_pat, term_gene_sig$term_lower)
] <- "Immune"

cat("\nBroadCategory assignment summary (after keyword rules):\n")
print(table(term_gene_sig$BroadCategory, useNA = "ifany"))

panther_map <- term_gene_sig %>%
  filter(!is.na(BroadCategory) & BroadCategory != "")

cat_df <- panther_map %>%
  distinct(GeneSymbol, BroadCategory) %>%
  count(BroadCategory, sort = TRUE) %>%
  rename(category = BroadCategory, count = n)

### =========================================================
### 3. Make barplot
### =========================================================

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

ggsave(
  "Functional_categories_panther_barplot.pdf",
  plot = p_cat,
  width = 8,
  height = 5
)

### =========================================================
### 4. Optional: inspect genes and GO terms by category
### =========================================================

cat("\nGO terms by category:\n")
print(split(unique(panther_map$PantherTermName), panther_map$BroadCategory))

cat("\nGenes by category:\n")
print(split(unique(panther_map$GeneSymbol), panther_map$BroadCategory))
