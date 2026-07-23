
# =========================================================

# run_GSEA.R
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

library(dplyr)
library(readr)
library(tibble)

# Load DEG tables
deg_E16 <- read_csv("/path/to/WT_1_7_vs_KO_6_E16_Wilcox.csv")
deg_P0  <- read_csv("/path/to/WT_1_7_vs_KO_6_P0_Wilcox.csv")

# Create ranking metric

deg_E16_rank <- deg_E16 %>%
  mutate(rank_metric = sign(avg_log2FC) * -log10(p_val_adj + 1e-300)) %>%
  arrange(desc(rank_metric))

deg_P0_rank <- deg_P0 %>%
  mutate(rank_metric = sign(avg_log2FC) * -log10(p_val_adj + 1e-300)) %>%
  arrange(desc(rank_metric))


# Convert to named vector (GSEA input)


geneList_E16 <- deg_E16_rank$rank_metric
names(geneList_E16) <- deg_E16_rank$gene
geneList_E16 <- sort(geneList_E16, decreasing = TRUE)

geneList_P0 <- deg_P0_rank$rank_metric
names(geneList_P0) <- deg_P0_rank$gene
geneList_P0 <- sort(geneList_P0, decreasing = TRUE)


# Save ranking files

write.csv(
  data.frame(gene = names(geneList_E16), rank = geneList_E16),
  "GSEA_rank_E16_UBC_WT_Mutant_DEGs_1_7_vs_6_v4.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = names(geneList_P0), rank = geneList_P0),
  "GSEA_rank_P0_UBC_WT_Mutant_DEGs_1_7_vs_6_v4.csv",
  row.names = FALSE
)


#run GSEA
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)

#convert gene entries
##E16
gene_df_E16 <- bitr(
  names(geneList_E16),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

geneList_E16 <- geneList_E16[gene_df_E16$SYMBOL]
names(geneList_E16) <- gene_df_E16$ENTREZID
geneList_E16 <- sort(geneList_E16, decreasing = TRUE)

#P0
gene_df_P0 <- bitr(
  names(geneList_P0),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

geneList_P0 <- geneList_P0[gene_df_P0$SYMBOL]
names(geneList_P0) <- gene_df_P0$ENTREZID
geneList_P0 <- sort(geneList_P0, decreasing = TRUE)


#GSEA for E16

# Clean ranked gene list


deg_E16_rank <- deg_E16 %>%
  arrange(desc(avg_log2FC))

geneList_E16 <- deg_E16_rank$avg_log2FC
names(geneList_E16) <- deg_E16_rank$gene

geneList_E16 <- geneList_E16[!is.na(geneList_E16)]
geneList_E16 <- geneList_E16[!duplicated(names(geneList_E16))]
geneList_E16 <- sort(geneList_E16, decreasing = TRUE)


library(clusterProfiler)
library(org.Mm.eg.db)

# Convert SYMBOL → ENTREZ
gene_df <- bitr(
  names(geneList_E16),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Match ranking values to converted genes
geneList_E16_entrez <- geneList_E16[gene_df$SYMBOL]
names(geneList_E16_entrez) <- gene_df$ENTREZID

# Remove duplicates and NA
geneList_E16_entrez <- geneList_E16_entrez[!is.na(names(geneList_E16_entrez))]
geneList_E16_entrez <- geneList_E16_entrez[!duplicated(names(geneList_E16_entrez))]

# Sort for GSEA
geneList_E16_entrez <- sort(geneList_E16_entrez, decreasing = TRUE)

geneList_E16_entrez <- geneList_E16[gene_df$SYMBOL]

names(geneList_E16_entrez) <- gene_df$ENTREZID

geneList_E16_entrez <- geneList_E16_entrez[!duplicated(names(geneList_E16_entrez))]

geneList_E16_entrez <- sort(geneList_E16_entrez, decreasing = TRUE)

gsea_E16 <- gseGO(
  geneList = geneList_E16_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  minGSSize = 5,
  maxGSSize = 1000,
  pvalueCutoff = 0.1,
  verbose = FALSE
)

dotplot(gsea_E16, showCategory = 20)


#run GSEA for P0

# Create ranking vectors


geneList_P0 <- deg_P0$avg_log2FC
names(geneList_P0) <- deg_P0$gene
geneList_P0 <- sort(geneList_P0, decreasing = TRUE)


# Convert SYMBOL → ENTREZ

gene_df_P0 <- bitr(
  names(geneList_P0),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Match ranking values
geneList_P0_entrez <- geneList_P0[gene_df_P0$SYMBOL]
names(geneList_P0_entrez) <- gene_df_P0$ENTREZID
geneList_P0_entrez <- geneList_P0_entrez[!duplicated(names(geneList_P0_entrez))]
geneList_P0_entrez <- sort(geneList_P0_entrez, decreasing = TRUE)


# Run GSEA (GO BP)


gsea_P0 <- gseGO(
  geneList = geneList_P0_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  minGSSize = 5,
  maxGSSize = 1000,
  pvalueCutoff = 1,
  verbose = FALSE
)


# Save results tables

write.csv(
  as.data.frame(gsea_E16),
  "GSEA_BP_E16_UBC_WT_Mutant_results_1_7_vs_6_v4.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(gsea_P0),
  "GSEA_BP_P0_UBC_WT_Mutant_results_1_7_vs_6_v4.csv",
  row.names = FALSE
)

