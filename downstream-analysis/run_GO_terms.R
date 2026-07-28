
# =========================================================

# run_GO_terms.R
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

# Load libraries

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(readr)
library(ggplot2)

# Read DEG CSV (P0) 

deg <- read_csv("/path/to/WT_1_7_vs_KO_6_E16_Wilcox.csv") #change here for the specific genotype/timepoint 


# Run GO enrichment on all DEGs

- deg %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1.0)

genes_GO <- deg_sig$gene

gene_GO_df <- bitr(
  genes_GO,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
  gene = gene_GO_df$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot overview GO
dotplot(ego, showCategory = 20)

# Save GO table
write.csv(
  as.data.frame(ego),
  "GO_BP_E16_results_1_7_vs_6.csv",
  row.names = FALSE
)


# Split genes by direction


deg_WT <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC < - 1)

deg_mut <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC > 1)


# Convert gene IDs

gene_WT <- bitr(
  deg_WT$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_mut <- bitr(
  deg_mut$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)


# GO enrichment WT

ego_WT <- enrichGO(
  gene = gene_WT$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)


write.csv(
  as.data.frame(ego_WT),
  "GO_BP_E16_WT_results_1_7_vs_6.csv",
  row.names = FALSE
)

# ===============================
# GO enrichment Mutant
# ===============================

ego_mut <- enrichGO(
  gene = gene_mut$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)


write.csv(
  as.data.frame(ego_mut),
  "GO_BP_E16_mutant_results_1_7_vs_6.csv",
  row.names = FALSE
)
