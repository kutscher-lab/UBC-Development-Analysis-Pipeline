# =========================
# Libraries
# =========================
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
library(patchwork)
library(cowplot)

library(Seurat)
#Read the RDS
combined <- readRDS("/path/to/filtered/clean/merged/object/with/resolution/07")
p_umap_clusters <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE)


#visualize and save the umap


p_umap_geno <- DimPlot(
  combined,
  group.by = "genotype",
  label = FALSE,
  repel = TRUE
)

p_umap_time <- DimPlot(
  combined,
  group.by = "timepoint",
  label = FALSE,
  repel = TRUE
)


# =========================
# Common theme (NO legend)
# =========================
umap_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# =========================
# CELLTYPE (D3 palette, 17 clusters)
# =========================
n_cell <- length(unique(combined$seurat_clusters))

cell_colors <- colorRampPalette(pal_d3()(10))(n_cell)

p_celltype <- DimPlot(
  combined,
  group.by = "seurat_clusters",
  pt.size = 0.4,
  label = FALSE
) +
  scale_color_manual(values = cell_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_celltype <- get_legend(
  DimPlot(combined, group.by = "seurat_clusters", pt.size = 0.4) +
    scale_color_manual(values = cell_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8)
    )
)

# =========================
# GENOTYPE (3 colors)
# =========================
geno_colors <- c("#eb681cff", "#4c99d9ff", "#d71070ff")

p_genotype <- DimPlot(
  combined,
  group.by = "genotype",
  pt.size = 0.4
) +
  scale_color_manual(values = geno_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_genotype <- get_legend(
  DimPlot(combined, group.by = "genotype", pt.size = 0.6) +
    scale_color_manual(values = geno_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 9)
    )
)

# =========================
# TIMEPOINT (2 colors)
# =========================
time_colors <- c("#1ceb42ff", "#e331aaff")

p_timepoint <- DimPlot(
  combined,
  group.by = "timepoint",
  pt.size = 0.4
) +
  scale_color_manual(values = time_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_timepoint <- get_legend(
  DimPlot(combined, group.by = "timepoint", pt.size = 0.6) +
    scale_color_manual(values = time_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 10)
    )
)

# =========================
# SAVE UMAPS (NO legends)
# =========================
ggsave("/path/to/Merged_Object_UMAP_celltype_unfiltered_v4_07.pdf", p_celltype, width = 10, height = 8, device = cairo_pdf)
ggsave("/path/to/Merged_Object_UMAP_genotype_unfiltered_v4_07.pdf", p_genotype, width = 10, height = 8, device = cairo_pdf)
ggsave("/path/to/Merged_Object_UMAP_timepoint_unfiltered_v4_07.pdf", p_timepoint, width = 10, height = 8, device = cairo_pdf)


# =========================
# COMBINE UMAP + LEGEND
# =========================

# Celltype (big legend → give more space)
p_celltype_with_legend <- plot_grid(
  p_celltype,
  legend_celltype,
  ncol = 2,
  rel_widths = c(1, 0.6)
)

# Genotype
p_genotype_with_legend <- plot_grid(
  p_genotype,
  legend_genotype,
  ncol = 2,
  rel_widths = c(1, 0.4)
)

# Timepoint
p_timepoint_with_legend <- plot_grid(
  p_timepoint,
  legend_timepoint,
  ncol = 2,
  rel_widths = c(1, 0.3)
)


# =========================
# SAVE WITH LEGENDS
# =========================

ggsave(
  "/path/to/Merged_Object_UMAP_celltype_with_legend_unfiltered_v4_07.pdf",
  p_celltype_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_genotype_with_legend_unfiltered_v4_07.pdf",
  p_genotype_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_timepoint_with_legend_unfiltered_v4_07.pdf",
  p_timepoint_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)


#remove off-target cell types for the visualization and analyses (for separate visualization only)

clusters_to_remove <- c(14, 13, 17, 18)

combined_filt <- subset(
  combined,
  subset = !(seurat_clusters %in% clusters_to_remove)
)


#check if removed

table(combined_filt$seurat_clusters)


# =========================
# Common theme (NO legend)
# =========================
umap_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# =========================
# CELLTYPE (D3 palette)
# =========================
n_cell <- length(unique(combined_filt$seurat_clusters))

cell_colors <- colorRampPalette(pal_d3()(10))(n_cell)

p_celltype <- DimPlot(
  combined_filt,
  group.by = "seurat_clusters",
  pt.size = 0.4,
  label = FALSE
) +
  scale_color_manual(values = cell_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_celltype <- get_legend(
  DimPlot(combined_filt, group.by = "seurat_clusters", pt.size = 0.4) +
    scale_color_manual(values = cell_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8)
    )
)


#keep the colors of unfiltered umap cell types
# =========================
# DEFINE COLORS ON FULL DATA (IMPORTANT)
# =========================
library(ggsci)
library(scales)

all_celltypes <- sort(unique(combined$seurat_clusters))

cell_colors <- setNames(
  colorRampPalette(pal_d3()(10))(length(all_celltypes)),
  all_celltypes
)

# check
table(combined_filt$seurat_clusters)


# =========================
# KEEP ONLY USED COLORS (BUT SAME MAPPING)
# =========================
cell_colors_filt <- cell_colors[
  unique(combined_filt$seurat_clusters)
]

# =========================
# PLOT
# =========================
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)

umap_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

p_celltype <- DimPlot(
  combined_filt,
  group.by = "seurat_clusters",
  pt.size = 0.4,
  label = FALSE
) +
  scale_color_manual(values = cell_colors_filt) +
  umap_theme +
  coord_fixed()

# =========================
# LEGEND
# =========================
legend_celltype <- get_legend(
  DimPlot(combined_filt, group.by = "seurat_clusters", pt.size = 0.4) +
    scale_color_manual(values = cell_colors_filt) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8)
    )
)


# =========================
# GENOTYPE (3 colors)
# =========================
geno_colors <- c("#eb681cff", "#4c99d9ff", "#d71070ff")

p_genotype <- DimPlot(
  combined_filt,
  group.by = "genotype",
  pt.size = 0.4
) +
  scale_color_manual(values = geno_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_genotype <- get_legend(
  DimPlot(combined_filt, group.by = "genotype", pt.size = 0.6) +
    scale_color_manual(values = geno_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 9)
    )
)

# =========================
# TIMEPOINT (2 colors)
# =========================
time_colors <- c("#1ceb42ff", "#e331aaff")

p_timepoint <- DimPlot(
  combined_filt,
  group.by = "timepoint",
  pt.size = 0.4
) +
  scale_color_manual(values = time_colors) +
  umap_theme +
  coord_fixed()

# --- legend (separate)
legend_timepoint <- get_legend(
  DimPlot(combined_filt, group.by = "timepoint", pt.size = 0.6) +
    scale_color_manual(values = time_colors) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 10)
    )
)

# =========================
# SAVE UMAPS (NO legends)
# =========================
ggsave(
  "/path/to/Merged_Object_UMAP_celltype_filtered_v4_07.pdf",
  p_celltype, width = 10, height = 8, device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_genotype_filtered_v4_07.pdf",
  p_genotype, width = 10, height = 8, device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_timepoint_filtered_v4_07.pdf",
  p_timepoint, width = 10, height = 8, device = cairo_pdf
)

# =========================
# COMBINE UMAP + LEGEND
# =========================

p_celltype_with_legend <- plot_grid(
  p_celltype,
  legend_celltype,
  ncol = 2,
  rel_widths = c(1, 0.6)
)

p_genotype_with_legend <- plot_grid(
  p_genotype,
  legend_genotype,
  ncol = 2,
  rel_widths = c(1, 0.4)
)

p_timepoint_with_legend <- plot_grid(
  p_timepoint,
  legend_timepoint,
  ncol = 2,
  rel_widths = c(1, 0.3)
)

# =========================
# SAVE WITH LEGENDS
# =========================
ggsave(
  "/path/to/Merged_Object_UMAP_celltype_with_legend_filtered_v4_07.pdf",
  p_celltype_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_genotype_with_legend_filtered_v4_07.pdf",
  p_genotype_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)

ggsave(
  "/path/to/Merged_Object_UMAP_timepoint_with_legend_filtered_v4_07.pdf",
  p_timepoint_with_legend,
  width = 10,
  height = 8,
  device = cairo_pdf
)

