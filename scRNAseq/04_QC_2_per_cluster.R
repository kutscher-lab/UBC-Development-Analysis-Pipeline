# =========================
# Libraries
# =========================
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
library(patchwork)
library(cowplot)
library(viridis)

# =========================
# Global settings
# =========================
set.seed(1234)

outdir <- "/path/to/outdir/to/QC_v2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS("/path/to/merged/object")

magma_rev <- rev(viridis::magma(256))

# =========================
# Helper functions
# =========================

save_plot_pdf <- function(plot_obj,
                          filename,
                          width = 8,
                          height = 6) {

  ggsave(
    filename = file.path(outdir, filename),
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height,
    useDingbats = FALSE
  )
}


make_qc_violin <- function(obj, group_var = "seurat_clusters") {

  p <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
    group.by = group_var,
    ncol = 2,
    pt.size = 0
  )

  return(p)
}


make_umap <- function(obj) {

  p <- DimPlot(
    obj,
    group.by = "seurat_clusters",
    label = TRUE
  )

  return(p)
}


make_dotplot <- function(obj, markers) {

  p <- DotPlot(
    obj,
    features = markers
  ) +
    RotatedAxis() +
    theme_bw() +
    labs(x = "Markers", y = "Cell types") +
    scale_color_gradientn(colors = magma_rev) +
    scale_size(
      range = c(0, 10),
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  return(p)
}


process_seurat <- function(obj,
                           dims_use = 1:20,
                           resolution = 0.5) {

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)

  obj <- RunUMAP(
    obj,
    dims = dims_use
  )

  obj <- FindNeighbors(
    obj,
    dims = dims_use
  )

  obj <- FindClusters(
    obj,
    resolution = resolution
  )

  return(obj)
}

# =========================
# Load object
# =========================

combined <- readRDS(
  "/path/to/merged/object/pre/annotation.RDS"
)

# =========================
# Initial QC
# =========================

p_umap_initial <- make_umap(combined)

save_plot_pdf(
  p_umap_initial,
  "initial_umap_clusters.pdf"
)

p_qc_initial <- make_qc_violin(combined)

save_plot_pdf(
  p_qc_initial,
  "initial_qc_violin.pdf",
  width = 12,
  height = 8
)

# =========================
# QC filtering
# =========================

combined_subset <- subset(
  combined,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 7000
) #arrange the filters based on the distriutions in each cluster

p_qc_filtered <- make_qc_violin(combined_subset)

save_plot_pdf(
  p_qc_filtered,
  "filtered_qc_violin.pdf",
  width = 12,
  height = 8
)

# =========================
# Remove low-quality clusters
# =========================
# Remove clusters 12 and 16
# 12 = low quality
# 16 = very low cell number / contamination

combined_sub_2 <- subset(
  combined_subset,
  idents = c(12, 16),
  invert = TRUE
)

# =========================
# PCA inspection
# =========================

p_elbow_50 <- ElbowPlot(
  combined_sub_2,
  ndims = 50
)

save_plot_pdf(
  p_elbow_50,
  "elbowplot_50pcs.pdf",
  width = 8,
  height = 6
)

# =========================
# Reprocess object
# =========================

combined_sub_2 <- process_seurat(
  combined_sub_2,
  dims_use = 1:20,
  resolution = 0.5
)

# =========================
# QC and clustering inspection
# =========================

p_umap_sub2 <- make_umap(combined_sub_2)

save_plot_pdf(
  p_umap_sub2,
  "combined_sub_2_umap.pdf"
)

p_qc_sub2 <- make_qc_violin(combined_sub_2)

save_plot_pdf(
  p_qc_sub2,
  "combined_sub_2_qc_violin.pdf",
  width = 12,
  height = 8
)

table(combined_sub_2$seurat_clusters)

# =========================
# Marker inspection
# =========================

markers <- c(
  "Mki67","Top2a","Sox2","Pcna","Atoh1","Otx2","Pax6",
  "Eomes","Lmx1a","Barhl1","Barhl2",
  "Calb2","Mgat5","Trpc3","Grm1","Grm8","En2",
  "Hcrtr2","Cxcr4",
  "Synpr","Neurod1","Zic1","Zic2","Kcnip4","Dcx",
  "Tlx3","En1","Pde1c","Cntn1","Pcsk9","Adarb2",
  "Gabra6","Rbfox3",
  "Gad1","Gad2","Hoxb2","Hoxb3",
  "Slc17a6","Slc17a7",
  "Aldh1l1","Gfap","Slc5a7","Gabrg3",
  "Scg2","Pax5","Slit2",
  "Wls","Olig2","Ptf1a","Gsx1",
  "Meis2","Kirrel2","Ttr"
)

Idents(combined_sub_2) <- "seurat_clusters"

combined_sub_2$seurat_clusters <- factor(
  combined_sub_2$seurat_clusters
)

p_dot_sub2 <- make_dotplot(
  combined_sub_2,
  markers
)

save_plot_pdf(
  p_dot_sub2,
  "combined_sub_2_marker_dotplot.pdf",
  width = 18,
  height = 8
)

table(
  combined_sub_2$seurat_clusters,
  combined_sub_2$sample
)

#before deciding, check the clusters per sample

clusters_check <- subset(
  combined_sub_2,
  idents = c(12, 15, 16)
)

# -------------------------
# QC violin plots by sample
# -------------------------

p_qc_sample <- VlnPlot(
  clusters_check,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mito"
  ),
  group.by = "seurat_clusters",
  split.by = "sample",
  pt.size = 0,
  ncol = 3
) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

save_plot_pdf(
  p_qc_sample,
  "clusters_12_15_16_qc_by_sample.pdf",
  width = 20,
  height = 10
)

p_qc_sample
# =========================
# Remove problematic clusters
# =========================
# Remove:
# 12 = low-quality GC-like cells
# 15 = WT-specific overlap cluster inside the main clusters and low QC
# 16 = Hox-positive contamination but QC seems fine (dont remove this)

combined_sub_3 <- subset(
  combined_sub_2,
  idents = c(12, 15),
  invert = TRUE
)

combined_sub_3 <- process_seurat(
  combined_sub_3,
  dims_use = 1:20,
  resolution = 0.5
)

# =========================
# Inspect reclustered object
# =========================

p_umap_sub3 <- make_umap(combined_sub_3)

save_plot_pdf(
  p_umap_sub3,
  "combined_sub_3_umap.pdf"
)

p_qc_sub3 <- make_qc_violin(combined_sub_3)

save_plot_pdf(
  p_qc_sub3,
  "combined_sub_3_qc_violin.pdf",
  width = 12,
  height = 8
)

p_dot_sub3 <- make_dotplot(
  combined_sub_3,
  markers
)

save_plot_pdf(
  p_dot_sub3,
  "combined_sub_3_marker_dotplot.pdf",
  width = 18,
  height = 8
)

table(
  combined_sub_3$seurat_clusters,
  combined_sub_3$sample
)

# =========================
# Optional higher-resolution clustering
# =========================

combined_sub_4 <- FindClusters(
  combined_sub_3,
  resolution = 0.6
)

p_umap_sub_4 <- make_umap(combined_sub_4)

save_plot_pdf(
  p_umap_sub_4,
  "combined_sub_4_umap_res0.6.pdf"
)

p_qc_sub4 <- make_qc_violin(combined_sub_4)

save_plot_pdf(
  p_qc_sub4,
  "combined_sub_4_qc_violin.pdf",
  width = 12,
  height = 8
)

p_dot_sub4 <- make_dotplot(
  combined_sub_4,
  markers
)

save_plot_pdf(
  p_dot_sub4,
  "combined_sub_4_marker_dotplot.pdf",
  width = 18,
  height = 8
)

table(
  combined_sub_4$seurat_clusters,
  combined_sub_4$sample
)


# =========================
# Higher-resolution clustering comparison
# =========================

resolutions <- c(0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {

  # -------------------------
  # Reclustering
  # -------------------------

  obj_res <- FindClusters(
    combined_sub_3,
    resolution = res
  )

  # -------------------------
  # UMAP
  # -------------------------

  p_umap <- make_umap(obj_res) +
    ggtitle(paste0("Resolution ", res)) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  save_plot_pdf(
    p_umap,
    paste0("combined_sub_resolution_", res, "_umap.pdf"),
    width = 10,
    height = 8
  )

  print(p_umap)

  # -------------------------
  # QC violin plots
  # -------------------------

  p_qc <- VlnPlot(
    obj_res,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mito"
    ),
    group.by = "seurat_clusters",
    split.by = "sample",
    pt.size = 0,
    ncol = 3
  ) +
    ggtitle(paste0("Resolution ", res)) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  save_plot_pdf(
    p_qc,
    paste0("combined_sub_resolution_", res, "_qc_violin.pdf"),
    width = 20,
    height = 10
  )

  print(p_qc)

  # -------------------------
  # Marker dotplot
  # -------------------------

  p_dot <- make_dotplot(
    obj_res,
    markers
  ) +
    ggtitle(paste0("Resolution ", res))

  save_plot_pdf(
    p_dot,
    paste0("combined_sub_resolution_", res, "_marker_dotplot.pdf"),
    width = 18,
    height = 8
  )

  print(p_dot)

  # -------------------------
  # Sample composition
  # -------------------------

  cat("\n=========================\n")
  cat("Resolution:", res, "\n")
  cat("=========================\n")

  print(
    table(
      obj_res$seurat_clusters,
      obj_res$sample
    )
  )

  print(
    prop.table(
      table(
        obj_res$seurat_clusters,
        obj_res$sample
      ),
      margin = 2
    )
  )
}


# =========================
# QC violin plots for multiple resolutions
# =========================

resolutions <- c(0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {

  # -------------------------
  # Reclustering
  # -------------------------

  obj_res <- FindClusters(
    combined_sub_3,
    resolution = res
  )

  # -------------------------
  # QC violin plots by cluster
  # -------------------------

  p_qc_cluster <- VlnPlot(
    obj_res,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mito"
    ),
    group.by = "seurat_clusters",
    pt.size = 0,
    ncol = 3
  ) +
    ggtitle(paste0("QC by cluster - Resolution ", res)) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  save_plot_pdf(
    p_qc_cluster,
    paste0("combined_sub_resolution_", res, "_qc_by_cluster.pdf"),
    width = 18,
    height = 10
  )

  print(p_qc_cluster)

  # -------------------------
  # QC violin plots split by sample
  # -------------------------

  p_qc_sample <- VlnPlot(
    obj_res,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mito"
    ),
    group.by = "seurat_clusters",
    split.by = "sample",
    pt.size = 0,
    ncol = 3
  ) +
    ggtitle(paste0("QC by sample - Resolution ", res)) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  save_plot_pdf(
    p_qc_sample,
    paste0("combined_sub_resolution_", res, "_qc_by_sample.pdf"),
    width = 20,
    height = 10
  )

  print(p_qc_sample)
}


# Reclustering at resolution 1.0
obj_res_1 <- FindClusters(
  combined_sub_3,
  resolution = 1.0
)

# Number of cells per cluster
table(obj_res_1$seurat_clusters)

table(
  obj_res_1$seurat_clusters,
  obj_res_1$sample
)


# =========================
# Save cluster cell counts for multiple resolutions
# =========================

resolutions <- c(0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {

  # -------------------------
  # Reclustering
  # -------------------------

  obj_res <- FindClusters(
    combined_sub_3,
    resolution = res
  )

  # -------------------------
  # Cells per cluster
  # -------------------------

  cluster_counts <- as.data.frame(
    table(obj_res$seurat_clusters)
  )

  colnames(cluster_counts) <- c(
    "cluster",
    "n_cells"
  )

  # Save overall counts
  write.csv(
    cluster_counts,
    file = file.path(
      outdir,
      paste0(
        "resolution_",
        res,
        "_cluster_cell_counts.csv"
      )
    ),
    row.names = FALSE
  )

  # -------------------------
  # Cells per cluster by sample
  # -------------------------

  cluster_sample_counts <- as.data.frame(
    table(
      obj_res$seurat_clusters,
      obj_res$sample
    )
  )

  colnames(cluster_sample_counts) <- c(
    "cluster",
    "sample",
    "n_cells"
  )

  # Save sample-specific counts
  write.csv(
    cluster_sample_counts,
    file = file.path(
      outdir,
      paste0(
        "resolution_",
        res,
        "_cluster_sample_cell_counts.csv"
      )
    ),
    row.names = FALSE
  )

  # -------------------------
  # Print to console
  # -------------------------

  cat("\n=========================\n")
  cat("Resolution:", res, "\n")
  cat("=========================\n")

  print(cluster_counts)

  print(cluster_sample_counts)
}

