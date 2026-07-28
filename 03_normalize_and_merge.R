#load the seurat objects (with singlets)
library(Seurat)

indir <- "/path/to/doublet/removed/filtered/objects"

files <- list.files(indir, pattern = "singlets_clean.RDS$", full.names = TRUE)

names(files) <- gsub("_singlets_clean.RDS", "", basename(files))

objs <- lapply(files, readRDS)

names(objs)


#rename the barcodes to assign to the sample name

objs <- Map(function(obj, nm) {
  RenameCells(obj, add.cell.id = nm)
}, objs, names(objs))

#check them
sapply(objs, function(x) head(colnames(x), 1))

#add clean meta data
#define the helper function
add_metadata <- function(obj, timepoint, genotype) {
  obj$timepoint <- factor(timepoint, levels = c("E16", "P0"))
  obj$genotype  <- genotype
  obj$condition <- genotype
  obj
}

#apply the function to all objects (for all the filtered datasets)
objs[["E16_ires_gfp"]] <- add_metadata(objs[["E16_ires_gfp"]], "E16", "ires_gfp")
objs[["E16_Tbr2_flx"]] <- add_metadata(objs[["E16_Tbr2_flx"]], "E16", "Tbr2_flx")
objs[["E16_Tbr2_WT"]]  <- add_metadata(objs[["E16_Tbr2_WT"]],  "E16", "Tbr2_WT")

objs[["P0_ires_gfp"]]  <- add_metadata(objs[["P0_ires_gfp"]],  "P0",  "ires_gfp")
objs[["P0_Tbr2_flx"]]  <- add_metadata(objs[["P0_Tbr2_flx"]],  "P0",  "Tbr2_flx")
objs[["P0_Tbr2_WT"]]   <- add_metadata(objs[["P0_Tbr2_WT"]],   "P0",  "Tbr2_WT")

#check
table(objs[["E16_ires_gfp"]]$timepoint)
table(objs[["E16_ires_gfp"]]$genotype)

table(objs[["P0_Tbr2_WT"]]$timepoint)
table(objs[["P0_Tbr2_WT"]]$genotype)

#apply globally
sapply(objs, function(x) unique(x$timepoint))
sapply(objs, function(x) unique(x$genotype))

#save the new metadata
outdir2 <- "path/to/new/metadata/for/all"
dir.create(outdir2, showWarnings = FALSE)

for (nm in names(objs)) {
  saveRDS(
    objs[[nm]],
    file = file.path(outdir2, paste0(nm, "_singlets_barcodes_metadata.RDS"))
  )
}

#Merge all objects
combined <- merge(
  x = objs[[1]],
  y = objs[-1]
)

combined

#Normalize and Find Variable Features
combined <- NormalizeData(combined, normalization.method = "LogNormalize")
combined <- FindVariableFeatures(combined, nfeatures = 3000)

#check
VariableFeaturePlot(combined)

#save the variable feature plot
library(ggplot2)
plotdir <- "path/to/plot/outdir"
dir.create(plotdir, showWarnings = FALSE)

p1 <- VariableFeaturePlot(combined)

ggsave(
  filename = file.path(plotdir, "Merged_Object_VariableFeaturePlot.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

#Scale the data (regress QC only)
combined <- ScaleData(
  combined,
  vars.to.regress = c("nCount_RNA", "percent.mito")
)

#check and save the QC Violin
p_qc <- VlnPlot(combined, features = c("nCount_RNA", "percent.mito"), ncol = 2)

ggsave(
  filename = file.path(plotdir, "Merged_Object_QC_violin_after_scaling.png"),
  plot = p_qc,
  width = 8,
  height = 4,
  dpi = 300
)

#run PCA
combined <- RunPCA(combined, npcs = 50)

#plot and save the ElbowPlot
p_elbow <- ElbowPlot(combined)

ggsave(
  filename = file.path(plotdir, "Merged_Object_PCA_ElbowPlot.png"),
  plot = p_elbow,
  width = 6,
  height = 5,
  dpi = 300
)

#Find Neighbors and run UMAP

combined <- FindNeighbors(combined, dims = 1:30)
combined <- RunUMAP(combined, dims = 1:30)

#UMAP by timepoint
p_umap_time_pre <- DimPlot(combined, group.by = "timepoint")

ggsave(
  filename = file.path(plotdir, "Merged_Object_UMAP_preclustering_timepoint.png"),
  plot = p_umap_time_pre,
  width = 7,
  height = 6,
  dpi = 300
)

#UMAP by genotype
p_umap_genotype_pre <- DimPlot(combined, group.by = "genotype")

ggsave(
  filename = file.path(plotdir, "Merged_Object_UMAP_preclustering_genotype.png"),
  plot = p_umap_genotype_pre,
  width = 7,
  height = 6,
  dpi = 300
)

#UMAP without category
p_umap_blank <- DimPlot(combined)

ggsave(
  filename = file.path(plotdir, "Merged_Object_UMAP_preclustering_blank.png"),
  plot = p_umap_blank,
  width = 7,
  height = 6,
  dpi = 300
)

#Find clusters
combined <- FindClusters(combined, resolution = 0.1) #decide the best resolution at the end, here trial-error

p_umap_clusters <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE)

ggsave(
  filename = file.path(plotdir, "Merged_Object_UMAP_clusters_res0.1.png"),
  plot = p_umap_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

#Save the merged clean object
out_combined <- "/path/to/Merged_Object_combined_res0.1_pre_annotation.RDS"

saveRDS(combined, file = out_combined)

#--------------------------Finding markers to assign cell types to the clusters----------------------

#load the merged object
 combined <- readRDS("path/to/Merged_Object_combined_res0.1_pre_annotation.RDS")

#Find Cluster Markers (res=0.1)
combined <- JoinLayers(combined)

markers_res01 <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)

write.csv(
  markers_res01,
  file = file.path(plotdir, "Merged_Object_ClusterMarkers_res0.1.csv"),
  row.names = FALSE
)

#save the merged object with joined layers
saveRDS(
  combined,
  file = file.path(plotdir, "Merged_Object_combined_res0.1_joinedlayers.RDS")
)


#inspect markers for the main clusters
library(dplyr)

top10 <- markers_res01 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top10

top10 %>% filter(cluster == 0)
top10 %>% filter(cluster == 1)
top10 %>% filter(cluster == 2)
top10 %>% filter(cluster == 3)
top10 %>% filter(cluster == 4)
top10 %>% filter(cluster == 5)
top10 %>% filter(cluster == 6)
top10 %>% filter(cluster == 7)
top10 %>% filter(cluster == 8)
top10 %>% filter(cluster == 9)


#dot plot for the top 10 gene per cluster
p_dot <- DotPlot(combined, features = unique(top10$gene)) + RotatedAxis()

ggsave(
  filename = file.path(plotdir, "Merged_Object_DotPlot_top10_markers_res0.1.png"),
  plot = p_dot,
  width = 20,
  height = 8,
  dpi = 300
)

#DotPlots with the Canonical Markers
canonical_markers <- c(
  "Mki67", "Top2a","Sox2", "Pcna", # progenitors
  "Atoh1", "Otx2", "Pax6", "Pax3", #early (RL also)
  "Eomes", "Lmx1a", "Barhl1", "Barhl2", 
  "Trpc3", "Calb2", "Mgat5", "Hcrtr2", "Grm8", "Grm1", "En2",  #UBC related
  "En1", "Tlx3",     #GC related
  "Aldh1l1", "Gfap",  # glia
  "Slc17a7","Slc17a6","Slc17a8", # glutamatergic (excitatory)
          "Slc6a2",
          "Gad1","Gad2", # gabaergic (inhib)
          "Slc6a5","Slc32a1", #VGAT (vesicular GABA transporter) #inhibitory
          "Garb2","Gabrg3",
          "Th","Slc18a2","Dbh","Slc6a3","Nr4a2","Ppp1r1b", #dopa
          "Tph2","Slc6a4", #serotonergic
          "Chat","Slc5a7","Slc18a3", #motor neurons
          "Hoxb1","Hoxb2","Hoxb3", "Hoxb6", "Hoxb8") #spinalcord

p_dot <- DotPlot(combined, features = canonical_markers) + RotatedAxis()

#custom dotplot
p_dot_custom <- DotPlot(
  combined,
  features = canonical_markers,
  dot.scale = 8   # increases overall bubble size contrast
) +
  RotatedAxis() +
  scale_size_continuous(
    limits = c(0, 100),
    range = c(1, 12)   # makes small vs large bubbles very different
  )

ggsave(
  filename = file.path(plotdir, "Merged_Seurat_Object_DotPlot_canonical_markers_res0.1.png"),
  plot = p_dot_custom,
  width = 18,
  height = 8,
  dpi = 300
)


#Find more clusters (change resolution)
combined <- FindClusters(combined, resolution = 0.3)

p_umap_clusters <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE)

ggsave(
  filename = file.path(plotdir, "Merged_Object_UMAP_clusters_res0.3.png"),
  plot = p_umap_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

p_dot <- DotPlot(combined, features = canonical_markers) + RotatedAxis()

#custom dotplot
p_dot_custom <- DotPlot(
  combined,
  features = canonical_markers,
  dot.scale = 8   # increases overall bubble size contrast
) +
  RotatedAxis() +
  scale_size_continuous(
    limits = c(0, 100),
    range = c(1, 12)   # makes small vs large bubbles very different
  )

ggsave(
  filename = file.path(plotdir, "Merged_Seurat_Object_DotPlot_canonical_markers_res0.3.png"),
  plot = p_dot_custom,
  width = 18,
  height = 8,
  dpi = 300
)


#Visualize the distribution of each cluster across genotypes
tab_genotype <- table(combined$seurat_clusters, combined$genotype)

df_genotype <- as.data.frame(tab_genotype)
colnames(df_genotype) <- c("cluster", "genotype", "count")

library(ggplot2)

p_bar_counts <- ggplot(df_genotype, aes(x = factor(cluster), y = count, fill = genotype)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Number of cells",
    title = "Cell distribution by genotype per cluster (res = 0.3)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(plotdir, "Merged_Seurat_Object_StackedBar_cluster_by_genotype_counts_res0.3.png"),
  plot = p_bar_counts,
  width = 9,
  height = 6,
  dpi = 300
)

#Visualize the distribution of cells for each cluster across each sample
tab_sample <- table(combined$seurat_clusters, combined$sample)

df_sample <- as.data.frame(tab_sample)
colnames(df_sample) <- c("cluster", "sample", "count")

head(df_sample)

library(ggplot2)

p_bar_sample <- ggplot(df_sample, aes(x = factor(cluster), y = count, fill = sample)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    color = "white"
  ) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Number of cells",
    title = "Cell distribution by sample per cluster (res = 0.3)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(plotdir, "Merged_Seurat_Object_StackedBar_cluster_by_sample_counts_res0.3.png"),
  plot = p_bar_sample,
  width = 10,
  height = 6,
  dpi = 300
)

