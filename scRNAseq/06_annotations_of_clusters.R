
# =========================
# Libraries
# =========================
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
library(patchwork)
library(cowplot)


#object background

#Read the RDS
combined <- readRDS("/path/to/filtered/merged/seurat/object")

p_umap_clusters <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE)

#create a clean annotation column
combined$celltype <- "Unassigned"

#assign the clusters to cell type (based on the results of marker sets from published data+ FindAllmarkers)
combined$celltype[combined$seurat_clusters %in% c(18)] <- "GABAergic neurons"
combined$celltype[combined$seurat_clusters %in% c(17)] <- "Spinal cord neurons"
combined$celltype[combined$seurat_clusters %in% c(16)] <- "Otx2- Kcnip4- GC_2"
combined$celltype[combined$seurat_clusters %in% c(15)] <- "Tlx3+ GC_2"
combined$celltype[combined$seurat_clusters %in% c(14)] <- "RL_VZ"
combined$celltype[combined$seurat_clusters %in% c(13)] <- "isth_N"
combined$celltype[combined$seurat_clusters %in% c(12)] <- "UBC_Hcrtr2"
combined$celltype[combined$seurat_clusters %in% c(11)] <- "Otx2- Kcnip4- GC_1"
combined$celltype[combined$seurat_clusters %in% c(8)] <- "GC/UBC diff_2"
combined$celltype[combined$seurat_clusters %in% c(9)] <- "GC/UBC diff_1"
combined$celltype[combined$seurat_clusters %in% c(3)] <- "GC/UBC pro_1"
combined$celltype[combined$seurat_clusters %in% c(10)] <- "GC/UBC pro_2"
combined$celltype[combined$seurat_clusters %in% c(7)] <- "UBC_defined"
combined$celltype[combined$seurat_clusters %in% c(6)] <- "UBC_cKO"
combined$celltype[combined$seurat_clusters %in% c(5)] <- "Tlx3+ GC_diff"
combined$celltype[combined$seurat_clusters %in% c(4)] <- "Tlx3+ GC_1"
combined$celltype[combined$seurat_clusters %in% c(2)] <- "Otx2+ Kcnip4+ GC_1"
combined$celltype[combined$seurat_clusters %in% c(1)] <- "UBC_diff"
combined$celltype[combined$seurat_clusters %in% c(0)] <- "Otx2+ Kcnip4+ GC_2"


#check if the clusters are assigned
table(combined$celltype)



#for figure 3 main 
markers <- c(
"Mki67","Top2a","Sox2","Pcna","Atoh1","Otx2","Pax6",
"Eomes","Lmx1a","Barhl1","Barhl2",
"Calb2","Mgat5","Trpc3","Grm1", "Grm8","En1","En2",
"Hcrtr2","Cxcr4",
"Synpr","Neurod1","Zic1","Zic2","Kcnip4","Dcx",
"Tlx3", "Pde1c", "Cntn1", "Pcsk9", "Adarb2", "Slc17a6", "Slc17a7")

celltype_order <- c(
"GC/UBC pro_1",
"GC/UBC pro_2",
"UBC_diff",
"UBC_Hcrtr2",
"UBC_defined",
"UBC_cKO",
"GC/UBC diff_1",
"GC/UBC diff_2",
"Tlx3+ GC_diff",
"Tlx3+ GC_1",
"Tlx3+ GC_2",
"Otx2+ Kcnip4+ GC_1",
"Otx2+ Kcnip4+ GC_2",
"Otx2- Kcnip4- GC_1",
"Otx2- Kcnip4- GC_2"
)


Idents(combined) <- "celltype"

combined <- subset(combined, idents = celltype_order)

combined$celltype <- factor(
  combined$celltype,
  levels = celltype_order
)

Idents(combined) <- "celltype"


library(Seurat)
library(ggplot2)
library(viridis)

magma_rev <- rev(viridis::magma(256))

p <- DotPlot(
  combined,
  features = markers
) +
  RotatedAxis() +
  theme_bw() +
  labs(x = "Markers", y = "Cell types")

p +
  scale_color_gradientn(colors = magma_rev) +
  scale_size(
    range = c(0, 12),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  )

final_plot <- p +
  scale_color_gradientn(colors = magma_rev) +
  scale_size(
    range = c(0, 12),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


#no grid
library(Seurat)
library(ggplot2)
library(viridis)

magma_rev <- rev(viridis::magma(256))

final_plot <- DotPlot(
  combined,
  features = markers
) +
  RotatedAxis() +
  theme_bw() +
  labs(x = "Markers", y = "Cell types") +
  scale_color_gradientn(colors = magma_rev) +
  scale_size(
    range = c(0, 12),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "/path/to/Merged_Object_Figure3_Main_dotplot_v4.pdf",
  plot = final_plot,
  width = 14,
  height = 6,
  units = "in",
  device = cairo_pdf  
)



#for supp
markers <- c(
"Mki67","Top2a","Sox2","Pcna","Atoh1","Otx2","Pax6", "Pax3","Aldh1l1", "Gfap", "Lmx1a",
"Slc17a17", "Slc17a8", "Slc5a7", "Hoxb2", "Hoxb3",  "Slc6a5", 
 "Gad1", "Gad2", "Gabrg3", "Scg2", "Slc17a6", "Pax5")

celltype_order <- c(
"RL_VZ",
"Spinal cord neurons",
"GABAergic neurons",
"isth_N"
)


#wo grid for supp

Idents(combined) <- "celltype"

combined <- subset(combined, idents = celltype_order)

combined$celltype <- factor(
  combined$celltype,
  levels = celltype_order
)

Idents(combined) <- "celltype"

magma_rev <- rev(viridis::magma(256))

final_plot <- DotPlot(
  combined,
  features = markers
) +
  RotatedAxis() +
  theme_bw() +
  labs(x = "Markers", y = "Cell types") +
  scale_color_gradientn(colors = magma_rev) +
  scale_size(
    range = c(0, 12),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(
  "/path/to/Merged_Object_Figure3_Supp_dotplot_v4.pdf",
  plot = final_plot,
  width = 10,
  height = 4,
  units = "in",
  device = cairo_pdf)

