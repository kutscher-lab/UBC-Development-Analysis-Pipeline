
# =========================================================

# preprocess_seurat_data.R
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

#The goal of this analysis to compare the regulatory network changes between a wild-type and a conditional knockout (cKO) of a TF of interest (in our case Eomes). The seurat object used for this analysis is the merged dataset from 2 genotypes and 2 timepoints for each genotype. First step is to separate the genotypes for WT and cKO for the downstream analyses.

#Preprocessing the Seurat object and generating a count matrix for SCENIC

# 1. Activate R / Seurat Environment in terminal on HPC
# module load R/4.4.3-GCCcore-14.1.0
# R

# 2. Load Libraries in R

library(Seurat)
library(SeuratObject)
library(Matrix)

# 3. Define Paths

outdir <- "/path/to/your/dir"

scenic_dir <- file.path(outdir, "SCENIC_separate")

dir.create(scenic_dir, showWarnings = FALSE)


# 4. Read Seurat Object

obj <- readRDS(
file.path(outdir, "seurat_obj.rds")
)


# 5. Check TF of interest Expression per genotype (optional)

AverageExpression(
obj,
features = "Eomes",
group.by = "genotype"
)

AverageExpression(
obj,
features = "Eomes",
group.by = "RNA_snn_res.0.7"
)


#In this next step, I am subsetting the MERGED data into separate objects based on their genotypes.

# 6. Create WT and FLOX Subsets

obj_wt <- subset(
obj,
subset = genotype == "Tbr2_WT"
)

obj_flox <- subset(
obj,
subset = genotype == "Tbr2_flx"
)

 
# 7. Save Subset Objects

saveRDS(
obj_wt,
file.path(scenic_dir, "obj_wt_res0.7.rds")
)

saveRDS(
obj_flox,
file.path(scenic_dir, "obj_flox_res0.7.rds")
)


#The preprocessing of the seurat object is done now. In the next part, the first step of the SCENIC starts, which requires the count matrices from the whole dataset for each genotype. At the end of this step, for each genotype, we should have: 1.matrix.mtx, 2. genes.tsv, and 3.barcodes.tsv files.

#1. matrix.csv: This contains only the numbers, the count matrix. scRNAseq matrices are sparse (containing zeros), so this mtx format stores the non-zero, expression values.

#2. genes.tsv: stores the gene names exported from the count matrix.

#3. barcodes.tsv: stores the cell IDs/barcodes exported from the count matrix.


# 8. Extract Counts Matrices

exprMat_wt <- GetAssayData(
obj_wt,
assay = "RNA",
layer = "counts"
)

exprMat_flox <- GetAssayData(
obj_flox,
assay = "RNA",
layer = "counts"
)


# 9. Create SCENIC Directories

wt_dir <- file.path(scenic_dir, "WT_SCENIC")

flox_dir <- file.path(scenic_dir, "FLOX_SCENIC")

dir.create(wt_dir, showWarnings = FALSE)

dir.create(flox_dir, showWarnings = FALSE)


# 10. Export WT Matrix Files

writeMM(
exprMat_wt,
file.path(wt_dir, "matrix.mtx")
)

write.table(
rownames(exprMat_wt),
file.path(wt_dir, "genes.tsv"),
quote = FALSE,
row.names = FALSE,
col.names = FALSE
)

write.table(
colnames(exprMat_wt),
file.path(wt_dir, "barcodes.tsv"),
quote = FALSE,
row.names = FALSE,
col.names = FALSE
)


# 10. Export FLOX Matrix Files


writeMM(
exprMat_flox,
file.path(flox_dir, "matrix.mtx")
)

write.table(
rownames(exprMat_flox),
file.path(flox_dir, "genes.tsv"),
quote = FALSE,
row.names = FALSE,
col.names = FALSE
)

write.table(
colnames(exprMat_flox),
file.path(flox_dir, "barcodes.tsv"),
quote = FALSE,
row.names = FALSE,
col.names = FALSE
)


#Analysis is continued in make_loom.py script

sessionInfo()
