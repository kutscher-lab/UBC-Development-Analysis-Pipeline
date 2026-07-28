
# =========================================================

# run_pairwise_DEGs.R
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

library(Seurat)

# INPUT / OUTPUT

input_rds <- "/path/to/seurat_object.rds"

outdir <- "/path/t/DEGs_pairwise"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outdir_wilcox <- file.path(outdir, "Wilcoxon")

dir.create(outdir_wilcox, recursive = TRUE, showWarnings = FALSE)

# LOAD OBJECT

message("Loading Seurat object...")
combined <- readRDS(input_rds)


# DEFINE COMPARISONS

comparisons <- list(

    WT_1_7_vs_KO_6 = list(
        wt = c(1, 7),
        ko = c(6)
    )
)

# DE FUNCTION

run_custom_DE <- function(
    seu,
    wt_clusters,
    ko_clusters,
    tp,
    comparison_name
){

    wt_cells <- WhichCells(
        seu,
        expression =
            genotype %in% c("Tbr2_WT", "ires_gfp") &
            timepoint == tp &
            seurat_clusters %in% wt_clusters
    )

    ko_cells <- WhichCells(
        seu,
        expression =
            genotype == "Tbr2_flx" &
            timepoint == tp &
            seurat_clusters %in% ko_clusters
    )

    message(
        paste(
            comparison_name,
            "|", tp,
            "| WT cells =", length(wt_cells),
            "| KO cells =", length(ko_cells)
        )
    )

    if(length(wt_cells) < 10 | length(ko_cells) < 10){
        message("Too few cells. Skipping.")
        return(NULL)
    }

print(length(wt_cells))
print(length(ko_cells))

print(head(wt_cells))
print(head(ko_cells))
    

   # Create temporary object containing only
   # cells used in this comparison

    tmp <- subset(
        seu,
        cells = c(wt_cells, ko_cells)
    )

    tmp$DE_group <- "WT"
    tmp$DE_group[colnames(tmp) %in% ko_cells] <- "KO"

    Idents(tmp) <- "DE_group"

    print(table(Idents(tmp)))


    # Differential expression
    

    deg <- FindMarkers(
        object = tmp,
        ident.1 = "KO",
        ident.2 = "WT",
        assay = "RNA",
        slot = "data",
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0.05,
        return.thresh = 1
    )

    deg$gene <- rownames(deg)

    deg$comparison <- comparison_name
    deg$timepoint <- tp
    deg$wt_cells <- length(wt_cells)
    deg$ko_cells <- length(ko_cells)

    return(deg)
}
    

# RUN ALL COMPARISONS

timepoints <- c("E16", "P0")

for(tp in timepoints){

    for(comp_name in names(comparisons)){

        comp <- comparisons[[comp_name]]

        res <- tryCatch(

            run_custom_DE(
                seu = combined,
                wt_clusters = comp$wt,
                ko_clusters = comp$ko,
                tp = tp,
                comparison_name = comp_name
            ),

            error = function(e){
                message("ERROR: ", e$message)
                NULL
            }
        )

        if(!is.null(res)){

            outfile <- file.path(
                outdir_wilcox,
                paste0(
                    comp_name,
                    "_",
                    tp,
                    "_Wilcox.csv"
                )
            )

            write.csv(
                res,
                outfile,
                row.names = FALSE
            )

            message("Saved: ", outfile)
        }
    }
}

message("All DE analyses completed successfully.")

