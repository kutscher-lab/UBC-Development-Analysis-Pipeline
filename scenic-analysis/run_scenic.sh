#!/bin/bash

# =========================================================

# run_scenic.sh
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

#In the next steps, we will need the mouse reference SCENIC/RcisTarget resources that were downloaded while installing SCENIC. Make sure they are in your directories.

# 1. Copy Reference Files to your SCENIC directories if necessary

cp mm_mgi_tfs.txt 
motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl 
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
/scenic_dir/WT_SCENIC/

cp mm_mgi_tfs.txt 
motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl 
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
/scenic_dir/FLOX_SCENIC/


#In the next step, we will create a bash script to submit it to the cluster as a job. This script will run the 1st major step of SCENIC: Gene Regulatory Network (GRN) inference. 
#GRN will try to find which TFs appear to regulate which genes based on the gene expression patterns across all cells. In this step, pySCENIC does not yet use the motif information or genome annotations. This is only predicted from the expression data.


# 2. WT GRN Script

cd /scenic_dir/WT_SCENIC

echo "Starting WT GRN"
date

/path/to/scenic_env/bin/pyscenic grn 
wt_scenic_input.loom 
mm_mgi_tfs.txt 
-o adjacencies.tsv 
--num_workers 8

echo "WT GRN complete"
date


# 3. FLOX GRN Script

cd /scenic_dir/FLOX_SCENIC

echo "Starting FLOX GRN"
date

/path/to/scenic_env/bin/pyscenic grn 
flox_scenic_input.loom 
mm_mgi_tfs.txt 
-o adjacencies.tsv 
--num_workers 8

echo "FLOX GRN complete"
date

# 4. Verify GRN Output


# After GRN jobs completed successfully:

ls -lh adjacencies.tsv

# Expected output:

# adjacencies.tsv present and non-empty




# 5. Major Step 2: Finding regulons

#WT
cd /scenic_dir/WT_SCENIC


echo "Starting WT ctx"
date

/path/to/pyscenic ctx 
adjacencies.tsv 
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
--annotations_fname motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl 
--expression_mtx_fname wt_scenic_input.loom 
--mode dask_multiprocessing 
--output regulons.csv 
--num_workers 8

echo "WT ctx complete"
date


# FLOX 


cd /scenic_dir/FLOX_SCENIC


echo "Starting FLOX ctx"
date

/path/to/pyscenic ctx 
adjacencies.tsv 
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather 
--annotations_fname motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl 
--expression_mtx_fname flox_scenic_input.loom 
--mode dask_multiprocessing 
--output regulons.csv 
--num_workers 8

echo "FLOX ctx complete"
date



#Expected ctx Output

# Main output: regulons.csv

# This contains:

# - motif-pruned regulons

# - high-confidence TF-target networks


# 6. Major Step 3: Calculating regulon activity scores


#WT

cd /path/to/scenic_dir/WT_SCENIC

echo "Starting WT AUCell"
date

/path/to/scenic_env/pyscenic aucell \
    wt_scenic_input.loom \
    regulons.csv \
    -o auc_mtx.csv \
    --num_workers 8

echo "WT AUCell complete"
date

cd /path/to/scenic_dir/FLOX_SCENIC

echo "Starting FLOX AUCell"
date

/path/to/scenic_env/pyscenic aucell \
    flox_scenic_input.loom \
    regulons.csv \
    -o auc_mtx.csv \
    --num_workers 8

echo "FLOX AUCell complete"
date

# Output: auc_mtx.csv

# This will provide: regulon activity scores per cell.
