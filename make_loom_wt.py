
# =========================================================

# make_loom_wt.py
# SCENIC Workflow for Separate WT and cKO Data
# 16 July 2026
# Author: Vuslat Akçay

# =========================================================

#So far in R, we preprocessed our data and created the inputs that SCENIC needs for the 1st step.
#SCENIC will start from gene x cell expression matrix to infer: 1. gene co-expression modules, 2. identify transcription factor regulons, and 3. score regulon activity in each cell.

#Make loom files for SCENIC for each genotype and run GRN inference with pySCENIC

#In the next step, I use a virtual environment (micromamba) for SCENIC on the cluster to run the python commands.

#1. Activate pySCENIC Environment

# micromamba activate /path/to/scenic_env


#2. Create WT Loom File

#cd /scenic_dir/WT_SCENIC


import loompy
import scipy.io
import numpy as np

mat = scipy.io.mmread("matrix.mtx").tocsc()

genes = np.array([x.strip() for x in open("genes.tsv")])
cells = np.array([x.strip() for x in open("barcodes.tsv")])

loompy.create(
"wt_scenic_input.loom",
mat,
row_attrs={"Gene": genes},
col_attrs={"CellID": cells}
)

print("WT loom file created successfully")
