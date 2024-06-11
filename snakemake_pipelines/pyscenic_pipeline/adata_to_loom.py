import loompy as lp
import numpy as np
import scanpy as sc
import sys


def prepare_for_scenic(adata, loom_out):
    
    row_attributes = {"Gene": np.array(adata.var.index)}
    col_attributes = {
        "CellID": np.array(adata.obs.index),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),}
    
    lp.create(loom_out, adata.X.transpose(), row_attributes, col_attributes)


def main():
    adata_file = sys.argv[1]
    loom_out = sys.argv[2]
    
    adata = sc.read_h5ad(adata_file)
    prepare_for_scenic(adata, loom_out)


if __name__ == "__main__":
    main()