import loompy as lp
import numpy as np
import scanpy as sc
import argparse


def prepare_for_scenic(adata, loom_out):
    
    row_attributes = {"Gene": np.array(adata.var.index)}
    col_attributes = {
        "CellID": np.array(adata.obs.index),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),}
    
    lp.create(loom_out, adata.X.transpose(), row_attributes, col_attributes)


def main():

    parser = argparse.ArgumentParser(description='Prepare loom file for pySCENIC')
    parser.add_argument('adata_file', help='h5ad file')
    parser.add_argument('loom_out', help='loom file')
    parser.add_argument('--layer', default='counts', help='Raw count layer to use')
    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata_file)
    # use the counts matrix (without log transformation or further processing)
    # https://github.com/aertslab/pySCENIC/issues/128
    if args.layer not in adata.layers.keys():
        raise ValueError(f"Layer {args.layer} not found in adata object")
    
    if args.layer != 'X':
        adata.X = adata.layers[args.layer]

    prepare_for_scenic(adata, args.loom_out)


if __name__ == "__main__":
    main()