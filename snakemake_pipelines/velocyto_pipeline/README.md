# single-cell RNA sequencing analysis

## Description

This pipeline is used to process single-cell RNA sequencing data. It is based on the [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10x Genomics. 

## Trajectory Analysis using 10x Genomics Single Cell Gene Expression Data

`velocyto.smk` use RNA velocity methods to analyze single cell gene expression data and infer the dynamic processes of cells. The pipeline is based on [Trajectory Analysis](https://www.10xgenomics.com/resources/analysis-guides/trajectory-analysis-using-10x-Genomics-single-cell-gene-expression-data).

The input folder should contain the output of the Cell Ranger pipeline:

```
.
├── sample1
│   └── sample1
│       ├── outs
│       │   ├── analysis
│       │   ├── filtered_feature_bc_matrix
│       │   └── raw_feature_bc_matrix
│       └── SC_RNA_COUNTER_CS
│           ├── CELLRANGER_PREFLIGHT
│           ├── CELLRANGER_PREFLIGHT_LOCAL
│           ├── fork0
│           ├── FULL_COUNT_INPUTS
│           ├── GET_AGGREGATE_BARCODES_OUT
│           ├── SC_MULTI_CORE
│           ├── _STRUCTIFY
│           └── WRITE_GENE_INDEX
└── sample2
    └── sample2
        ├── outs
        │   ├── analysis
        │   ├── filtered_feature_bc_matrix
        │   └── raw_feature_bc_matrix
        └── SC_RNA_COUNTER_CS
            ├── CELLRANGER_PREFLIGHT
            ├── CELLRANGER_PREFLIGHT_LOCAL
            ├── fork0
            ├── FULL_COUNT_INPUTS
            ├── GET_AGGREGATE_BARCODES_OUT
            ├── SC_MULTI_CORE
            ├── _STRUCTIFY
            └── WRITE_GENE_INDEX
```