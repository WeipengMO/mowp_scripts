# Whole-Genome Bisulfite Sequencing Data Processing Pipeline

Whole-genome bisulfite sequencing (WGBS) is used to discover methylation patterns at single-base resolution. Bisulfite treatment is used to convert unmethylated cytosines into uracils, but leaves methylated cytosines unchanged. After mapping bisulfite sequencing reads against a C-->U transformed genome, this pipeline can extract the CpG, CGH and CHH methylation patterns genome-wide.


## Inputs
| File format | Information contained in file | File description | 
|---------- |---------- |---------- |
| fastq | reads | Gzipped DNA-sequencing reads | 


## Output
| File format | Information contained in file | File description |
|---------- |---------- |---------- |
| bam | alignments | Produced by mapping reads to the genome |
| bigWig | Methylation sites coverage | Read coverage at CpG/CHG/CHH sites |
| *.methratio.txt.gz | methylation ratio | Prodcued by bsmapz `methratio.py` |

## Notes

For metaplot visualization, please refer to [metaplot_bs-seq.ipynb]('./notebooks/metaplot_bs-seq.ipynb')

## File details

```
.
├── aligned_data
│   └── col0_wgbs.bam
├── bw_files
│   ├── col0_wgbs.methratio.cg.bw
│   ├── col0_wgbs.methratio.chg.bw
│   └── col0_wgbs.methratio.chh.bw
├── config.yml
├── lambda_dna
│   └── Lambda.fa
├── methratio
│   ├── col0_wgbs_methratio.txt.gz
│   └── col0_wgbs_wiggle.txt.gz
├── notebooks
│   └── metaplot_bs-seq.ipynb
├── raw_data
│   ├── col0_wgbs_1.fastq.gz
│   └── col0_wgbs_2.fastq.gz
├── README.md
├── script
│   ├── filter_to_plot.py
│   ├── get_bedgraph.py
│   ├── metaplot_bs_seq.py
│   └── methylKit_DMR_v2.Rs
├── Snakefile
└── snakemake
    ├── bsmapz_pe.smk
    ├── bsmapz_se.smk
    ├── call_DMR.smk
    └── methylation_ratio.smk
```