# HiC Data Processing Pipeline

3D-chromatin pipelines such as HiC and ChIA-PET are assays that mapping the organization of genome and to examine the spatial proximity of chromosomal loci. HiC assay provides the three-dimensional architecture of whole genomes by coupling proximity-based ligation with massively parallel sequencing. 


## Inputs
| File format | Information contained in file | File description | 
|---------- |---------- |---------- |
| fastq | reads | G-zipped reads, paired-ended or single ended, stranded or unstranded. | 
| fasta | genome indices (bowtie2) | Indices are dependent on the assembly being used for mapping | 

## Output
| File format | Information contained in file | File description |
|---------- |---------- |---------- |
| cool | genomic matrix data | HDF5 as the container format |
| mcool | multi-resolution | There is a special `@` syntax for using these files, which allow you to set the resolution which should be used for this analysis |

## File details

```
.
├── config.yaml
├── hindIII_hic_SRR1504819.hic.bam
├── hindIII_hic_SRR1504819.hic.q1.resolution_1000.cool
├── hindIII_hic_SRR1504819.hic.q1.resolution_1000.mcool
├── hindIII_hic_SRR1504819.hic.q1.sorted.rmdup.pairsam.gz
├── hindIII_hic_SRR1504819.hic.q1.sorted.rmdup.pairsam.gz.px2
├── notebooks
│   └── cooltools_visualization.ipynb
├── raw_data  # input folder
│   ├── hindIII_hic_SRR1504819_1.fastq.gz
│   └── hindIII_hic_SRR1504819_2.fastq.gz
├── README.md
└── Snakefile
```