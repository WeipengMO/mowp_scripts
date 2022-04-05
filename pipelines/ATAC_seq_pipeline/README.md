# ATAC-seq Pipelines

ATAC-seq (Assay for Transposase Accessible Chromatin with high-throughput sequencing) is viewed as an alternative to DNase-seq and MNase-seq; it probes DNA accessibility with hyperactive Tn5 transposase, which inserts sequencing adapters into accessible regions of chromatin.

## Inputs
| File format | Information contained in file | File description | 
|---------- |---------- |---------- |
| fastq | reads | G-zipped reads, paired-ended or single ended, stranded or unstranded. | 
| fasta | genome indices (bowtie2) | Indices are dependent on the assembly being used for mapping | 

## Output
| File format | Information contained in file | File description |
|---------- |---------- |---------- |
| bigWig | raw signal (bw_files), fold change over control (bw_compare) | signal coverage tracks |
| narrowPeak | peaks | Peak calls for each replicate individually. |
| png | figure | metaplot around the gene body |

## Notes
Due to issue with MASC2 installation, it is recommended to use singularity to run conda:

You can also pass additional arguments to singularity, including bind points, like this:

`snakemake --use-singularity --singularity-args "-B /path/outside/container/:/path/inside/container/"`

More details, you can see the [Snakemake + docker example, how to use volumes](https://stackoverflow.com/questions/52742698/snakemake-docker-example-how-to-use-volumes)