# Structural variation Pipeline

This workflow provides an easy way to call structural variants in human genomic data.

The pipeline performs the following steps:

- Maps reads using lra
- Calls variants using cuteSV
- Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)

## Inputs

| File format | Information contained in file | File description | 
|---------- |---------- |---------- |
| fastq | reads | basecalled nanopore reads | 

## Output
| File format | Information contained in file | File description |
|---------- |---------- |---------- |
| vcf | Variant Call Format file | [Introduction_to_Variant_Call_Format_(vcf)_files](https://github.com/epi2me-labs/tutorials/blob/master/tutorials/Introduction_to_Variant_Call_Format_(vcf)_files.ipynb) |
| bigWig | read coverage | signal coverage tracks |

## File details

```
.
├── aligned_data
│   ├── GM24385.nf7.chr20_af_minimap2.sorted.bam
│   └── GM24385.nf7.chr20_af_minimap2.sorted.bam.bai
├── config.yml
├── cuteSV
│   ├── GM24385.nf7.chr20_af_minimap2.cuteSV.vcf
│   ├── GM24385.nf7.chr20_af_minimap2.mosdepth.global.dist.txt
│   └── GM24385.nf7.chr20_af_minimap2.mosdepth.summary.txt
├── raw_data
│   └── GM24385.nf7.chr20_af_minimap2.fastq.gz
├── README.md
├── scripts
│   └── get_filter_calls_command.py
├── Snakefile
└── snakemake
    ├── call_methylation.smk
    ├── call_sv.smk
    ├── lra.smk
    ├── nanonome_data_parse.smk
    └── winnowmap.smk

```

## Notes

This pipeline is adapted from [wf-human-sv](https://github.com/epi2me-labs/wf-human-sv).
