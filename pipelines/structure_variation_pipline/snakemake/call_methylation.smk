'''
Author       : windz
Date         : 2022-04-02 14:34:44
LastEditTime : 2022-04-02 15:24:34
LastEditors  : windz
FilePath     : call_methylation.smk
'''

rule index:
    input: 
        fq='{sample}.fastq.gz',
        f5=config['fast5dir']
    output:
        '{sample}.fastq.gz.index.readdb'
    conda:
        "methy"
    shell:
        '''
nanopolish index -d {input.f5} {input.fq}
        '''


rule call_methylation:
    input:
        fq='raw_data/{sample}.fastq.gz',
        index='raw_data/{sample}.fastq.gz.index.readdb',
        bam='aligned_data/{sample}.sorted.bam'
    output:
        "mcall/{sample}.{mod}.tsv"
    params:
        genome=config['genome'],
    threads:
        48
    conda:
        "methy"
    shell:
        '''
nanopolish call-methylation -q {wildcards.mod} -t {threads} -r {input.fq} -g {params.genome} -b {input.bam} > {output}
        '''