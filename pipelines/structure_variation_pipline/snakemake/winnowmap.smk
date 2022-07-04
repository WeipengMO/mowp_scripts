# Make reference index  
rule ref_index:
    output:
        "genome.k15.txt"
    params:
        meryl=config['meryl'],
        genome=config['genome']
    threads: 32
    shell:
        '''
{params.meryl} count threads={threads} k=15 output genome.merylDB {params.genome}
{params.meryl} print greater-than distinct=0.9998 genome.merylDB > {output}
        '''


##Align with winnowmap
rule align:
    input:
        fastq='raw_data/{sample}.fastq.gz',
        refidx=rules.ref_index.output
    output:
        bam ="aligned_data/{sample}.sorted.bam",
        bai ="aligned_data/{sample}.sorted.bam.bai"
    threads: 58
    params:
        genome=config['genome']
    conda:
        'genome'
    shell:
        '''
winnowmap -t {threads} -W {input.refidx} --MD -ax map-ont {params.genome} {input.fastq} | samtools view -b -u -F 256 | samtools sort -o {output.bam}
samtools index -@ {threads} {output.bam}
        '''