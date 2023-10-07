reference = config['genome']

rule lra_align:
    input:
        'raw_data/{sample}.fastq.gz'
    output:
        bam='aligned_data/{sample}.sorted.bam',
        bai='aligned_data/{sample}.sorted.bam.bai',
    threads: 48
    conda:
        'genome'
    shell:
        '''
seqtk seq -A {input} | lra align -ONT -p s -t 48 {reference} - | samtools calmd -u - {reference} | samtools sort -@ {threads} -o {output.bam} -
samtools index -@ {threads} {output.bam}
        '''