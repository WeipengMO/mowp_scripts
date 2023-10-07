'''
Author       : windz
Date         : 2020-05-18 20:38:48
LastEditTime : 2021-12-08 14:20:05
'''


rule run_fastp:
    input:
        fq1='raw_data/{sample_name}_1.fastq.gz',
        fq2='raw_data/{sample_name}_2.fastq.gz',
    output:
        fq1=temp('raw_data/{sample_name}.1.clean.fastq.gz'),
        fq2=temp('raw_data/{sample_name}.2.clean.fastq.gz')
    params:
        html='raw_data/{sample_name}.html',
        json='raw_data/{sample_name}.json'
    threads: 16
    shell:
        '''
fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} -h {params.html} -j {params.json}
        '''


rule run_hisat2_pe:
    input:
        fq1='raw_data/{sample_name}.1.clean.fastq.gz',
        fq2='raw_data/{sample_name}.2.clean.fastq.gz',
    output:
        bam=temp('aligned_data/{sample_name}.sorted.bam'),
        bai=temp('aligned_data/{sample_name}.sorted.bam.bai')
    params:
        genome=config['genome'],
        max_intronlen=config['max_intronlen'],
    threads: 30
    shell:
        '''
hisat2 -x {params.genome} -p {threads} --min-intronlen 20 --max-intronlen {params.max_intronlen} --dta --time -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -O bam -o {output.bam} - && samtools index {output.bam}
        '''


rule MarkDuplicates:
    input:
        'aligned_data/{sample_name}.sorted.bam'
    output:
        bam='aligned_data/{sample_name}.sorted.rmdup.bam',
        bai='aligned_data/{sample_name}.sorted.rmdup.bam.bai'
    threads: 8
    shell:
        '''
java -jar /public/apps/picard_2.20.2/picard.jar MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I={input} O={output.bam} M={output.bam}.markdump.txt && samtools index {output.bam}
        '''
