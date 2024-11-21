configfile: 
    "config.yml"

suffix = config['suffix']
sep1 = config['sep1']
sep2 = config['sep2']

print(suffix, sep1, sep2)

rule run_fastp:
    input:
        fq1='raw_data/{sample_name}'+sep1+suffix,
        fq2='raw_data/{sample_name}'+sep2+suffix,
    output:
        fq1=temp('raw_data/{sample_name}.1.clean.fastq.gz'),
        fq2=temp('raw_data/{sample_name}.2.clean.fastq.gz')
    params:
        html='raw_data/{sample_name}.html',
        json='raw_data/{sample_name}.json'
    log:
        'logs/fastp/{sample_name}.log'
    threads: 8
    shell:
        '''
fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} -h {params.html} -j {params.json} &> {log}
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
    threads: 30
    log:
        'logs/hisat2/{sample_name}.log'
    conda: 'rnaseq'
    shell:
        '''
hisat2 -x {params.genome} -p {threads} --dta --time -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -O bam -o {output.bam} - && samtools index {output.bam}
        '''


rule MarkDuplicates:
    input:
        'aligned_data/{sample_name}.sorted.bam'
    output:
        bam='aligned_data/{sample_name}.sorted.rmdup.bam',
        bai='aligned_data/{sample_name}.sorted.rmdup.bam.bai'
    threads: 8
    params:
        picard=config['picard_path']
    conda: 'rnaseq'
    log:
        'logs/picard/{sample_name}.log'
    shell:
        '''
java -jar {params.picard} MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I={input} O={output.bam} M={output.bam}.markdump.txt &> {log} && samtools index {output.bam}
        '''
