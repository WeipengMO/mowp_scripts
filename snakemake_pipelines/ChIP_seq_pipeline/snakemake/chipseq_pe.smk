configfile: 
    "config.yml"

suffix = config['suffix']
sep = config['sep']

rule run_fastp:
    input:
        fq1='raw_data/{sample_name}'+f'{sep[:-1]}1{suffix}',
        fq2='raw_data/{sample_name}'+f'{sep[:-1]}2{suffix}',
    output:
        fq1=temp('raw_data/{sample_name}_R1.clean.fq.gz'),
        fq2=temp('raw_data/{sample_name}_R2.clean.fq.gz')
    params:
        html='raw_data/{sample_name}.html',
        json='raw_data/{sample_name}.json'
    log:
        'logs/{sample_name}.fastp.log'
    threads: 16
    shell:
        '''
fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} -h {params.html} -j {params.json} &> {log}
        '''
        
# paired-end
rule run_bowtie2_pe:
    input:
        fq1='raw_data/{sample_name}_R1.clean.fq.gz',
        fq2='raw_data/{sample_name}_R2.clean.fq.gz',
    output:
        bam=temp('aligned_data/{sample_name}.sorted.bam'),
        bai=temp('aligned_data/{sample_name}.sorted.bam.bai')
    params:
        genome=config['genome']
    log:
        'logs/{sample_name}.bowtie2.log'
    threads: 30
    conda: 'chipseq'
    shell:
        '''
bowtie2 -t -p {threads} --dovetail -X 1000 -x {params.genome} -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -O bam -o {output.bam} - &> {log}
samtools index -@ {threads} {output.bam}
        '''


rule MarkDuplicates:
    input:
        'aligned_data/{sample_name}.sorted.bam'
    output:
        bam='aligned_data/{sample_name}.sorted.rmdup.bam',
        bai='aligned_data/{sample_name}.sorted.rmdup.bam.bai'
    threads: 8
    log:
        'logs/{sample_name}.MarkDuplicates.log'
    params:
        picard=config['picard_path']
    shell:
        '''
java -jar {params.picard} MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I={input} O={output.bam} M={output.bam}.markdump.txt &> {log}
samtools index -@ 10 {output.bam}
        '''


rule bamCoverage:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        'bw_files/{sample_name}.sorted.rmdup.CPM.bw'
    threads: 16
    log:
        'logs/{sample_name}.bamCoverage.log'
    params:
        gsize=config['gsize']
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 10 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors {threads} &> {log}
        '''


rule computeMatrix:
    input:
        'bw_files/{sample_name}.sorted.rmdup.CPM.bw',
    output:
        matrix=temp('deeptools_profile/{sample_name}.matrix.gz'),
        png='deeptools_profile/{sample_name}.scale.png'
    params:
        config['bed']
    threads: 16
    log:
        'logs/{sample_name}.computeMatrix.log'
    shell:
        '''
computeMatrix scale-regions -b 1000 -a 1000 -R {params} -S {input} --skipZeros -o {output.matrix} -p {threads} &> {log}
plotProfile -m {output.matrix} -out {output.png} &> {log}
        '''


