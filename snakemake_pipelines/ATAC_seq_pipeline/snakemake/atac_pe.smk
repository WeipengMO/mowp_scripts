configfile: 
    "config.yml"

suffix = config['suffix']

rule run_fastp:
    input:
        fq1='raw_data/{sample_name}_1'+f'{suffix}',
        fq2='raw_data/{sample_name}_2'+f'{suffix}',
    output:
        fq1=temp('raw_data/{sample_name}_R1.clean.fq.gz'),
        fq2=temp('raw_data/{sample_name}_R2.clean.fq.gz')
    params:
        html='raw_data/{sample_name}.html',
        json='raw_data/{sample_name}.json'
    threads: 16
    shell:
        '''
fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} -h {params.html} -j {params.json}
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
    threads: 30
    conda: 'chipseq'
    shell:
        '''
bowtie2 -t -p {threads} --very-sensitive -X 2000 -x {params.genome} -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -O bam -o {output.bam} -
samtools index -@ {threads} {output.bam}
        '''


rule MarkDuplicates:
    input:
        'aligned_data/{sample_name}.sorted.bam'
    output:
        bam='aligned_data/{sample_name}.sorted.rmdup.bam',
        bai='aligned_data/{sample_name}.sorted.rmdup.bam.bai'
    params:
        picard=config['picard_path']
    threads: 8
    shell:
        '''
java -jar {params.picard} MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I={input} O={output.bam} M={output.bam}.markdump.txt
samtools index -@ 10 {output.bam}
        '''


rule bamCoverage:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        'bw_files/{sample_name}.sorted.rmdup.CPM.bw'
    threads: 16
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 10 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors {threads}
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
    shell:
        '''
computeMatrix scale-regions -b 1000 -a 1000 -R {params} -S {input} --skipZeros -o {output.matrix} -p {threads}
plotProfile -m {output.matrix} -out {output.png}
        '''


rule macs2_callpeak:
    input:
        treatment='aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        treat='macs2_result/{sample_name}_peaks.narrowPeak',
    threads: 1
    params:
        name='{sample_name}',
        out_dir='macs2_result/',
        gsize=config['gsize']
    conda:
        'chipseq'
    shell:
        '''
macs2 callpeak -t {input.treatment} -f BAM -g {params.gsize} -n {params.name} --nomodel --shift -100 --extsize 200 --outdir {params.out_dir}
        '''


rule run_Genrich:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        bam=temp('aligned_data/{sample_name}.sorted.name.bam'),
        peak='genrich_result/{sample_name}_peaks.narrowPeak',
    params:
        exclude_chroms=config['exclude_chroms']
    threads: 1
    shell:
        '''
samtools sort -n -@ 30 -O bam -o {output.bam} {input}
Genrich -t {output.bam} -o {output.peak} -j  -y  -r  -e {params.exclude_chroms}  -v
        '''