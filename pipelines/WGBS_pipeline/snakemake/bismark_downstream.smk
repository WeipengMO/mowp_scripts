rule deduplicate_bismark:
    input:
        'bismark_aligned_data/{sample_name}_bismark_bt2.bam',
    output:
        temp('bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.bam')
    threads: 1
    shell:
        '''
export PATH=/public/home/mowp/softwares/Bismark:$PATH

deduplicate_bismark --bam {input} --output_dir bismark_aligned_data/
    '''


rule sort_bam:
    input:
        'bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.bam'
    output:
        bam='bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.sorted.bam',
        bai='bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.sorted.bam.bai'

    threads: 16
    shell:
        '''
export PATH=/public/home/mowp/softwares/Bismark:$PATH

samtools sort -@ {threads} -o {output.bam} {input}
samtools index -@ {threads} {output.bam}
    '''


rule bismark_bamCoverage:
    input: 
        'bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.sorted.bam'
    output:
        'bismark_aligned_data/bw_covearge/{sample_name}.cov.bw'
    threads: 16
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing None --skipNonCoveredRegions --numberOfProcessors {threads}
        '''


rule MethylDackel_extract:
    input: 
        'bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.sorted.bam'
    output:
        temp('bismark_aligned_data/bw_methy/{sample_name}_CpG.bedGraph'),
        temp('bismark_aligned_data/bw_methy/{sample_name}_CHG.bedGraph'),
        temp('bismark_aligned_data/bw_methy/{sample_name}_CHH.bedGraph'),
        # 'bismark_aligned_data/bw_methy/{sample_name}.cytosine_report.txt'
    threads: 16
    params:
        name = 'bismark_aligned_data/bw_methy/{sample_name}',
        genome = config['genome']
    shell:
        '''
MethylDackel extract -@ {threads} --CHG --CHH -o {params.name} {params.genome} {input}
        '''


rule MethylDackel_bdg_to_bw:
    input: 
        'bismark_aligned_data/bw_methy/{sample_name}_{methy_type}.bedGraph',
    output:
        bdg=temp('bismark_aligned_data/bw_methy/{sample_name}_{methy_type}.bdg'),
        bw='bismark_aligned_data/bw_methy/{sample_name}_{methy_type}.bw',
    threads: 1
    params:
        gsize = config['gsize']
    shell:
        '''
tail -n +2 {input} | awk '{{print$1"\\t"$2"\\t"$3"\\t"$4/100}}' > {output.bdg}
bedGraphToBigWig {output.bdg} {params.gsize} {output.bw}
        '''