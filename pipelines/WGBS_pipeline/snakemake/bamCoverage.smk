rule bamCoverage:
    input: 
        'aligned_data/{sample_name}.tmpSrt.bam',
    output:
        'bw_coverage/{sample_name}.{strand}.bw',
    threads: 16
    params:
        strand = lambda wildcard: 'forward' if wildcard.strand == 'reverse' else 'reverse'
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing None --filterRNAstrand {params.strand} --skipNonCoveredRegions --numberOfProcessors {threads}
        '''


rule bamCoverage_bismark:
    input: 
        'bismark_aligned_data/{sample_name}_bismark_bt2.deduplicated.sorted.bam',
    output:
        'bismark_aligned_data/bw_coverage/{sample_name}.{strand}.bw',
    threads: 16
    params:
        strand = lambda wildcard: 'forward' if wildcard.strand == 'reverse' else 'reverse'
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing None --filterRNAstrand {params.strand} --skipNonCoveredRegions --numberOfProcessors {threads}
        '''