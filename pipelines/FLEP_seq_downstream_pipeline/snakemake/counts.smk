rule run_stringtie:
    input:
        'aligned_data/{sample_name}.sorted.tagged.bam'
    output:
        gene_abund='stringtie/{sample_name}.gene_abund.tsv',
        gtf='stringtie/{sample_name}.gtf'
    threads: 32
    params:
        '/public/home/mowp/db/Arabidopsis_thaliana/gff3/Araport11_GFF3_genes_transposons.201606.gff'
    shell:
        '''
stringtie -A {output.gene_abund} -L -e --rf -B -p {threads} -G {params} -o {output.gtf} {input}
        '''


rule featureCounts:
    input:
        'aligned_data/{sample_name}.sorted.tagged.bam'
    output:
        'featureCounts/{sample_name}.counts_genes.txt'
    threads: 1
    params:
        gtf=config['gtf']
    shell:
        '''
featureCounts -T {threads} --primary -F GTF -g gene_id -t gene -a {params.gtf} -o {output} -L {input}
        '''


rule bamCoverage_strand:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        'bw_files/{sample_name}.sorted.rmdup.{strand}.bw'
    threads: 16
    params:
        strand=lambda wildcard: 'reverse' if wildcard.strand == 'fwd' else 'forward' 
    shell:
        '''
bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors {threads} --filterRNAstrand reverse
        '''