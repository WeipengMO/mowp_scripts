reference = config['genome']

rule cuteSV:
    input:
        'aligned_data/{sample}.sorted.bam'
    output:
        'cuteSV/{sample}.cuteSV.vcf'
    threads: 48
    conda:
        'genome'
    params:
        min_sv_length=30,  # Minimum size of SV to be reported
        max_sv_length=1000000,  # Maximum size of SV to be reported. All SVs are reported when using -1
        min_read_length=10000,  # Ignores reads that only report alignments with not longer than bp
        min_read_mapping_quality=20,
        min_read_support=3,
        output='{sample}.cuteSV.vcf',
        workdir='cuteSV',
    shell:
        '''
cuteSV \
    --threads {threads} \
    --report_readid \
    --genotype \
    -l {params.min_sv_length} \
    -L {params.max_sv_length} \
    -r {params.min_read_length} \
    -q {params.min_read_mapping_quality} \
    -s {params.min_read_support} \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 \
    {input} \
    {reference} \
    {output} \
    .
        '''


rule mosdepth:
    input:
        bam='aligned_data/{sample}.sorted.bam'
    output:
        'cuteSV/{sample}.mosdepth.global.dist.txt'
    threads: 48
    params:
        prefix='cuteSV/{sample}'
    conda:
        'genome'
    shell:
        '''
mosdepth \
    -n \
    -x \
    -t {threads} \
    {params.prefix} \
    {input.bam}
        '''
# -b --by <bed|window>       optional BED file or (integer) window-sizes.


rule filterCalls:
    input:
        depth_bedfile='cuteSV/{sample}.mosdepth.global.dist.txt',
        vcf='cuteSV/{sample}.cuteSV.vcf'
    output:
        '{sample}.filter.sh'
    params:
        prefix='{sample}',
        min_sv_length = 30,
        max_sv_length = 1000000,
        min_read_length = 1000,
        min_read_mapping_quality = 20,
        min_read_support = 3,
        min_read_support_limit=5,
        sv_types='INS DEL DUP INV BND'
    shell:
        '''
python scripts/get_filter_calls_command.py \
    --vcf {input.vcf} \
    --depth_bedfile {input.depth_bedfile} \
    --min_sv_length {params.min_sv_length} \
    --max_sv_length {params.max_sv_length} \
    --min_read_support {params.min_read_support} \
    --min_read_support_limit {params.min_read_support_limit} \
    --sv_types {params.sv_types} > {output}
        '''