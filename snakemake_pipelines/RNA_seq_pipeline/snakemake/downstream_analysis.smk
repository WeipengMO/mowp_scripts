rule run_stringtie:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        'stringtie/{sample_name}/{sample_name}.gene.gtf'
    threads: 16
    params:
        gff=config['ann_gff'],
        gene_abund='stringtie/{sample_name}/{sample_name}.gene_abund.tab'
    shell:
        '''
stringtie -A {params.gene_abund} -e --rf -B -p {threads} --rf -G {params.gff} -o {output} {input}
        '''

# 提取表达量
rule extract_rpkm:
    input:
        'stringtie/{sample_name}/{sample_name}.gene.gtf'
    output:
        mrna='stringtie/{sample_name}/{sample_name}.RNA.rpkm.txt',
        gene='stringtie/{sample_name}/{sample_name}.gene.rpkm.txt'
    params:
        label='{sample_name}',
        path='stringtie/{sample_name}/'
    threads: 1
    shell:
        '''
export PATH=/public/home/mowp/anaconda3/envs/R/bin/:$PATH
Rscript script/extract_rpkm_from_ballgown.R {params.label} {params.path} {output.mrna} {output.gene}
        '''


# 整合stringtie结果，提取基因/转录本count
# prepDE.py 为stringtie的一个脚本
# /public/home/mowp/softwares/bio/bin/prepDE.py
sample_name = [os.path.basename(fn.split('_1')[0]) for fn in glob('raw_data/*_1.fastq.gz')]
rule extract_gene_count:
    input:
        [f'stringtie/{sample_name}/{sample_name}.RNA.rpkm.txt' for sample_name in sample_name]
    output:
        gene_count_matrix='stringtie/gene_count_matrix.csv',
        transcript_count_matrix='stringtie/transcript_count_matrix.csv'
    threads: 1
    params:
        in_dir='stringtie'
    shell:
        '''
prepDE.py -g {output.gene_count_matrix} -t {output.transcript_count_matrix} -i {params.in_dir}
        '''


rule AScaller:
    input:
        'aligned_data/{sample_name}.sorted.rmdup.bam'
    output:
        'irratio/{sample_name}.irratio.txt'
    threads: 1
    params:
        intron_pos='supplementary_data/intron_pos.repr.bed'
    shell:
        '''
python script/ASCaller.py -i {input} -o {output} --file_intron_pos {params.intron_pos} --strand_flag 0 --min_overlap 6
        '''