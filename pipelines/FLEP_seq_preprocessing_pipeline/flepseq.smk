'''
Author: 
    windz
Usage:
    snakemake -j 15 --cluster 'bsub -J {rulename} -n {threads} -o log/%J.stdout -e log/%J.stderr -R span[hosts=1]' -s flepseq.smk
'''

# create log dir
if not os.path.exists('log'):
    os.mkdir('log')


configfile: 'config.yml'
sample=config['sample']
genome_data=config['genome']


rule all:
    input:
        expand('results/{sample}.read.info.txt', sample=sample),
        expand('results/{sample}.read.splicing_kinetics.txt', sample=sample),
        expand('results/{sample}.read.rna.ir.stat', sample=sample),


rule fastq_to_fasta:
    input:
        'guppy_out/pass/',
    output:
        'basecalled_data/{sample}.fasta'
    threads: 1
    shell:
        '''
python script/fastqdir2fasta.py --indir {input} --out {output}
        '''


rule mapping_to_genome:
    input:
        'basecalled_data/{sample}.fasta'
    output:
        bam='aligned_data/{sample}.sorted.bam'
    params:
        genome=genome_data
    threads: 36
    shell:
        '''
minimap2 -t {threads} -ax splice --secondary=no -G 12000 {params.genome} {input} | samtools sort -@ {threads} -o {output.bam} -
samtools index -@ {threads} {output.bam}
        '''


rule find_3linker:
    input:
        bam='aligned_data/{sample}.sorted.bam',
        fasta='basecalled_data/{sample}.fasta'
    output:
        'aligned_data/{sample}.adapter.result.txt'
    threads: 36
    shell:
        '''
python script/adapterFinder.py --inbam {input.bam} --inseq {input.fasta} --out {output} --threads {threads} --mode 1
        '''


rule polyacaller:
    input:
        adapter_result='aligned_data/{sample}.adapter.result.txt',
        sequencing_summary='sequencing_summary.txt',
        fast5_dir='guppy_out/workspace',
    output:
        'aligned_data/{sample}.polyA_tail.result.txt'
    threads: 36
    shell:
        '''
python script/PolyACaller.py --inadapter {input.adapter_result} --summary {input.sequencing_summary}  --fast5dir {input.fast5_dir} --out {output} --threads {threads}
        '''


rule identify_read_info:
    input:
        'aligned_data/{sample}.sorted.bam'
    output:
        'aligned_data/{sample}.read_info.result.txt'
    params:
        'genome_data/exon_intron_pos.repr.bed'
    threads: 1
    shell:
        '''
python script/extract_read_info.py --inbed {params} --inbam {input} --out {output}
        '''


rule merge_info:
    input:
        read_info='aligned_data/{sample}.read_info.result.txt',
        adapter='aligned_data/{sample}.adapter.result.txt',
        polya='aligned_data/{sample}.polyA_tail.result.txt'
    output:
        read_info_result='results/{sample}.read.info.txt'
    threads: 1
    params:
        'Nanopore'
    shell:
        '''
export PATH=/scem/work/mowp/anaconda3/envs/R/bin/:$PATH
Rscript script/merge_read_info.R --type {params} --inreadinfo {input.read_info} --inadapter {input.adapter} --inpolya {input.polya} --out {output.read_info_result}
        '''


rule splicing_kinetics:
    input:
        read_info='results/{sample}.read.info.txt'
    output:
        splicing_data='results/{sample}.read.intron.pos.splicing.txt',
        splicing_kinetics='results/{sample}.read.splicing_kinetics.txt',
        figure='results/{sample}.read.splicing_kinetics.pdf'
    threads: 1
    params:
        inbed='genome_data/exon_intron_pos.repr.bed',
        select_intron='genome_data/select_introns.txt'
    shell:
        '''
python script/prepare_data_for_splice_kinetics.py --inreadinfo {input.read_info} --inbed {params.inbed} --out {output.splicing_data}
export PATH=/scem/work/mowp/anaconda3/envs/R/bin/:$PATH
Rscript script/plot_intron_splicing_kinetics.R --inrelpos {output.splicing_data} --inreadinfo {input.read_info} --inintron {params.select_intron} --out {output.splicing_kinetics} --pdf {output.figure}
        '''


# Calculate intron retention ratio of polyadenylated transcripts
rule intron_retention_ratio:
    input:
        splicing_data='results/{sample}.read.intron.pos.splicing.txt',
        read_info='results/{sample}.read.info.txt',
    output:
        rna_ir='results/{sample}.read.rna.ir.stat',
        intron_ir='results/{sample}.read.intron.ir.stat',
    threads: 1
    shell:
        '''
export PATH=/scem/work/mowp/anaconda3/envs/R/bin/:$PATH
Rscript script/cal_polya_transcript_ir.R --inrelpos {input.splicing_data} --inreadinfo {input.read_info} --outrna {output.rna_ir} --outintron {output.intron_ir}
        '''

