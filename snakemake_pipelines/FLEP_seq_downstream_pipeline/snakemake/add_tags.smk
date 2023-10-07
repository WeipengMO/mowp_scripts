rule run_minimap2:
    input:
        reads='basecalled_data/{sample_name}.fastq.gz',
    output:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai'
    threads: 58
    params:
        ref=config['genome'],
    shell:
        '''
minimap2 -t {threads} -ax splice -G 12000 --secondary=no {params.ref} {input.reads} | samtools view -@ {threads} -F 2308 -hb - | samtools sort -@ {threads} -O bam -o {output.bam} -
samtools index -@ {threads} {output.bam}
        '''


rule add_tags_to_bam:
    input:
        infile='aligned_data/{sample_name}.sorted.bam',
        read_info='results/{sample_name}.read_info.result.txt',
        adapter_info='results/{sample_name}.adapter.result.txt',
        polya_info='results/{sample_name}.polyA_tail.result.txt'
    output:
        bam='aligned_data/{sample_name}.sorted.tagged.bam',
        bai='aligned_data/{sample_name}.sorted.tagged.bam.bai'
    threads: 1
    shell:
        '''
python script/add_tag_to_bam.py -i {input.infile} --read_info {input.read_info} --adapter_info {input.adapter_info} --polya_info {input.polya_info} -o {output.bam}
        '''


rule get_elongating_reads:
    input:
        'aligned_data/{sample_name}.sorted.tagged.bam'
    output:
        bam='elongating_data/{sample_name}.elongating.bam',
        bai='elongating_data/{sample_name}.elongating.bam.bai'
    threads: 1
    shell:
        '''
python script/get_elongating_reads.py -i {input} -o {output.bam}
samtools index -@10 {output.bam}
        '''


rule get_polyadenylated_reads:
    input:
        'aligned_data/{sample_name}.sorted.tagged.bam'
    output:
        bam='polyadenylated_data/{sample_name}.polyadenylated.bam',
        bai='polyadenylated_data/{sample_name}.polyadenylated.bam.bai'
    threads: 1
    shell:
        '''
python script/get_polyadenylated_reads.py -i {input} -o {output.bam}
samtools index -@10 {output.bam}
        '''


