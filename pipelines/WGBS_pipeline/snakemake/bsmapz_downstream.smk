rule bsmapz_sort_bam:
    input:
        'bsmapz_aligned_data/{sample_name}.bam'
    output:
        bam='bsmapz_aligned_data/{sample_name}.sorted.bam',
        bai='bsmapz_aligned_data/{sample_name}.sorted.bam.bai'

    threads: 16
    shell:
        '''
samtools sort -@ {threads} -o {output.bam} {input}
samtools index -@ {threads} {output.bam}
    '''


rule bsmapz_MethylDackel_extract:
    input: 
        'bsmapz_aligned_data/{sample_name}.sorted.bam'
    output:
        temp('bsmapz_aligned_data/bw_methy/{sample_name}_CpG.bedGraph'),
        temp('bsmapz_aligned_data/bw_methy/{sample_name}_CHG.bedGraph'),
        temp('bsmapz_aligned_data/bw_methy/{sample_name}_CHH.bedGraph'),
    threads: 16
    params:
        name = 'bsmapz_aligned_data/bw_methy/{sample_name}',
        genome = config['genome']
    shell:
        '''
MethylDackel extract -@ {threads} --CHG --CHH -o {params.name} {params.genome} {input}
        '''


rule bsmapz_MethylDackel_bdg_to_bw:
    input: 
        'bsmapz_aligned_data/bw_methy/{sample_name}_{methy_type}.bedGraph',
    output:
        bdg=temp('bsmapz_aligned_data/bw_methy/{sample_name}_{methy_type}.bdg'),
        bw='bsmapz_aligned_data/bw_methy/{sample_name}_{methy_type}.bw',
    threads: 1
    params:
        gsize = config['gsize']
    shell:
        '''
tail -n +2 {input} | awk '{{print$1"\\t"$2"\\t"$3"\\t"$4/100}}' > {output.bdg}
bedGraphToBigWig {output.bdg} {params.gsize} {output.bw}
        '''