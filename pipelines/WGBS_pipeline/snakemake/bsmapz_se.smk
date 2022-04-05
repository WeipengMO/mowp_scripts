rule run_bsmapz_se:
    input:
        fq1='raw_data/{sample_name}.fastq.gz'
    output:
        'aligned_data/{sample_name}.bam',
    params:
        genome=config['genome'],
        mismatch=0.08,
        base_quality=33,  # quality coding format
        gap_size=0,
        best_hits_num=1,
        report_repeat=0,
        lib_type=0  # #-n 0 对应Lister建库方式（比对到两条正链，notags）；-n 1 对应Cokus建库方式（比对到所有的四种链，tags）
    threads: 32
    shell:
        '''
export PATH=/public/home/mowp/anaconda3/envs/py2/bin/:$PATH
export PATH=/public/home/mowp/softwares/bsmapz/bin/:$PATH

bsmapz -a {input.fq1} -d {params.genome} -p {threads} -v {params.mismatch} -z {params.base_quality} -n {params.lib_type} -u -g {params.gap_size} -w {params.best_hits_num} -r {params.report_repeat} -o {output}
    '''