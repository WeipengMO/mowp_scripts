'''
Author: 
    windz
Usage:
    snakemake -j 15 --cluster 'bsub -J {rulename} -n {threads} -gpu "num=2" -o %J.stdout -e %J.stderr -R span[hosts=1]' -s basecalling.smk
'''

configfile: 'config.yml'

split_num = config['split_num']
dirname = [i for i in range(split_num)]

sample = config['sample']
fast5_path = config['fast5_path']

rule all:
    input:
        expand('guppy_out/{dirname}/sequencing_summary.txt', dirname=dirname),
        'sequencing_summary.txt',


rule split_fast5_files:
    input:
        fast5_path
    output:
        f'{fast5_path}/split_files_into_dirs.log'
    threads: 1
    params:
        config['split_num']
    shell:
        '''
python script/split_files_into_dirs.py --dir_path {input} --suffix fast5 --split_num {params}
        '''


rule run_guppy:
    input:
        f'{fast5_path}/split_files_into_dirs.log',  # log文件，只是用来等上一步完成
    output:
        'guppy_out/{dirname}/sequencing_summary.txt'
    threads: 9
    params:
        fast5_path=config['fast5_path'],
        dir='{dirname}',
        out='guppy_out/{dirname}/',
        model=config['model'],
        cuda=config['cuda']
    shell:
        '''
guppy_basecaller -i {params.fast5_path}/{params.dir} -s {params.out} -c {params.model} --recursive --fast5_out --disable_pings --qscore_filtering --device {params.cuda}
        '''


rule merge_sequencing_summary:
    input:
        expand('guppy_out/{dirname}/sequencing_summary.txt', dirname=dirname)
    output:
        'sequencing_summary.txt'
    shell:
        '''
head -n 1 {input[0]} > {output}
ls {input} | while read i
do
tail -n -1 $i >> {output}
done
        '''