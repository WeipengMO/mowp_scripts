# 从bam文件中提取甲基化信息
rule extract_methratio:
    input:
        'aligned_data/{sample_name}.bam',
    output:
        out1=temp("methratio/{sample_name}_methratio.txt"),
        out2=temp("methratio/{sample_name}_wiggle.txt"),
        out3="methratio/{sample_name}_methratio.txt.gz",
        out4="methratio/{sample_name}_wiggle.txt.gz",
        sorted_bam='aligned_data/{sample_name}.tmpSrt.bam',
    params:
        genome=config['genome'],
        wig_bin=1,
        min_depth=1
    threads: 32
    shell:
        '''
export PATH=/public/home/mowp/anaconda3/envs/py2/bin/:$PATH
export PATH=/public/home/mowp/softwares/bsmapz/bin/:$PATH

methratio.py -o {output.out1} -d {params.genome} --wig {output.out2} --wig-bin {params.wig_bin} -u -z -r -m {params.min_depth} {input} -N {threads} 

pigz -p 10 -c {output.out1} > {output.out3}
pigz -p 10 -c {output.out2} > {output.out4}
        '''


# 对methratio.py的输出文件进行过滤，去掉CT_count=0的位点，分为CG,CHG和CHH进行保存
rule get_bedgraph:
    input: 
        'methratio/{sample_name}_methratio.txt.gz',
    output:
        temp('bdg_files/{sample_name}.methratio.cg.bdg'),
        temp('bdg_files/{sample_name}.methratio.chg.bdg'),
        temp('bdg_files/{sample_name}.methratio.chh.bdg')
    params:
        output_path='bdg_files'
    shell:
        '''
python script/get_bedgraph.py -i {input} -o {params.output_path}
        '''


# 将bedGraph文件转化成bigwig文件
rule convert_bdg_to_bw:
    input: 
        cg='bdg_files/{sample_name}.methratio.cg.bdg',
        chg='bdg_files/{sample_name}.methratio.chg.bdg',
        chh='bdg_files/{sample_name}.methratio.chh.bdg'
    output:
        cg='bw_files/{sample_name}.methratio.cg.bw',
        chg='bw_files/{sample_name}.methratio.chg.bw',
        chh='bw_files/{sample_name}.methratio.chh.bw'
    params:
        gsize=config['gsize']
    shell:
        '''
bedGraphToBigWig {input.cg} {params.gsize} {output.cg}
bedGraphToBigWig {input.chh} {params.gsize} {output.chh}
bedGraphToBigWig {input.chg} {params.gsize} {output.chg}
        '''