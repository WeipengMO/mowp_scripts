# fork from Yiming, not test!

# 对methratio.py的输出文件进行过滤，去掉CT_count=0的位点，分为CG,CHG和CHH进行保存
chr_ID = ['1', '2', '3', '4', '5']
meth_type = ['cg','chg','chh']
rule filter_to_plot:
    input:
        "methratio/{sample_name}_methratio.txt.gz"
    output:
        out1="methratio/{sample_name}_methratio_filter.txt.gz",
        out2="methratio/{sample_name}_methratio_filter.cg.gz",
        out3="methratio/{sample_name}_methratio_filter.chg.gz",
        out4="methratio/{sample_name}_methratio_filter.chh.gz",

        out5=expand("methratio/{{sample_name}}_methratio_filter.{methy}.{chr_id}", chr_id=chr_ID, methy=meth_type),
        out6=expand("DMR/input/{{sample_name}}_methratio_filter.{methy}.{chr_id}.4DMR", chr_id=chr_ID, methy=meth_type)
    threads: 1
    shell:
        '''
python script/filter_to_plot.py --input {input} --out_total {output.out1} --out_CG {output.out2} --out_CHG {output.out3} --out_CHH {output.out4} --chrMCid Pt Mt
        '''


# DMR的一些参数
mutant = config['mutant']
wildtype = config['wildtype']
Bin_size1=[200]

Diff_cg=70
Diff_chg=50
Diff_chh=10

# CG DMR calling
rule CG_DMR_calling:
    input:
        int1="/DMR/input/{mutant}_methratio_fil.cg.{chr_id}.4DMR",
        int2="/DMR/input/{widetype}_methratio_fil.cg.{chr_id}.4DMR"
    output:
        out1=expand("/DMR/output/{bin_size}/{{mutant}}_vs_{{widetype}}.cg.{{chr_id}}.4DMR.obj.RData", bin_size=Bin_size1),
        out2=expand("/DMR/output/{bin_size}/{{mutant}}_vs_{{widetype}}.cg.{{chr_id}}.4DMR.fil.obj.RData", bin_size=Bin_size1),
        out3=expand("/DMR/output/{bin_size}/{{mutant}}_vs_{{widetype}}.cg.{{chr_id}}.4DMR.{bin_size}.Diff{diff_cg}p.RData", bin_size=Bin_size1, diff_cg=Diff_cg)
    params:
        bin_size=Bin_size1[0],
        min_cov=4,
        assembly="tair",  ## 没什么用处
        context='cg',
        diff_meth=Diff_cg
    threads: 1
    shell:
        '''
R --slave --args {params.bin_size} {params.min_cov} {input.int1} {input.int2} {params.assembly} {params.context} {params.diff_meth} {output.out1} {output.out2} {output.out3} < script/methylKit_DMR_v2.R
        '''