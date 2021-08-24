#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-08-24 16:29:53
LastEditTime : 2021-08-24 19:24:15
LastEditors  : windz
FilePath     : /public/home/mowp/test/BS_Seq/script/filter_to_plot.py
'''

import argparse
import gzip
import re

def main():
    description="""
    input: 输入文件 Nip_methratio.txt.gz，脚本methratio.py的输出结果
    out_total：输出文件1，对输入文件Nip_methratio.txt.gz进行整理，输出，'chrBase','chr','base','strand','coverage','freqC','freqT'
    out_CG, out_CHG, out_CHH: 将输出文件1 out_total 安装CG，CHG和CHH进行分割
    """
    parser = argparse.ArgumentParser(description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/Nip_methratio.txt.gz")
    parser.add_argument("--out_total", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/Nip_methratio_fil.txt.gz")
    parser.add_argument("--out_CG", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/Nip_methratio_fil.cg.gz")
    parser.add_argument("--out_CHG", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/Nip_methratio_fil.chg.gz")
    parser.add_argument("--out_CHH", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/Nip_methratio_fil.chh.gz")
    parser.add_argument("--chrMCid", nargs='+', default="chrC chrM") ## 直接传入列表 ['chrC','chrM']
    #parser.add_argument("--out_chrid", default="/public/home/yuym/taiyi/data/BS_Seq/rice/all_T1_4DMR_20190411/methratio/chrID")

    args = parser.parse_args()

    input_ = args.input
    out_total = args.out_total
    out_CG = args.out_CG
    out_CHG = args.out_CHG
    out_CHH = args.out_CHH
    chrMCid = args.chrMCid

    def read_raw_meth(input_, out_total, out_CG, out_CHG, out_CHH, chrMCid):
        chr_ = set()
        with gzip.open(input_, 'rt') as f, gzip.open(out_total, 'wt') as w, gzip.open(out_CG, 'wt') as w1, gzip.open(out_CHG, 'wt') as w2, gzip.open(out_CHH, 'wt') as w3:
            title=['chrBase','chr','base','strand','coverage','freqC','freqT']
            w.write('\t'.join(title)+'\n')
            for i in f:
                i = i.rstrip('\n').split('\t')
                if i[0] == "chr": continue
                chr_.add(i[0])
                if i[7] != 0:
                    chrBase = f"{i[0]}.{i[1]}.{i[3]}"
                    T_ct = str(int(i[7]) - int(i[6]))
                    w.write(f"{chrBase}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[7]}\t{i[6]}\t{T_ct}\n")
                    if i[3] == "CG":
                        w1.write(f"{chrBase}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[7]}\t{i[6]}\t{T_ct}\n")
                    elif i[3] == "CHG":
                        w2.write(f"{chrBase}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[7]}\t{i[6]}\t{T_ct}\n")
                    elif i[3] == "CHH":
                        w3.write(f"{chrBase}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[7]}\t{i[6]}\t{T_ct}\n")
        for id_ in chrMCid:
            if id_ in chr_: chr_.remove(id_)
        return(chr_)

    def split_chr_file(chr_id, out_file):
        meth_dict = {}
        with gzip.open(out_file, 'rt') as f:
            for i in f:
                i = i.rstrip('\n')
                chr_ = i.split('\t')[1]
                try:
                    meth_dict[chr_].append(i)
                except:
                    meth_dict[chr_] = [i]
        for x in chr_id:
            file_ = re.sub(".gz", f".{x}", out_file)
            file_dmr = f"{file_}.4DMR"
            file_dmr = re.sub('methratio/','DMR/input/',file_dmr)
            with open(file_, 'wt') as w1, open(file_dmr, 'wt') as w2:
                w2.write(f"chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
                for y in meth_dict[x]:
                    w1.write(f"{y}\n")
                    m = y.split('\t')
                    if m[3] == "+":
                        m[3] = "F"
                    elif m[3] == "-":
                        m[3] = "R"
                    m[5] = str(round(int(m[5])/int(m[4])*100,2))
                    m[6] = str(round(int(m[6])/int(m[4])*100,2))
                    w2.write('\t'.join(m)+'\n')

    chr_id = read_raw_meth(input_, out_total, out_CG, out_CHG, out_CHH, chrMCid)

    split_chr_file(chr_id, out_CG)
    split_chr_file(chr_id, out_CHG)
    split_chr_file(chr_id, out_CHH)
    
if __name__ == "__main__":
    main()