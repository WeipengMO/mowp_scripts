#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-01-22 17:36:35
LastEditors  : windz
LastEditTime : 2021-08-03 09:29:48
FilePath     : /tools/metaplot/metaplot.py
'''

import numpy as np
import pyBigWig
import pyranges as pr
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

# 加载画图包
from matplotlib import pyplot as plt


# get last pa sites
# get single polya site gene
last_pa_bed = '/public/home/mowp/workspace/termination/cbRNA_pool/polya_sites/cbRNA.last_polya_cluster_summit.bed'
last_pa = pr.read_bed(last_pa_bed, as_df=True)
last_pa.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ratio']
last_pa.loc[:, 'Chromosome'] = last_pa.loc[:, 'Chromosome'].astype('str')

mask = last_pa.loc[:, 'Name'].str.contains('_1')
single_pa_site_gene = last_pa[mask].loc[:, 'Name'].map(lambda x: x.split('_')[0])

last_pa['Name'] = last_pa['Name'].map(lambda x: x.split('_')[0])
last_pa = last_pa.set_index(['Name'])


# get gene model
gene_model_bed = '/public/home/mowp/db/Arabidopsis_thaliana/isoform/araport11.gene.bed'
gene_model = pr.read_bed(gene_model_bed, as_df=True)
gene_model = gene_model.set_index(['Name'])


# load TSS
tss_bed = '/public/home/mowp/data/public_data/epigentics_data/cage_seq_2020/major_TSS.bed'
tss_bed = pr.read_bed(tss_bed, as_df=True)
tss_bed = tss_bed.set_index(['Name'])


def get_target_site(site_type, gene_id):
    if site_type == 'PAS':
        # 获取基因PAS
        return last_pa.at[gene_id, 'End']
    elif site_type == 'TSS':
        # 获取基因TSS
        # array(['1', 3629, 5899, '.', '+'], dtype=object)
        try:
            values = tss_bed.loc[gene_id, :].values
            if values[4] == '+':
                return values[1]
            else:
                return values[2]
        except KeyError:
            values = gene_model.loc[gene_id, :].values
            if values[4] == '+':
                return values[1]
            else:
                return values[2]
            
    elif site_type == 'TES':
        # 获取基因TES
        # array(['1', 3629, 5899, '.', '+'], dtype=object)
        values = gene_model.loc[gene_id, :].values
        if values[4] == '+':
            return values[2]
        else:
            return values[1]
    else:
        raise KeyError

        
def get_cov(infile, gene_id, before, after):
    '''
    计算sites(加上上下游)区域覆盖度的频数
    '''
    chrom, start, end, _, strand = gene_model.loc[gene_id]
    if chrom in {'Pt', 'Mt'}:
        return None
    
    tss_site = get_target_site('TSS', gene_id)
    pas_site = get_target_site('PAS', gene_id)

    bw = pyBigWig.open(infile)
    try:
        if strand == '+':
            tss_cov = bw.values(chrom, tss_site-before, tss_site+after)
            pas_cov = bw.values(chrom, pas_site-before, pas_site+after)
        else:
            tss_cov = bw.values(chrom, tss_site-after, tss_site+before)[::-1]
            pas_cov = bw.values(chrom, pas_site-after, pas_site+before)[::-1]
    except RuntimeError:
        return
    
    tss_cov = np.nan_to_num(tss_cov)
    pas_cov = np.nan_to_num(pas_cov)
    sum_cov = sum(tss_cov)+sum(pas_cov)
    if sum_cov > 0:
        return tss_cov, pas_cov, sum_cov

    
def get_meta_result(infile, gene_list, before=2000, after=2000, threads=64):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list)/threads)
        results = e.map(get_cov, repeat(infile), gene_list, repeat(before), repeat(after), chunksize=chunksize)
    
    tss_cov = np.zeros(before+after)
    pas_cov = np.zeros(before+after)
    sum_cov = 0

    for res in results:
        if res is not None:
            tss_cov += res[0]
            pas_cov += res[1]
            #sum_cov += res[2]
            sum_cov += 1
            
    return tss_cov, pas_cov, sum_cov


def metaplot(ax, cov, n, before=2000, after=2000, target_site=0, label=None, ylabel=None):
    ax.plot(cov/n, label=label)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.set_xlabel('Relative distance')
    
    ax.set_xticks([0, before, before+after])
    ax.set_xticklabels([f'-{before//1000} Kb', target_site, f'{after//1000} Kb'], rotation=90)
    
    ax.axvline(before, ls='--', color='#555555')
    if label is not None:
        ax.legend(frameon=False)