#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-01-22 17:36:35
LastEditors  : windz
LastEditTime : 2021-08-04 15:21:15
FilePath     : /tools/metaplot/metaplot.py
'''

import numpy as np
import pyBigWig
import pyranges as pr
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import pysam

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

# For bw file
def get_target_site(site_type: str, gene_id: str) -> int:
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

        
def get_bw_cov(infile: str, gene_id: str, before: int, after: int):
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

    
def get_bw_meta_result(infile, gene_list, before=2000, after=2000, threads=64):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list)/threads)
        results = e.map(get_bw_cov, repeat(infile), gene_list, repeat(before), repeat(after), chunksize=chunksize)
    
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


# For bam file
STRAND_TO_BOOL = {'-': True, '+': False}

def get_bam_cov(infile: str, gene_id: str, before: int, after: int):
    """
    BAM file for tagged FLEP-seq data
    Ignore splicing junction
    """
    chrom, *_, strand = gene_model.loc[gene_id]
    strand_boo = STRAND_TO_BOOL[strand]
    if chrom in {'Pt', 'Mt'}:
        return None
    
    n = 0
    cov_list = []
    read_set = set()
    for site_type in ['TSS', 'PAS']:
        target_site = get_target_site(site_type, gene_id)
        cov = np.zeros(before+after)

        if strand == '+':
            start = target_site-before
            end = target_site+after
        else:
            start = target_site-after
            end = target_site+before
        
        if start <= 0:
            return

        with pysam.AlignmentFile(infile, 'rb') as inbam:
            for read in inbam.fetch(chrom, start, end):
                # 判断是否跟基因是同个方向，针对于链特异文库
                read_strand = read.is_reverse
                if strand_boo is not read_strand:
                    continue
                
                read_gene_id = read.get_tag('gi')
                if read_gene_id not in {gene_id, 'None'}:
                    continue
                    
                if strand == '+':
                    read_five_end = read.reference_start
                    read_three_end = read.reference_end
                    cov_start = read_five_end-start if read_five_end-start >= 0 else 0
                    cov_end = read_three_end-start if read_three_end-start <= before+after else end-start
                else:
                    read_five_end = read.reference_end
                    read_three_end = read.reference_start
                    cov_start = end-read_five_end if end-read_five_end >= 0 else 0
                    cov_end = end-read_three_end if end-read_three_end <= before+after else end-start

                cov[cov_start: cov_end] += 1
                
                if read.query_name not in read_set:
                    n += 1
                    read_set.add(read.query_name)
        
        cov_list.append(cov)
    
    if max(cov_list[1]) > 20:
        # 只返回与PAS有>20条reads overlap的基因
        return cov_list[0], cov_list[1], n, gene_id


def get_bam_meta_result(infile, gene_list, before=2000, after=2000, threads=64):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list)/threads)
        results = e.map(get_bam_cov, repeat(infile), gene_list, repeat(before), repeat(after), chunksize=chunksize)
    
    tss_cov, pas_cov = np.zeros(before+after), np.zeros(before+after)
    n = 0
    for res in results:
        if res is not None:
            tss_cov += res[0]/res[2]
            pas_cov += res[1]/res[2]
            n += 1
            
    return tss_cov, pas_cov, n


def get_bam_total_readcounts(infile: str):
    """
    获取bam文件中所有reads数
    """
    return eval('+'.join([line.split('\t')[2] for line in pysam.idxstats(infile).rstrip().split('\n')]))



def plot(ax, cov, n, before=2000, after=2000, target_site=0, label=None, ylabel=None):
    """
    画metaplot
    """
    ax.plot(cov/n, label=label)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_xticks([0, before, before+after])
    ax.set_xticklabels([f'-{before//1000} Kb', target_site, f'{after//1000} Kb'], rotation=90)
    
    ax.axvline(before, ls='--', color='#555555')
    if label is not None:
        ax.legend(frameon=False)