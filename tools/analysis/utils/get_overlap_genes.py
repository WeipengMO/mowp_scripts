#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-23 14:45:12
LastEditTime : 2022-03-23 21:23:44
LastEditors  : windz
FilePath     : get_overlap_genes.py
'''

import pandas as pd


def get_overlap_genes(gene_info: str) -> set:
    '''
    Get overlap genes set and non-overlap gene counts

    Gene1 ---------------------------- (keep)
    Gene2       -----------            (exclude)

    Args:
        gene_info: gene model intersect with itself
    
    Return:
        overlap_gene: gene that inside other gene
        total_gene_counts: non-overlap gene counts
    '''

    # non-overlap gene set
    gene_intersect = pd.read_csv(
        gene_info, sep='\t', 
        names = [
            'chrom1', 'start1', 'end1', 'gene_id1', 'score1', 'strand1',
            'chrom2', 'start2', 'end2', 'gene_id2', 'score2', 'strand2', 'cov' 
        ]
        )

    total_gene_counts = len(set(gene_intersect['gene_id1']))
    gene_intersect.query('gene_id1 != gene_id2', inplace=True)
    
    # 包含在其他基因内部的同向基因
    overlap_gene = set()
    for item in gene_intersect.query('start1 <= start2 and end1 >= end2').itertuples():
        overlap_gene.add(item.gene_id2)

    total_gene_counts += -len(overlap_gene)  # gff文件中所有基因的个数，减去overlap基因个数
    return overlap_gene, total_gene_counts


gene_info = '/public1/mowp/workspace/cstf_mutant/supplementary_data/gene.intersect'
exclude_genes, total_gene_counts = get_overlap_genes(gene_info)