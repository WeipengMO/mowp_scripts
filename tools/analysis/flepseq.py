
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-20 18:03:56
LastEditTime : 2022-03-20 19:22:08
LastEditors  : windz
FilePath     : flepseq.py
'''

import pysam
import numpy as np
import pandas as pd

def get_read_polya_len(infile: str, min_len: int = 15):
    '''
    Get read polyA length from BAM file
    Args:
        infile: path to BAM file
        min_len: minimum of the polyA tail length
    
    Return:
        polya_len_list: polya length of all reads in np.array
    '''
    polya_len_list = []
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            # require unique mapped read
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            polya_len = read.get_tag('pa')
            if polya_len >= min_len:
                polya_len_list.append(polya_len)
        
        polya_len_list = np.array(polya_len_list)

    return polya_len_list


def get_gene_polya_len(infile, min_len=15, min_count=10, exclude_chr=('Pt', 'Mt')):
    '''
    Get gene median polyA length
    Args:
        infile: path to BAM file
        min_len: minimum of the polyA tail length
        min_count: minimum read counts support the gene polyA length
        exclude_chr: the chromosomes to be exclude
    
    Return:
        gene_polya_len: gene polya length in np.array
                        (('gene_id', length), ...)
        gene_polya_len_raw: raw polyA length of the gene, dict
    '''
    from collections import defaultdict

    exclude_chr = set(exclude_chr)
    gene_polya_len_raw = defaultdict(lambda : []) 
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            if read.reference_name in exclude_chr:
                continue
                
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
                
            polya_len = read.get_tag('pa')
            gene_id = read.get_tag('gi')
            if gene_id == 'None':
                continue
                
            if polya_len >= min_len:
                gene_polya_len_raw[gene_id].append(polya_len)
    
    gene_polya_len = []
    for gene_id in gene_polya_len_raw:
        if len(gene_polya_len_raw[gene_id]) >= min_count:
            gene_polya_len.append(np.median(gene_polya_len_raw[gene_id]))
    
    gene_polya_len = np.array(gene_polya_len)
            
    return gene_polya_len, gene_polya_len_raw


def get_idxstats(inbam: str):
    '''
    Get total read number from BAM file
    Args:
        inbam: path to BAM file
    
    Return:
        total_count: total count in BAM file
        mapped_count: mapped count in BAM file
        unmapped_count: unmapped count in BAM file
        count_per_chr: mapped count per chromosomes
    '''
    count_per_chr = {}
    total_count, mapped_count, unmapped_count = 0, 0, 0
    for l in pysam.idxstats(inbam).split('\n')[:-1]:
        l = l.rstrip().split('\t')
        total_count += eval('+'.join(l[2:]))
        if l[0] != '*':
            count_per_chr[l[0]] = int(l[2])
            mapped_count += int(l[2])
        else:
            unmapped_count += int(l[3])
    
    return total_count, mapped_count, unmapped_count, count_per_chr