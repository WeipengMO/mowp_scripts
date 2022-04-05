
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-20 18:03:56
LastEditTime : 2022-03-28 19:34:58
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


def get_gene_polya_len_with_splice_info(
    infile: str, 
    min_len: int = 15, 
    exclude_chr: tuple = ('Pt', 'Mt'),
):
    '''
    Get spliced and unspliced reads polyA length
    Args:
        infile: path to BAM file
        min_len: minimum read count to support the median polyA length
        exclude_chr: the chromosomes to be exclude
    
    Return:
        splice_res, incompletely_splice_res: gene polya length info with splice info
    '''
    from collections import defaultdict


    exclude_chr = set(exclude_chr)
    gene_polya_len_unspliced = defaultdict(lambda : [])
    gene_polya_len_spliced = defaultdict(lambda : []) 
    
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
                if read.get_tag('rn') == 0:
                    gene_polya_len_spliced[gene_id].append(polya_len)
                else:
                    gene_polya_len_unspliced[gene_id].append(polya_len)
                    
    spliced_res, incompletely_spliced_res = {}, {}
    for k, v in gene_polya_len_spliced.items():
        if len(v) >= min_len:
            spliced_res[k] = np.median(v)
            
    for k, v in gene_polya_len_unspliced.items():
        if len(v) >= min_len:
            incompletely_spliced_res[k] = np.median(v)
            
    return spliced_res, incompletely_spliced_res

    
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


def get_gene_counts(
    bed_intersect: str = None, 
    inbam : str = None, repr_gene : str = None, 
    exclude_chr: set = {'Pt', 'Mt'}
    ) -> dict:
    '''
    计算有reads覆盖的基因个数
    Args:
        bed_intersect: bamfile intersect with gene model
            # bedtools intersect -abam {input} -b {params.repr_gene} -s -wo -split -bed > {output}

        **如果没有提供bed_intersect结果，可以调用pybedtools进行计算，二选一**
        bam_file: path to BAM file
        repr_gene: gene model in BED6 format
        exclude_chr: chromosomes to be excluded
    
    Return:
        gene_counts: read count per gene, Counter
    '''
    from collections import Counter
    from utils.get_overlap_genes import exclude_genes

    
    NAMES=[
    'Chromosome', 'Start', 'End', 'read_id', 'Score', 'Strand', 
    'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts', 
    'geneChromosome', 'geneStart', 'geneEnd', 'gene_id', 'geneScore', 'geneStrand', 
    'cov'
    ]

    USECOLS = [
        'Chromosome', 'Start', 'End', 'read_id', 'Strand',
        'geneStart', 'geneEnd', 'gene_id', 'geneStrand', 'cov',
    ]
    
    if bed_intersect is not None:
        df = pd.read_csv(
            bed_intersect, 
            sep='\t', 
            names=NAMES,
            usecols=USECOLS,
            header=None,
            dtype={"Chromosome": str}
            )
    elif inbam is not None and repr_gene is not None:
        from pybedtools import BedTool

        bam = BedTool(inbam)
        if not bam._isbam:
            raise ValueError('inbam is not BAM file')
        res = bam.intersect(repr_gene, s=True, wo=True, split=True, bed=True)
        df = res.to_dataframe(
            disable_auto_names=True, 
            names=NAMES, 
            usecols=USECOLS, 
            dtype={"Chromosome": str}
            )
    else:
        raise ValueError('Please no input found')

    gene_counts = Counter()
    for item in df.itertuples():
        if item.Chromosome in exclude_chr:
            continue
        if item.Strand != item.geneStrand:
            continue
        if item.gene_id in exclude_genes:
            continue

        if item.Strand == '+' and item.Start >= item.geneStart-100:
            gene_counts[item.gene_id] += 1
        elif item.Strand == '-' and item.End <= item.geneEnd+100:
            gene_counts[item.gene_id] += 1
    
    return gene_counts