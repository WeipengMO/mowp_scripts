#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2021-12-24 15:09:14
LastEditTime : 2021-12-24 15:42:26
LastEditors  : windz
FilePath     : /atx1_mutant/scripts/get_splice_info.py
Usage        : python script/get_splice_info.py --inbam elongating_data/20211014_atx1_seedling_part.elongating.bam --outbam 20211014_atx1_seedling_part.test.bam --repr_gene /public/home/mowp/db/Arabidopsis_thaliana/representative_gene_model/araport11.representative.gene_model.bed
'''


import numpy as np
import pybedtools
from collections import namedtuple
import pysam
from tqdm import tqdm
import click
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s: %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p')


def is_overlap(block, intron):
    minn = max(block[0], intron[0])
    maxn = min(block[1], intron[1])
    if intron[1]-intron[0] < 0:
        return False
        
    if (maxn-minn)/(intron[1]-intron[0]) >= .5:
        return True
    else:
        return False
        

def get_splice_stats(blocks, introns):
    '''
    计算reads中没有剪切的intron个数
    返回:
        没有被剪切的intron个数
        跨过的总intron个数
    '''
    read_start = blocks[0][0]
    read_end = blocks[-1][-1]

    block_total_count = len(blocks)
    blocks = iter(blocks)
    introns = iter(introns)
    
    unsplice_intron = []  # 用于储存未剪切intron编号
    
    block_count, unsplice_count, span_intron_count, intron_number = 0, 0, 0, 0
    block = next(blocks)
    block_count += 1
    intron = next(introns)
    intron_number += 1  # 当前intron序号
    span_intron = set()

    while True:
        '''
        忽略read没有跨过intron的情况
        read: _________________________
        gene: ________....._______..........._______
        ''' 
        if block[1] > intron[1] and span_intron_count < intron_number:
            span_intron_count = intron_number  # read 当前跨过的intron_id
            
        if not(block_count == block_total_count and block[1] < intron[1]) and is_overlap(block, intron):
            unsplice_count += 1
            unsplice_intron.append(span_intron_count)

        # if read_start < intron[0] and read_end > intron[1]:
        if read_start < intron[1] and read_end > intron[1]:
            span_intron.add(intron_number)
        try:
            if intron[0] >= block[1]:
                block = next(blocks)
                block_count += 1
            elif intron[0] <= block[1]:
                intron = next(introns)
                intron_number += 1
        except StopIteration:
            break
    
    span_intron = sorted(span_intron)
    return unsplice_count, span_intron, unsplice_intron


def get_read_info_dict(inbam: str, repr_gene: str):
    '''Get splice info from bed intersect results

    Args: 
        inbam: input BAM file
        repr_gene: representative gene model, BED12 file
    
    Returns: 
        read_info_dict: 
            key: read_id,
            value: cov, gene_id, unsplice_count, span_intron, unsplice_intron, relative_pos_to_3ss
    '''

    logging.info('[Start] bedtools intersect')
    a = pybedtools.BedTool(inbam)
    b = pybedtools.BedTool(repr_gene)
    a_and_b = a.intersect(b, wo=True, s=True, split=True, bed=True)
    logging.info('[Finish] bedtools intersect')

    gene_info_dict = {}
    gene_info_tuple = namedtuple('gene_info', ['gene_start', 'gene_end', 'gene_exons', 'gene_introns'])
    read_info_dict = {}

    for n, item in enumerate(a_and_b, 1):
        if n % 100000 == 0:
            logging.info(f'process {n} items')
            
        chromosome, readStart, readEnd, read_id, _, strand, _, _, _, readBlockCount, readBlockSizes, readBlockStarts, geneChromosome, geneStart, geneEnd, geneName, _, _, _, _, _, geneBlockCount, geneBlockSizes, geneBlockStarts, cov = item

        # get current gene information
        gene_id = geneName.split('.')[0]
        if gene_id not in gene_info_dict:
            geneBlockSizes = np.fromstring(geneBlockSizes, sep=',', dtype='int')
            geneBlockStarts = np.fromstring(geneBlockStarts, sep=',', dtype='int')
            geneStart = int(geneStart)
            geneEnd = int(geneEnd)

            if strand == '+':
                gene_exons = np.array([(i, j) for i, j in zip(geneStart+geneBlockStarts, geneStart+geneBlockStarts+geneBlockSizes)])
                gene_introns = np.array([(gene_exons[i][1], gene_exons[i+1][0]) for i in range(len(gene_exons)-1)])
            else:
                gene_exons = np.array([(-j, -i) for i, j in zip(geneStart+geneBlockStarts, geneStart+geneBlockStarts+geneBlockSizes)][::-1])
                gene_introns = np.array([(gene_exons[i][1], gene_exons[i+1][0]) for i in range(len(gene_exons)-1)])

            current_gene_info = gene_info_tuple(geneStart, geneEnd, gene_exons, gene_introns)
            gene_info_dict[gene_id] = current_gene_info
        else:
            current_gene_info = gene_info_dict[gene_id]

        
        # get_read_splice_info
        cov = int(cov)
        if read_id not in read_info_dict:
            read_info_dict[read_id] = None  # cov, gene_id, unsplice_count, span_intron, unsplice_intron, relative_pos_to_3ss
        elif read_info_dict[read_id][0] > cov:
            continue

        # skip intronless gene
        if len(current_gene_info.gene_introns) == 0:
            read_info_dict[read_id] = [cov, gene_id, 0, None, None, None]
            continue

        readStart = int(readStart)
        readEnd = int(readEnd)
        readBlockSizes = np.fromstring(readBlockSizes, sep=',', dtype='int')
        readBlockStarts = np.fromstring(readBlockStarts, sep=',', dtype='int')

        if strand == '+':
            read_blocks = [(i, j) for i, j in zip(readStart+readBlockStarts, readStart+readBlockStarts+readBlockSizes)]
            relative_pos_to_3ss = readEnd - np.array(current_gene_info.gene_introns)[:, 1]
        else:
            read_blocks = [(-j, -i) for i, j in zip(readStart+readBlockStarts, readStart+readBlockStarts+readBlockSizes)][::-1]
            relative_pos_to_3ss=-np.array(current_gene_info.gene_introns)[:, 1] - readStart
        
        unsplice_count, span_intron, unsplice_intron = get_splice_stats(read_blocks, current_gene_info.gene_introns)

        read_info_dict[read_id] = [cov, gene_id, unsplice_count, span_intron, unsplice_intron, relative_pos_to_3ss]

    return read_info_dict


def add_tags(infile: str, outfile: str, read_info_dict: dict, min_mapq: int = None):
    '''Add splice info tags to bam

    Args:
        infile: BAM file
        outfile: BAM file
        read_info_dict: generate from `get_read_info_dict`
        min_mapq: Minimum mapping quality score
    '''
    inbam = pysam.AlignmentFile(infile, 'rb')
    outbam = pysam.AlignmentFile(outfile, 'wb', template=inbam)
    
    for read in tqdm(inbam, desc='Add splice info tags to bam'):
        if min_mapq is not None and read.mapping_quality < min_mapq:
            continue
        if read.query_name in read_info_dict:
            cov, gene_id, unsplice_count, span_intron, unsplice_intron, relative_pos_to_3ss = read_info_dict[read.query_name]
            span_intron = ':'.join(map(str, span_intron)) if span_intron is not None else ''
            unsplice_intron = ':'.join(map(str, unsplice_intron)) if unsplice_intron is not None else ''
            relative_pos_to_3ss = ':'.join(map(str, relative_pos_to_3ss)) if relative_pos_to_3ss is not None else ''

            read.set_tag('gi', gene_id)
            read.set_tag('sn', span_intron)
            read.set_tag('rn', unsplice_count)
            read.set_tag('ri', unsplice_intron)
            read.set_tag('rp', relative_pos_to_3ss)

        else:
            read.set_tag('gi', 'None')
            read.set_tag('sn', 'None')
            read.set_tag('rn', 'None')
            read.set_tag('ri', 'None')
            read.set_tag('rp', 'None')
        
        outbam.write(read)
    
    inbam.close()
    outbam.close()

    pysam.index(outfile, '-@ 10')


def main(inbam: str, outbam: str, repr_gene: str, min_mapq: int = None):
    read_info_dict = get_read_info_dict(inbam, repr_gene)
    add_tags(inbam, outbam, read_info_dict, min_mapq=min_mapq)


@click.command()
@click.option('--inbam', required=True)
@click.option('--outbam', required=True)
@click.option('--repr_gene', required=True)
@click.option('--min_mapq', required=False, default=None)
def main_click(inbam, outbam, repr_gene, min_mapq):
    main(inbam, outbam, repr_gene, min_mapq)


if __name__ == '__main__':
    main_click()