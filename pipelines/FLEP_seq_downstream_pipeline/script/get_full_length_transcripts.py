#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2020-04-15 15:26:26
LastEditTime : 2021-10-09 16:43:49
Description  : get full length transcripts from bam
'''


import pysam
import pandas as pd
import click


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('--full_len', required=True)
@click.option('--first_exon_path', required=True)
def main(infile, full_len, first_exon_path):
    #first_exon_path = '/public/home/mowp/test/nanopore_cdna/supplementary_data/representative_gene_model/representative_gene_first_exon.tsv'
    first_exon_df = pd.read_csv(first_exon_path, sep='\t')
    first_exon_df.set_index(['gene_id'], inplace=True)
    first_exon_dict = first_exon_df.to_dict(orient='index')

    #infile = '/public/home/mowp/test/nanopore_cdna/aligned_data/fhh.tagged.mm2.sorted.bam'
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        full_len_bam = pysam.AlignmentFile(full_len, 'wb', template=inbam)
        for read in inbam:
            read_gene_id = read.get_tag('gi')

            if read_gene_id in first_exon_dict:
                # 过滤与基因方向不一致的reads
                if first_exon_dict[read_gene_id]['strand'] == '+' and read.is_reverse:
                    continue
                if first_exon_dict[read_gene_id]['strand'] == '-' and not read.is_reverse:
                    continue

                if (first_exon_dict[read_gene_id]['strand'] == '+' and 
                    read.reference_start <= first_exon_dict[read_gene_id]['exon_end']):
                    full_len_bam.write(read)
                elif (first_exon_dict[read_gene_id]['strand'] == '-' and
                    read.reference_end >= first_exon_dict[read_gene_id]['exon_start']):
                    full_len_bam.write(read)

    
    full_len_bam.close()
    pysam.index(full_len)

    
if __name__ == "__main__":
    main()