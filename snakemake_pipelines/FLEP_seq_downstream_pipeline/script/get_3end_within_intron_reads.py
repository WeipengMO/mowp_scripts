#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-01-04 19:49:41
LastEditTime : 2022-01-17 11:09:46
LastEditors  : windz
FilePath     : /data/Zhaijx/mowp/data/nanopore_data/Repository/script/get_3end_within_intron_reads.py
'''

import sys


try:
    infile = sys.argv[1]  # 20211028_col_seedling_total_R1.read.info.txt
    outfile = sys.argv[2]  # 20211028_col_seedling_total_R1.intron.read.info.txt
except IndexError:
    print('Input or Output not detected')

if infile.endswith('.gz'):
    import gzip
    f = gzip.open(infile, 'rt')
else:
    f = open(infile, 'r')
    
with open(outfile, 'w') as o:
    line = next(f)
    o.write(line)

    for line in f:
        line_ = line.rstrip().split('\t')
        if line_[24] == 'intron':
            o.write(line)

f.close()