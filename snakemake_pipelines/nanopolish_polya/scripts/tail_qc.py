#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import click

@click.command()
@click.option("--infile", "-i", help='input file')
@click.option("--outprefix", help='prefix of output file')
def main(infile, outprefix):
    df = pd.read_csv(infile, sep="\t")

    plt.figure(figsize=(3, 3))
    plt.pie(
        df.qc_tag.value_counts().sort_index(),
        labels=df.qc_tag.value_counts().sort_index().index,
        autopct='%.2f%%')
    plt.savefig(outprefix+'.tail_qc.png')

    plt.figure(figsize=(3, 3))
    df.query('qc_tag == "PASS"')['polya_length'].hist(bins=100, range=(0, 1000))
    plt.xlabel('polya_length')
    plt.ylabel('count')
    plt.savefig(outprefix+'.polya_len.png')


if __name__ == '__main__':
    main()