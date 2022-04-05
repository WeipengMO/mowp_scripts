#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-21 10:56:29
LastEditTime : 2022-03-23 23:34:33
LastEditors  : windz
FilePath     : plot.py
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def color_pal(key):
    '''
    Custom color palette
    
    To preview: 
        for pal in color_palette.values():
            sns.palplot(pal) 
    '''

    color_palette = {
        'set1': ['#6e9ece', '#e6928f', '#4e9595', '#84574d', '#8d6ab8', '#efdbb9', '#76ba80'], 
        'set2': ['#46a1cd', '#ce3c35', '#4258a1', '#57b058', '#7b4b99', '#f2df52', '#a9a9a9'],
        'default': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    }

    return color_palette[key]


def plot_scatter_plot(
    data1: pd.DataFrame,
    data2: pd.DataFrame,
    on: str = 'gene_id',
    key: str = 'data',
    data1_label: str = None,
    data2_label: str = None,
    data1_font_style: str = None,
    data2_font_style: str = None,
    s: int = 5,
    p_values: float = .001,
    figsize: tuple = (3, 3),
    xlim: tuple = (0, 1000),
    ylim: tuple = (0, 1000)
):
    '''
    Plot scatterplot with p_values
    Args:
        data1, data2: input DataFrame with raw data
        on: Column or index level names to join on
        key: Column names of the raw data
        s: dot size
        p_values: threshold of the significant

    Return:
        data: merge Dataframe
        ax: pre-existing axes for the plot
    '''
    from scipy.stats import mannwhitneyu

    data = pd.merge(data1, data2, on=on)
    data.dropna(inplace=True)

    data['p_values'] = data.apply(lambda x: mannwhitneyu(
        x[f'{key}_x'], x[f'{key}_y'])[1], axis=1)
    data['data_x'] = data[f'{key}_x'].map(np.median)
    data['data_y'] = data[f'{key}_y'].map(np.median)

    data['color'] = data['p_values'].map(lambda x: 1 if x < p_values else 0)

    plt.figure(figsize=figsize)
    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values > @p_values'),
        s=s,
        c='#D2D2D2',
        alpha=.5
    )

    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values < @p_values and data_y > data_x'),
        s=s,
        c='#3182BD',
        alpha=.5
    )

    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values < @p_values and data_y < data_x'),
        s=s,
        c='#C00000',
        alpha=.5
    )

    plt.plot(xlim, ylim, ls='--', color='#555555')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(data1_label, style=data1_font_style)
    plt.ylabel(data2_label, style=data2_font_style)
    plt.xticks(np.linspace(xlim[0], xlim[1], 5, dtype=int))
    plt.yticks(np.linspace(ylim[0], ylim[1], 5, dtype=int))

    ax = plt.gca()
    ax.patch.set_visible(False)

    sns.despine(top=True, right=True)

    return data, ax


def boxplot_with_jitter(
    data: list, 
    figsize: tuple = (3, 3), 
    widths = .5, 
    subsample = 1, 
    labels=None
):
    '''
    带有扰动的boxplot
    Args: 
        data: boxplot data, eg. [[...], [...], ...]

    Return:
        ax: pre-existing axes for the plot
    '''
    from random import sample

    if subsample > 1:
        raise ValueError(f'{subsample=}, must <= 1')

    plt.figure(figsize=figsize)
    plt.boxplot(
        data,
        labels=labels,
        widths=widths,
        showfliers=False
    )

    ylim = plt.ylim()
    for i in range(len(data)):
        y = sample(list(data[i]), int(subsample*len(data[i])))
        x = np.random.normal(1+i, 0.04, size=len(y))
        plt.scatter(x, y, s=1, color='grey', alpha=.5)

    ylim = plt.ylim(ylim)

    sns.despine(top=True, right=True)
    ax = plt.gca()

    return ax


def step_histplot(
    data: list,
    bins = 'auto',
    labels = None,
    stat = 'density',
    xlabel = None,
    figsize=(4, 3),
    bbox_to_anchor: tuple = (1, 0, .5, 1)
):
    if labels is not None:
        labels = iter(labels)
    else:
        labels = iter([None]*len(data))

    plt.figure(figsize=figsize)
    for _data in data:
        ax = sns.histplot(
            _data,
            fill=False,
            stat=stat,
            element='step',
            bins=bins,
            label=next(labels)
        )

    ax.legend(frameon=False, bbox_to_anchor=bbox_to_anchor)
    plt.xlabel(xlabel)
    sns.despine(top=True, right=True)

    return ax


def read_count_per_gene(
    gene_counts: dict,
    bs = .2,
    max_counts = 3,
    all_gene_counts = None,
    ylabel = 'Gene counts',
    ):
    if all_gene_counts is None:
        # 如果没有给出基因总数，则默认注释文件里面的基因总数
        from utils.get_overlap_genes import _, total_gene_counts
        all_gene_counts = total_gene_counts

    df = pd.DataFrame(gene_counts.items(), columns=['gene_id', 'counts'])
    read_count_per_gene = list(np.log10(df['counts']+1))
    gene_with_null_read = all_gene_counts-len(df)  # 在文库中没有被测到的基因总数
    read_count_per_gene.extend([0]*gene_with_null_read)
    read_count_per_gene = np.array(read_count_per_gene)

    counts = read_count_per_gene

    i = 0
    x, y = [], []
    while i <= max_counts:
        mask = (counts <= i) & (counts > i-bs)
        y.append(len(counts[mask]))
        x.append(i)
        i += bs

    plt.figure(figsize=(4, 3))
    plt.bar(x, y, width=bs)
    plt.xlabel('$\log_{10}\mathrm{(read\ counts + 1)}$')
    plt.ylabel(ylabel)
    sns.despine(top=True, right=True)

    ax = plt.gca()

    return ax