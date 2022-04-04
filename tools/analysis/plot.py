#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-21 10:56:29
LastEditTime : 2022-03-21 13:31:37
LastEditors  : windz
FilePath     : plot.py
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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


def plot_boxplot_with_jitter(
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
    figsize=(4, 3)
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

    ax.legend(frameon=False)
    plt.xlabel(xlabel)
    sns.despine(top=True, right=True)

    return ax