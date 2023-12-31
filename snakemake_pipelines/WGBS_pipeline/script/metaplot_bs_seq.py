#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-08-26 10:29:53
LastEditTime : 2022-04-07 16:14:16
LastEditors  : windz
FilePath     : metaplot_bs_seq.py
'''


import numpy as np
import pyBigWig
import pyranges as pr
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import scipy


def get_bin_cov(methyratio: list, bins: int, threshold=4):
    methyratio = np.array(methyratio)
    results = []
    for i in range(0, len(methyratio), bins):
        bin_methyratio = methyratio[i:i + bins]
        if len(bin_methyratio[~np.isnan(bin_methyratio)]) >= threshold:
            # we selected for bins with at least 4 cytosines
            # that are each covered by at least 4 reads
            mean_methyratio = np.nanmean(bin_methyratio)
        else:
            mean_methyratio = np.nan
        results.append(mean_methyratio)

    return results


# site-point
def get_site_cov(
    infile: str, 
    chrom: str, site: int, strand: str,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    chrom_prefix: str = '',
    exclude_chr = None
    ):
    '''
    Args:
        infile: path to bigWig file
        chrom: chromosome name
        site: target site
        strand: '+' or '-'
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''

    if chrom in exclude_chr:
        return
    chrom = chrom_prefix + chrom
        
    bwfile = pyBigWig.open(infile)
    if strand == '+':
        start = site - before
        end = site + after
    else:
        start = site - after
        end = site + before
    if start < 0 or end > bwfile.chroms()[chrom]:
        # remove Invalid interval
        return

    if strand == '+':
        values = bwfile.values(
            chrom, 
            start,
            end)
        cov = get_bin_cov(values, bins)

    elif strand == '-':
        values = bwfile.values(
            chrom, 
            start,
            end)[::-1]
        cov = get_bin_cov(values, bins)

    return cov


def get_meta_site_result(
    infile: str, 
    site_info: list,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    chrom_prefix: str = '',
    exclude_chr = None,
    threads=64):
    '''
    Args:
        infile: path to bigWig file
        site_info: [(chrom, site1, site2, strand), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    chrom = site_info[:, 0]
    site = site_info[:, 1]
    strand = site_info[:, 2]
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info) / threads)
        results = e.map(
            get_site_cov,
            repeat(infile),
            chrom,
            site,
            strand,
            repeat(before),
            repeat(after),
            repeat(bins),
            repeat(chrom_prefix),
            repeat(exclude_chr),
            chunksize=chunksize)

    cov = []
    for res in results:
        if res is not None:
            cov_ = res
            cov.append(cov_)
    
    cov = np.nanmean(cov, axis=0)
    return cov


def set_ax(ax, b1, a1, b2, a2, bins, ylabel=None):
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_ylabel(ylabel)

    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].yaxis.set_ticks_position('none')

    ax[0].set_xticks([0, b1//bins, (a1+b1)//bins])
    ax[0].set_xticklabels([f'-{b1//1000} kb', 'TSS', f'{a1//1000} kb'], rotation=90)

    ax[1].set_xticks([0, b2//bins, (a2+b2)//bins])
    ax[1].set_xticklabels([f'-{b2//1000} kb', 'PAS', f'{a2//1000} kb'], rotation=90)

    ax[0].axvline(b1//bins, ls='--', color='#555555')
    ax[1].axvline(b2//bins, ls='--', color='#555555')


# reference scale
def get_scale_cov(
    infile: str, 
    chrom: str, site1: int, site2: int, strand: str,
    before: int = 1000, after : int = 1000, regionbody : int = 1000, 
    bins: int = 100,
    chrom_prefix: str = '',
    exclude_chr = None
    ):
    '''
    Args:
        infile: path to bigWig file
        chrom: chromosome name
        site1: 5' site
        site2: 3' site
        strand: '+' or '-'
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    site1 = int(site1)
    site2 = int(site2)
    chrom = str(chrom)
    if site2 - site1 < bins:
        return

    if exclude_chr is not None and chrom in exclude_chr:
        return

    chrom = chrom_prefix + chrom

    bwfile = pyBigWig.open(infile)
    if strand == '+':
        start = site1 - before
        end = site2 + after
    else:
        start = site1 - after
        end = site2 + before
    if start < 0 or end > bwfile.chroms()[chrom]:
        # remove Invalid interval
        return
    
    if strand == '+':
        cov_5 = bwfile.values(chrom, start, site1)
        cov_5 = get_bin_cov(cov_5, bins)
        cov_3 = bwfile.values(chrom, site2, end)
        cov_3 = get_bin_cov(cov_3, bins)
        # gene_body_region
        cov_gb = bwfile.values(chrom, site1, site2)
        cov_gb = scipy.ndimage.zoom(
            cov_gb,
            regionbody / len(cov_gb),
            order=0,
            mode='nearest')
        cov_gb = get_bin_cov(cov_gb, bins)

    elif strand == '-':
        cov_5 = bwfile.values(chrom, site2, end)[::-1]
        cov_5 = get_bin_cov(cov_5, bins)
        cov_3 = bwfile.values(chrom, start, site1)[::-1]
        cov_3 = get_bin_cov(cov_3, bins)
        # gene_body_region
        cov_gb = bwfile.values(chrom, site1, site2)[::-1]
        cov_gb = scipy.ndimage.zoom(
            cov_gb,
            regionbody / len(cov_gb),
            order=0,
            mode='nearest')
        cov_gb = get_bin_cov(cov_gb, bins)
    else:
        raise ValueError('strand must be "-" or "+"')
    
    cov = np.concatenate([cov_5, cov_gb, cov_3])

    return cov


def get_meta_scale_result(
    infile: str, 
    site_info: list,
    before: int = 1000, after : int = 1000, regionbody : int = 1000, 
    bins: int = 100,
    chrom_prefix: str = '',
    exclude_chr = None,
    threads=64):
    '''
    Args:
        infile: path to bigWig file
        site_info: [(chrom, site1, site2, strand), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    site_info = np.array(site_info)
    chrom = site_info[:, 0]
    site1 = site_info[:, 1]
    site2 = site_info[:, 2]
    strand = site_info[:, 3]
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info) / threads)
        results = e.map(
            get_scale_cov,
            repeat(infile),
            chrom,
            site1,
            site2,
            strand,
            repeat(before),
            repeat(after),
            repeat(regionbody),
            repeat(bins),
            repeat(chrom_prefix),
            repeat(exclude_chr),
            chunksize=chunksize)

    cov = []
    n = 0
    for res in results:
        if res is not None:
            cov_ = res
            cov.append(cov_)
            n += 1
    print(f'number of sites: {n}')
    cov = np.nanmean(cov, axis=0)
    return cov

