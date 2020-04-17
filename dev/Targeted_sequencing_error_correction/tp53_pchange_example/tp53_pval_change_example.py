# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:44:03 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
import multiprocessing
from joblib import Parallel, delayed

# =============================================================================
# Constants
# =============================================================================

gene = 'TP53'
ploidy = 2
n_cores = multiprocessing.cpu_count()

# =============================================================================
# Helpers
# =============================================================================

def simulate_lr(copy_change, tc, ploidy, xbar, std):
    cn = ploidy + copy_change
    print(tc)
    adj_tc = max(0, min(0.99999,tc * (2**np.random.normal(0,0.1)) + np.random.normal(0,0.02)))
    sim_lr = math.log2(1 + adj_tc*(cn/ploidy - 1))
    pval = 2*(1-stats.norm.cdf(abs((sim_lr - xbar)/std)))
    return pval;

def get_str(tc):
    tc = int(tc*100)
    n  = str(tc)
    if tc < 10:
        return '00%s' % n
    elif tc < 100:
        return '0%s' % n
    else:
        return n;

def visualize_pval_dists(tc, lr_dists, p_value):
    fig,ax = plt.subplots(figsize = (5,3))
    
    hist, bins = np.histogram(lr_dists[-1], bins=100)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    ax.hist(lr_dists[-1], bins = logbins, color = '#9CC5E9')
    
    hist, bins = np.histogram(lr_dists[-2], bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),100)
    ax.hist(lr_dists[-2], bins = logbins, color = '#3F60AC')
    
    ax.text(x = 0, y = 750, s = 'Tumor content = %i\np-value: %f' % (int(tc*100), p_value), horizontalalignment = 'left')
    
    ax.set_ylim(0,800)
        
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expected pvalue')
    
    ax.set_xscale('log')
    
    fig.tight_layout()
    str_tc = get_str(tc)
    fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/tp53_pchange_example/p_value_pngs/%s.png' % str_tc)
    return;
    
def simulate_tc_distributions(tc, ploidy, xbar, std):
    pval_dists = pd.DataFrame()
    for copy_change in [-1,-2]:
        print(tc)
        pval_dists[copy_change] = Parallel(n_jobs = n_cores)(delayed(simulate_lr)(copy_change, tc, ploidy, xbar, std) for i in np.arange(0,10000,1))
    stat,p = stats.ttest_ind(pval_dists[-1], pval_dists[-2], equal_var = False)
    print(p)
    p = round(p,3)
    visualize_pval_dists(tc, pval_dists, p)
    return;

def get_gene_stats(gene):
    tmp = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t')
    tmp = tmp[tmp['GENE'] == gene]
    tmp = tmp.drop_duplicates('GENE').reset_index()
    xbar = tmp['lr_mean'].at[0]
    std = tmp['lr_std'].at[0]
    return xbar, std

# =============================================================================
# Main
# =============================================================================

''' Create a visualiztion of the change in expected log ratio of copy neutral,
heterozygous loss, and biallelic loss as tumor content increases. This example
will use TP53 as the gene of choice. Additionally, this visualiztion is for the
Ghent M1RP cohort. Other cohorts may appear differently
'''

xbar, std = get_gene_stats(gene)

for tc in np.arange(0,1.01, 0.01):
    simulate_tc_distributions(tc, ploidy, xbar, std)