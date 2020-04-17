# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:25:17 2020

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

def simulate_lr(copy_change, tc, ploidy):
    cn = ploidy + copy_change
    print(tc)
    adj_tc = max(0, min(0.99999,tc * (2**np.random.normal(0,0.1)) + np.random.normal(0,0.02)))
    sim_lr = math.log2(1 + adj_tc*(cn/ploidy - 1))
    return sim_lr;

def get_str(tc):
    tc = int(tc*100)
    n  = str(tc)
    if tc < 10:
        return '00%s' % n
    elif tc < 100:
        return '0%s' % n
    else:
        return n;

def visualize_lr_dists(tc, lr_dists, p_value):
    fig,ax = plt.subplots(figsize = (5,3))
    ax.hist(lr_dists[-1], bins = 100, color = '#9CC5E9')
    ax.hist(lr_dists[-2], bins = 100, color = '#3F60AC')
    
    ax.text(x = 0, y = 750, s = 'Tumor content = %i\np-value: %f' % (int(tc*100), p_value), horizontalalignment = 'right')
    
    ax.set_ylim(0,800)
    if ax.get_xlim()[0] < -8:
        ax.set_xlim(-8, 0)
    else:
        ax.set_xlim(right = 0)
    
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expected log-ratio')
    
    fig.tight_layout()
    str_tc = get_str(tc)
    fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/tp53_pchange_example/lr_pngs/%s.png' % str_tc)
    return;
    
def simulate_tc_distributions(tc, ploidy):
    lr_dists = pd.DataFrame()
    for copy_change in [-1,-2]:
        print(tc)
        lr_dists[copy_change] = Parallel(n_jobs = n_cores)(delayed(simulate_lr)(copy_change, tc, ploidy) for i in np.arange(0,10000,1))
    stat,p = stats.ttest_ind(lr_dists[-1], lr_dists[-2], equal_var = False)
    print(p)
    p = round(p,3)
    visualize_lr_dists(tc, lr_dists, p)
    return;

# =============================================================================
# Main
# =============================================================================

''' Create a visualiztion of the change in expected log ratio of copy neutral,
heterozygous loss, and biallelic loss as tumor content increases. This example
will use TP53 as the gene of choice. Additionally, this visualiztion is for the
Ghent M1RP cohort. Other cohorts may appear differently
'''

for tc in np.arange(0,1.01, 0.01):
    simulate_tc_distributions(tc, ploidy)