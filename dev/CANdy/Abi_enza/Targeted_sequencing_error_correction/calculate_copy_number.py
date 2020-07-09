# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:25:10 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import multiprocessing
from joblib import Parallel, delayed
import math
import sys

# =============================================================================
# Constants
# =============================================================================
if len(sys.argv) == 1:
    cohort = 'Abienza'
    alpha = 0.001
else:
    cohort = sys.argv[1]
    alpha = float(sys.argv[2])
n_cores = multiprocessing.cpu_count()
n_sims = 1000
amp_ddel_alpha = 0.01

# =============================================================================
# Helpers
# =============================================================================

def get_noise_tc(tc):
    return max(0, min(0.9999,tc * (2**np.random.normal(0,0.1)) + np.random.normal(0,0.02)))
    
def get_lr_dist(row, cn_change):
    adj_tc = get_noise_tc(row['ctDNA%'])
    ploidy = row['ploidy']
    cn = max(0, ploidy + cn_change)
    sim_lr = math.log2(1 + adj_tc*(cn/ploidy - 1)) + row['lr_mean']
    return sim_lr;

def get_copy_number(row):
    direction = -1 if row['Log_ratio'] < 0 else 1
    lr_dist = np.array(Parallel(n_jobs = n_cores)(delayed(get_lr_dist)(row, direction) for i in range(n_sims)))
    if direction == 1: 
        pval = stats.norm.cdf(row['Log_ratio'], loc = lr_dist.mean(), scale = np.std(lr_dist))
    else:
        pval = stats.norm.sf(row['Log_ratio'], loc = lr_dist.mean(), scale = np.std(lr_dist))
    if direction == 1:
        if pval < amp_ddel_alpha and row['Log_ratio'] > lr_dist.mean():
            lr_dist2 = np.array(Parallel(n_jobs = n_cores)(delayed(get_lr_dist)(row, direction*2) for i in range(n_sims)))
            pval2 = stats.norm.cdf(row['Log_ratio'], loc = lr_dist2.mean(), scale = np.std(lr_dist2))
            if pval2 > amp_ddel_alpha or row['Log_ratio'] > lr_dist2.mean() or row['Log_ratio'] - row['lr_mean'] > -1.1:
                return 4;
            else:
                return 3
        elif row['Log_ratio'] - row['lr_mean'] > 1.1:
            return 4;
        else:
            return 3
    elif direction == -1:
        if pval < amp_ddel_alpha and row['Log_ratio'] < lr_dist.mean():
            lr_dist2 = np.array(Parallel(n_jobs = n_cores)(delayed(get_lr_dist)(row, direction*2) for i in range(n_sims)))
            pval2 = stats.norm.sf(row['Log_ratio'], loc = lr_dist2.mean(), scale = np.std(lr_dist2))
            if pval2 > amp_ddel_alpha or row['Log_ratio'] < lr_dist2.mean() or row['Log_ratio'] - row['lr_mean'] < -1.1:
                return 0;
            else:
                return 1
        elif row['Log_ratio'] - row['lr_mean'] < -1.1:
            return 0;
        else:
            return 1
    

# =============================================================================
# Main
# =============================================================================

'''
Add cn_call column to copy number table. For each row where p_value < alpha, 
create a distribution of the expected p-values for 1 copy change. If p_value 
lies below that distribution, consider the change 2 (or more in case of gains). 
Those changes will be given values of -2 or 4. Other significant p values will
be given the copy number value 1, 3. Copy neutral will be given copy number 2. 
If the tumor content is < the min_tc required to reach the current alpha value, 
copy number is assigned -1.
'''

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Targeted_sequencing_error_correction/CANdy_%s_ctDNA_cna.tsv' % cohort, sep = '\t')

cn['Adjusted_copy_num'] = np.nan

for index, row in cn.iterrows():
    if row['p_val'] < alpha / 73:
        cn.at[index, 'Adjusted_copy_num'] = get_copy_number(row)
    elif row['ctDNA%'] < row['min_tc_1loss']:
        cn.at[index, 'Adjusted_copy_num'] = -1
    else:
        cn.at[index, 'Adjusted_copy_num'] = 2
        
cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Targeted_sequencing_error_correction/CANdy_%s_ctDNA_cna.tsv' % cohort, sep = '\t', index = None)