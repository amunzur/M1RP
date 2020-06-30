# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:15:00 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import time
import multiprocessing
from joblib import Parallel, delayed
import sys

'''
Calculate adjusted Log_ratio, std, and p_value
'''
cohort = 'Abienza'

# =============================================================================
# Helpers
# =============================================================================

def noisy_tcs(tcs):
    return np.array(tcs.apply(lambda x: x * (2**np.random.normal(0,0.1)) + np.random.normal(0, 0.02)).clip(lower = 0.001, upper = 1))

# =============================================================================
# Import data
# =============================================================================

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Abi-enza_logratios_melted.tsv', sep = '\t')
gene_noise = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/normal_samples_cn/abienza_normal_samples_cn.tsv', sep = '\t')
tumor_fraction = pd.read_csv('C:/Users/amurtha/Dropbox/2017 - Abi-enza second manuscript/Data freeze/ctdna_fractions.tsv', sep = '\t', header = None, names = ['Sample ID','ctDNA%'])

tumor_fraction['ctDNA%'] = tumor_fraction['ctDNA%'] / 100

# =============================================================================
# Merge sample noise and gene noise
# =============================================================================

cn = cn.merge(gene_noise, on = 'GENE', how = 'left')
cn = cn.merge(tumor_fraction, on = 'Sample ID', how = 'left')

# =============================================================================
# Add ploidy 
# =============================================================================

cn['ploidy'] = 2
cn.loc[cn['CHROMOSOME'] == 'chrX', 'ploidy'] = 1

# =============================================================================
# Calculate Copies
# =============================================================================

cn['copies'] = cn['ploidy'] * (1 + (2**(cn['Log_ratio']-cn['lr_mean'])-1)/cn['ctDNA%'])

# =============================================================================
# Calculate z scores, pvalues
# =============================================================================

cn['p_val'] = cn.apply(lambda x: stats.norm.sf(abs(x['Log_ratio']-x['lr_mean']),loc = 0, scale = x['lr_std']), axis = 1)

cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/abi_enza/Targeted_sequencing_error_correction/CANdy_%s_ctDNA_cna.tsv' % cohort, sep = '\t', index = None)

