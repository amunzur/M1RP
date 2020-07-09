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

if len(sys.argv) > 1:
    cohort = sys.argv[1]
    abridged = bool(sys.argv[3])

else:
    cohort = 'M1RP'
    abridged = False

# =============================================================================
# Helpers
# =============================================================================

def noisy_tcs(tcs):
    return np.array(tcs.apply(lambda x: x * (2**np.random.normal(0,0.1)) + np.random.normal(0, 0.02)).clip(lower = 0.001, upper = 1))

# =============================================================================
# Import data
# =============================================================================

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_FFPE_cna.tsv' % cohort, sep = '\t')
gene_noise = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/normal_samples_cn/normal_samples_cn.tsv', sep = '\t')
tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction = tumor_fraction.set_index('Sample ID')
tumor_fraction['Final tNGS_TC'] = tumor_fraction['Final tNGS_TC'].str.split('%').str[0].astype(np.float64) / 100

# =============================================================================
# Merge sample noise and gene noise
# =============================================================================

cn = cn.merge(gene_noise, on = 'GENE', how = 'left')
cn = cn.merge(tumor_fraction[['Final tNGS_TC']], left_on = 'Sample ID',right_index = True, how = 'left')

# =============================================================================
# Add ploidy 
# =============================================================================

cn['ploidy'] = 2
cn.loc[cn['CHROMOSOME'] == 'chrX', 'ploidy'] = 1

# =============================================================================
# Calculate Copies
# =============================================================================

cn['copies'] = cn['ploidy'] * (1 + (2**(cn['Log_ratio']-cn['lr_mean'])-1)/cn['Final tNGS_TC']).clip(lower = 0)

# =============================================================================
# Calculate z scores, pvalues
# =============================================================================

cn['p_val'] = cn.apply(lambda x: stats.norm.sf(abs(x['Log_ratio']-x['lr_mean']),loc = 0, scale = x['lr_std']), axis = 1)

cn['Patient ID'] = cn['Sample ID'].str.split('_').str[1]

if abridged == True:
    cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_%s_FFPE_cna_abridged.tsv' % cohort, sep = '\t', index = None)
else:
    cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_%s_FFPE_cna.tsv' % cohort, sep = '\t', index = None)
