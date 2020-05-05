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

'''
Calculate adjusted log_ratio, std, z_score, and p_value
'''

# =============================================================================
# Helpers
# =============================================================================

def noisy_tcs(tcs):
    return np.array(tcs.apply(lambda x: x * (2**np.random.normal(0,0.1)) + np.random.normal(0, 0.02)).clip(lower = 0.001, upper = 1))

'''
def calculate_z_scores(cn, copies_mean, copies_std, zipped_lr_adjusted_std):
    noisy_tc = noisy_tcs(cn['mut_TC'])
    noisy_lr = np.array(list(map(lambda x: np.random.normal(x[0], x[1]), zipped_lr_adjusted_std)))
    sim_copies = cn['ploidy'] * (1 + (2**noisy_lr-1)/noisy_tc)
    print(time.time()-start_time)
    sim_df = pd.DataFrame({0:copies_mean, 1:copies_std, 2:sim_copies})
    z_scores = sim_df.apply(lambda x: (x[2]-x[0]) / x[1], axis = 1)
    print(time.time()-start_time)
    return z_scores
'''

# =============================================================================
# Import data
# =============================================================================

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/gene_cna_melted.tsv', sep = '\t')
sample_noise = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/sample_noise.tsv', sep = '\t')
gene_noise = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/normal_samples_cn/normal_samples_cn.tsv', sep = '\t')
tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction = tumor_fraction.set_index('Sample ID')
tumor_fraction['mut_TC'] = tumor_fraction['mut_TC'].str.split('%').str[0].astype(np.float64) / 100


# =============================================================================
# Merge sample noise and gene noise
# =============================================================================

cn = cn.merge(sample_noise, on = 'Sample ID', how = 'left')
cn = cn.merge(gene_noise, on = 'GENE', how = 'left')
cn = cn.merge(tumor_fraction[['mut_TC']], left_on = 'Sample ID',right_index = True, how = 'left')

# =============================================================================
# Keep relevant samples
# =============================================================================

# cn = cn[cn['mut_TC'] > 0.2]

# =============================================================================
# Calulate adjusted noise
# =============================================================================

# Noise
avg_std = sample_noise['Sample_noise'].mean()
cn['Adjusted_std'] = cn['lr_std'] #* (cn['Sample_noise'] / avg_std)

# =============================================================================
# Add ploidy 
# =============================================================================

cn['ploidy'] = 2
cn.loc[cn['CHROMOSOME'] == 'chrX', 'ploidy'] = 1

# =============================================================================
# Calculate Copies
# =============================================================================

cn['copies'] = cn['ploidy'] * (1 + (2**cn['log_ratio']-1)/cn['mut_TC'])

# =============================================================================
# Calculate z scores, pvalues
# =============================================================================

cn['z_score'] = cn.apply(lambda x: (x['log_ratio'] - x['lr_mean']) / x['Adjusted_std'], axis = 1)
cn['p_val'] = cn['z_score'].apply(lambda z: 1 - stats.norm.cdf(abs(z)))

cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t', index = None)

'''
# =============================================================================
# Simulate copy number 100000 times
# =============================================================================

n_sim = 100000
copies_mean=np.array(cn['copy_mean'])
copies_std=np.array(cn['copy_std'])
zipped_lr_adjusted_std = np.array(list(zip(cn['log_ratio'],cn['Adjusted_std'])))
num_cores = multiprocessing.cpu_count()
print(num_cores)

start_time = time.time()
result = Parallel(n_jobs=num_cores)(delayed(calculate_z_scores)(cn,copies_mean,copies_std,zipped_lr_adjusted_std) for r in range(int(n_sim)))

cn['z_total'] = 0.
for r in result:
    cn['z_total'] = cn['z_total'] + r
print(time.time() - start_time)
cn['z_total'] = cn['z_total']/n_sim

cn['p_val'] = cn['z_total'].apply(lambda z: stats.norm.cdf(z))
'''