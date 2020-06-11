# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:08:53 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np


# =============================================================================
# Import data
# =============================================================================

sample_noise = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/sample_noise.tsv', sep = '\t')

tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction = tumor_fraction.set_index('Sample ID')
tumor_fraction['mut_TC'] = tumor_fraction['mut_TC'].str.split('%').str[0].astype(np.float64) / 100.

# =============================================================================
# Merge sample tumor fraction onto sample noise
# =============================================================================

sample_noise = sample_noise.merge(tumor_fraction[['mut_TC']], left_on = 'Sample', right_index = True, how = 'left')

# =============================================================================
# Keep only FFPE samples
# =============================================================================

sample_noise['Sample type'] = sample_noise['Sample'].str.split('_').str[2]
sample_noise['Sample type'] = sample_noise['Sample type'].str.extract('(\D+)')
sample_noise = sample_noise[sample_noise['Sample type'].isin(['MB','PB','MLN','RP'])]

# =============================================================================
# Create df with no TC negative samples
# =============================================================================

sample_noise = sample_noise[~sample_noise['mut_TC'].isnull()]
sample_noise_TCpos = sample_noise[sample_noise['mut_TC'] > 0]

# =============================================================================
# Calculate r values
# =============================================================================

r_all, p_all = stats.pearsonr(sample_noise['mut_TC'], sample_noise['Sample_noise'])
r_pos, p_pos = stats.pearsonr(sample_noise_TCpos['mut_TC'], sample_noise_TCpos['Sample_noise'])

# =============================================================================
# Create scatter plot
# =============================================================================

fig,[ax1, ax2] = plt.subplots(nrows = 2, sharex = True)

ax1.scatter(sample_noise['mut_TC'], sample_noise['Sample_noise'], s = 8, color = 'k')
ax2.scatter(sample_noise_TCpos['mut_TC'], sample_noise_TCpos['Sample_noise'], s = 8, color = 'k')

ax1.text(x = 1.0, y = 0.15, s = 'r = %s, p = %s' % (str(round(r_all, 2)), str(round(p_all, 4))), horizontalalignment = 'right')
ax2.text(x = 1.0, y = 0.15, s = 'r = %s, p = %s' % (str(round(r_pos, 2)), str(round(p_pos, 4))), horizontalalignment = 'right')

ax1.set_ylim(0.14,0.35)
ax2.set_ylim(0.14,0.35)

ax2.set_xlabel('Tumor content')
ax1.set_ylabel('Sample noise')
ax2.set_ylabel('Sample noise')

fig.tight_layout()

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Figures/tc_sampleNoise_scatter/tc_sampleNoise_scatter.pdf')