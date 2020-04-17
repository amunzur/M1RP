# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:04:03 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats


# =============================================================================
# Visualize the graph of theoretical cutoff of cancer cell fraction as std of the gene increases
# =============================================================================

min_p = 0.001 #1/1000 chance that the difference is due to noise
ploidy = 2
lr_mean = 0

# =============================================================================
# Set up dataframe
# =============================================================================
ccf_min = pd.DataFrame({'std':np.arange(0.001,0.4, 0.0001)})

# =============================================================================
# Single copy loss
# =============================================================================

copy_change = -1
ccf_min['Single copy loss'] = (2**(-1*stats.norm.ppf(1 - min_p / 2) * ccf_min['std'] + lr_mean) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Single copy gain
# =============================================================================

copy_change = 1
ccf_min['Single copy gain'] = (2**(1*stats.norm.ppf(1 - min_p / 2) * ccf_min['std'] + lr_mean) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Two copy loss
# =============================================================================

copy_change = -2
ccf_min['Two copy loss'] = (2**(-1*stats.norm.ppf(1 - min_p / 2) * ccf_min['std'] + lr_mean) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Two copy gain
# =============================================================================

copy_change = 2
ccf_min['Two copy gain'] = (2**(1*stats.norm.ppf(1 - min_p / 2) * ccf_min['std'] + lr_mean) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Set up graph
# =============================================================================

fig, [ax1,ax2] = plt.subplots(nrows = 2)

ax1.plot(ccf_min['std'],ccf_min['Single copy loss'],  color = 'blue', clip_on = True, label = 'Single copy loss')
ax1.plot(ccf_min['std'],ccf_min['Single copy gain'], color = 'red', clip_on = True, label = 'Single copy gain')
ax1.plot(ccf_min['std'],ccf_min['Two copy loss'], color = 'k', clip_on = True, label = 'Two copy loss')
ax1.plot(ccf_min['std'],ccf_min['Two copy gain'], color = 'green', clip_on = True, label = 'Two copy gain')
ax1.plot(ccf_min['std'], [1]*len(ccf_min), color = 'grey', linestyle = 'dashed')

ax1.set_xlabel('Standard deviation')
ax1.set_ylabel('Tumor fraction')
ax1.set_ylim(0,1)

ax1.legend(fontsize = 8)

# =============================================================================
# Visualize the graph of theoretical cutoff of cancer cell fraction as std of the gene increases
# =============================================================================

min_p = 0.001 #1/100 chance that the difference is due to noise
ploidy = 2
std = 0.2

# =============================================================================
# Set up dataframe
# =============================================================================
ccf_min = pd.DataFrame({'lr_mean':np.arange(-0.4,0.4, 0.0001)})

# =============================================================================
# Single copy loss
# =============================================================================

copy_change = -1
ccf_min['Single copy loss'] = (2**(-1*stats.norm.ppf(1 - min_p) * std + ccf_min['lr_mean']) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Single copy gain
# =============================================================================

copy_change = 1
ccf_min['Single copy gain'] = (2**(1*stats.norm.ppf(1 - min_p) * std + ccf_min['lr_mean']) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Two copy loss
# =============================================================================

copy_change = -2
ccf_min['Two copy loss'] = (2**(-1*stats.norm.ppf(1 - min_p) * std + ccf_min['lr_mean']) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Two copy gain
# =============================================================================

copy_change = 2
ccf_min['Two copy gain'] = (2**(1*stats.norm.ppf(1 - min_p) * std + ccf_min['lr_mean']) - 1)/(((ploidy+copy_change)/ploidy) - 1)

# =============================================================================
# Set up graph
# =============================================================================

ax2.plot(ccf_min['lr_mean'],ccf_min['Single copy loss'],  color = 'blue', clip_on = True, label = 'Single copy loss')
ax2.plot(ccf_min['lr_mean'],ccf_min['Single copy gain'], color = 'red', clip_on = True, label = 'Single copy gain')
ax2.plot(ccf_min['lr_mean'],ccf_min['Two copy loss'], color = 'k', clip_on = True, label = 'Two copy loss')
ax2.plot(ccf_min['lr_mean'],ccf_min['Two copy gain'], color = 'green', clip_on = True, label = 'Two copy gain')
ax2.plot(ccf_min['lr_mean'], [1]*len(ccf_min), color = 'grey', linestyle = 'dashed')

ax2.set_xlabel('Log-ratio mean')
ax2.set_ylabel('Tumor fraction')
ax2.set_ylim(0,1)


fig.tight_layout()

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/std_lrMean_tf_lines (alpha=%f).pdf' % min_p)