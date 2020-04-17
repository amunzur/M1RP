 # -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:23:02 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed


alpha = 0.001
cn_neutral_color = '#ababab'
biallelic_deletion_color = '#3F60AC'
amplificaiton_color = '#EE2D24'
deletion_color = '#9CC5E9'


# =============================================================================
# Helpers
# =============================================================================

def visualzie_cn(sample, tmp):
    fig,ax = plt.subplots(figsize =(15,5))
    
    ax.scatter(tmp['GENE'], tmp['log_ratio'] - tmp['lr_mean'], color = tmp['color'], s = 50, zorder = 100)
    
    for index, row in tmp.iterrows():
        ax.plot([row['GENE'],row['GENE']], (-2,2), color = row['color'], linestyle = 'dashed', alpha = 0.6, zorder = 0, linewidth = 0.5)
    
    ax.tick_params(bottom = False)
    
    ax.set_xticklabels(tmp['GENE'], rotation = 90, fontsize = 10)
    
    ax.set_xlim(-0.4, len(tmp)-0.4)
    ax.set_ylim(-2,2)
    
    ax.grid(axis = 'y', linestyle = 'dashed', linewidth = 0.3, zorder = 0, color = '0.5')
    fig.tight_layout()
    plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Figures/cn_noise_based_calls/%s.pdf' % sample)
    return;
    
# =============================================================================
# Import data
# =============================================================================

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t')

# =============================================================================
# Add color column based on p-value and direction of change
# =============================================================================


cn['color'] = cn_neutral_color
cn.loc[(cn['p_val'] <= alpha)&(cn['log_ratio'] >= cn['lr_mean']), 'color'] = amplificaiton_color
cn.loc[(cn['p_val'] <= alpha)&(cn['log_ratio'] <= cn['lr_mean']), 'color'] = deletion_color
cn.loc[cn['mut_TC'] <= cn['min_tc_1loss'], 'color'] = '#efefef'

# =============================================================================
# Split by sample and plot
# =============================================================================

num_cores = multiprocessing.cpu_count()
Parallel(n_jobs = num_cores)(delayed(visualzie_cn)(sample, cn[cn['Sample ID'] == sample]) for sample in cn['Sample ID'].unique().tolist())