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
import sys

colors = {
    4:'#EE2D24',
    3:'#F59496',
    2:'#E6E7E8',
    1:'#9CC5E9',
    0:'#3F60AC',
    -1:'#efefef',
    }

if len(sys.argv)>1:
    cohort = sys.argv[1]
else:
    cohort = 'M1RP'

# =============================================================================
# Helpers
# =============================================================================

def visualzie_cn(sample, tmp):
    fig,ax = plt.subplots(figsize =(15,5))
    
    ax.scatter(tmp['GENE'], tmp['Log_ratio'] - tmp['lr_mean'], color = tmp['color'], s = 50, zorder = 100)
    
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

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Targeted_sequencing_error_correction/CANdy_%s_ctDNA_cna.tsv' % cohort, sep = '\t')

# =============================================================================
# Add color column based on p-value and direction of change
# =============================================================================

cn['color'] = ''
for index, row in cn.iterrows():
    cn.at[index, 'color'] = colors.get(row['Adjusted_copy_num'])

# =============================================================================
# Split by sample and plot
# =============================================================================

num_cores = multiprocessing.cpu_count()
Parallel(n_jobs = num_cores)(delayed(visualzie_cn)(sample, cn[cn['Sample ID'] == sample]) for sample in cn['Sample ID'].unique().tolist())