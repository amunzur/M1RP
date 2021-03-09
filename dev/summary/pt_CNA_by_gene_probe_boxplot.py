# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 13:30:06 2021

@author: amurtha
"""

# =============================================================================
# Create by patient boxplot 
# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import string
import numpy as np
import matplotlib.gridspec as gridspec
import math
import os
import scipy.stats as stats

# =============================================================================
# Constants
# =============================================================================

sample_cat = {'RP':'Primary',
              'PB':'Primary',
              'MLN': 'Met',
              'cfDNA': 'cfDNA'}

# =============================================================================
# Import samples, cn data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_allSamples_cna_curated.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] > 0.40]

tc['Sample Category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits).apply(lambda x: sample_cat.get(x))
tc = tc[tc['Sample Category'] != 'cfDNA']

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Bring in copy number data
# =============================================================================



for y in np.arange(0,43,1):
    pt = 'ID'+str(y+1)
    pt_samples = tc[tc['Patient ID'] == pt]
    pt_samples = pt_samples.set_index('Sample ID')
    pt_samples = pt_samples.sort_values(['Sample Category','Final tNGS_TC'], ascending=False)    
    max_samples = pt_samples['Sample Category'].value_counts().max()
    
    if len(pt_samples['Sample Category'].value_counts()) <= 1:
            continue;
    
    fig = plt.figure()
    gs = gridspec.GridSpec(1,2*3,width_ratios=[pt_samples['Sample Category'].value_counts()['Primary'],pt_samples['Sample Category'].value_counts()['Met']]*3)
    
    for x_0, gene in enumerate(['TP53','RB1','PTEN']):
        boxplots_pri = []
        boxplots_met = []
        for sample in pt_samples.index.tolist():
            if not os.path.exists('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample):
                continue;
            sample_cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample, sep = '\t')
            cn_row = cn.loc[(sample,gene)]
            sample_cn = sample_cn[sample_cn['CHROM'] == cn_row['CHROMOSOME']]
            sample_cn = sample_cn[(sample_cn['START'] >= cn_row['START']-1000) & (sample_cn['END'] <= cn_row['END']+1000)]
            sample_tc = pt_samples.at[sample,'Final tNGS_TC']
            sample_cn['CN_corrected'] = (2**(sample_cn[sample]+1) + 2*sample_tc - 2)/sample_tc
            if pt_samples.at[sample,'Sample Category'] == 'Primary':
                boxplots_pri.append(sample_cn['CN_corrected'])
            else:
                boxplots_met.append(sample_cn['CN_corrected'])
    
        ax_pri = fig.add_subplot(gs[0,x_0*2])
        ax_met = fig.add_subplot(gs[0,x_0*2+1], sharey = ax_pri)
        
        ax_pri.set_ylim(0,4)
        
        ax_pri.boxplot(boxplots_pri, showfliers=False)
        ax_pri.set_xlabel(gene+'\nprimary.', fontsize = 8)
        
        ax_pri.tick_params(bottom = False, labelbottom = False, labelsize = 8)
        ax_met.tick_params(bottom = False, labelbottom = False, labelsize = 8)
        
        ax_met.boxplot(boxplots_met, showfliers=False)
        ax_met.set_xlabel(gene+'\nmet.', fontsize = 8)
        
        if len(boxplots_pri) > 1:
            p_pri = stats.f_oneway(*boxplots_pri)[1]
            ax_pri.text(1, 3.9, 'p=%.4f' % p_pri, fontsize = 8)
        if len(boxplots_met) > 1:
            p_met = stats.f_oneway(*boxplots_met)[1]
            ax_met.text(1, 3.9, 'p=%.4f' % p_met, fontsize = 8)
        
        if x_0 == 0:
            ax_pri.set_ylabel("TC-corrected copy-number", fontsize = 8)
        fig.suptitle(pt)
        
        
    fig.tight_layout()
    fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Probe level CNA analysis/%s_boxplots.pdf' % pt)