# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 13:06:49 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt

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
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

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
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] >= 0.40]

tc['Sample Category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits).apply(lambda x: sample_cat.get(x))
tc = tc[tc['Sample Category'] != 'cfDNA']

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Bring in copy number data
# =============================================================================



for y in np.arange(0,43,1):
# for y in [1]:
    pt = 'ID'+str(y+1)
    pt_samples = tc[tc['Patient ID'] == pt]
    pt_samples = pt_samples.set_index('Sample ID')
    pt_samples = pt_samples.sort_values(['Sample Category','Final tNGS_TC'], ascending=False)    
    max_samples = pt_samples['Sample Category'].value_counts().max()
    
    if len(pt_samples['Sample Category'].value_counts()) <= 1:
            continue;
    
    fig = plt.figure(figsize = (6.5,2.5))
    gs = gridspec.GridSpec(1,3)
    
    for x_0, gene in enumerate(['TP53','RB1','PTEN']):
        no_copy_change = []
        hetz_del = []
        cats = []
        samples = []
        for sample, cat in zip(pt_samples.index.tolist(), pt_samples['Sample Category'].tolist()):
            if not os.path.exists('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample):
                continue;
            sample_cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample, sep = '\t')
            cn_row = cn.loc[(sample,gene)]
            sample_cn = sample_cn[sample_cn['CHROM'] == cn_row['CHROMOSOME']]
            sample_cn = sample_cn[(sample_cn['START'] >= cn_row['START']-1000) & (sample_cn['END'] <= cn_row['END']+1000)]
            sample_tc = pt_samples.at[sample,'Final tNGS_TC']
            pred_hetDel = math.log2(2-sample_tc)-1
            
            sample_cn['No change'] = sample_cn[sample]**2
            no_copy_change.append(sample_cn['No change'].sum() / len(sample_cn))
            
            sample_cn['Hetz del'] = (sample_cn[sample] - pred_hetDel)**2
            hetz_del.append(sample_cn['Hetz del'].sum() / len(sample_cn))
            
            samples.append(sample)
            cats.append(cat)
            
        df = pd.DataFrame({'Sample ID':samples, 'Diploid': no_copy_change, 'hetz loss': hetz_del, 'Sample Category':cats})
        df['color'] = 'red'
        df.loc[df['Sample Category'] == 'Primary', 'color'] = 'blue'
        ax = fig.add_subplot(gs[0,x_0])
        ax.scatter(df['Diploid'],df['hetz loss'], c = df['color'], s = 8, zorder = 10)
        ax.set_title(gene, fontsize = 8)
        ax.set_xlabel("Var. from CN2 (LR = 0)")
        ax.set_ylabel("Var. from pred. hetz. loss")
        ax.plot([0,1],[0,1], ls = 'dashed', lw = 0.5, c = 'k', zorder = 0)
        
        p1 = Polygon([(0,0),(1,0),(1,1)], color = '#9CC5E9', lw = 0)
        p2 = Polygon([(0,0),(0,1),(1,1)], color = '#EFEFEF', lw = 0)
        
        ax.add_patch(p1)
        ax.add_patch(p2)
        
        if x_0 != 2:
            ax.set_xlim(0, 0.5)
            ax.set_ylim(0, 0.5)
        else:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            
        if x_0 == 0:
            ax.legend([Line2D([0],[0],lw = 0, marker = 'o', markerfacecolor='red', markeredgewidth=0),Line2D([0],[0],lw = 0, marker = 'o', markerfacecolor='blue', markeredgewidth=0)], ['Metastatic','Primary'], fontsize=6)
                        
    fig.tight_layout()
    fig.suptitle(pt)
    fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Probe level CNA analysis/%s_TSloss_variancePrediction.pdf' % pt)