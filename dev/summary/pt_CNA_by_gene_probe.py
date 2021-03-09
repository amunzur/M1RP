# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 11:31:07 2021

@author: amurtha
"""

# =============================================================================
# Copy number comparison between discordard tumor supressors
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

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Bring in copy number data
# =============================================================================

pt_gene = {'ID2': ['TP53','RB1','PTEN'],
           'ID11': ['PTEN','RB1']}

for pt in pt_gene.keys():
    pt_samples = tc[tc['Patient ID'] == pt]
    pt_samples = pt_samples.set_index('Sample ID')
    max_samples = pt_samples['Sample Category'].value_counts().max()
    
    fig = plt.figure(figsize = (2.5*len(pt_gene.get(pt))+1,1.5*max_samples))
    gs = gridspec.GridSpec(max_samples,len(pt_gene.get(pt)*2))
    for x_0, gene in enumerate(pt_gene.get(pt)):
        met_count = 0
        primary_count = 0
        for sample in pt_samples.index.tolist():
            if not os.path.exists('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample):
                continue;
            if pt_samples.at[sample, 'Sample Category'] == 'Primary':
                y = primary_count
                x = x_0*2
                primary_count += 1
            else: 
                y = met_count
                met_count += 1
                x = x_0*2+1
            ax = fig.add_subplot(gs[y,x])
            ax.set_xlabel(sample, fontsize = 8)
            sample_cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample, sep = '\t')
            cn_row = cn.loc[(sample,gene)]
            sample_cn = sample_cn[sample_cn['CHROM'] == cn_row['CHROMOSOME']]
            sample_cn = sample_cn[(sample_cn['START'] >= cn_row['START']-1000) & (sample_cn['END'] <= cn_row['END']+1000)]
            
            ax.scatter(sample_cn['START'], sample_cn[sample], c = 'k', s = 7)
            ax.plot([cn_row['START'] - 1000, cn_row['END'] + 1000], [cn_row['Log_ratio'], cn_row['Log_ratio']], ls = 'dashed', lw = 0.8, color = 'k')
            
            sample_tc = pt_samples.at[sample,'Final tNGS_TC']
            pred_hetDel = math.log2(2-sample_tc)-1
            pred_homoDel = math.log2(1-sample_tc)
            ax.plot([cn_row['START'] - 1000, cn_row['END'] + 1000], [pred_hetDel, pred_hetDel], lw = 0.8, color = '#9CC5E9')
            ax.plot([cn_row['START'] - 1000, cn_row['END'] + 1000], [pred_homoDel, pred_homoDel], lw = 0.8, color = '#3F60AC')
            
            ax.tick_params(bottom = False, labelbottom = False, labelsize = 8)
            
            ax.set_xlim(cn_row['START'] - 1000, cn_row['END'] + 1000)
            ax.set_ylim(-1.5, 0.5)
            ax.set_yticks(np.arange(-1.5,1,0.5))
            
            if y == 0:
                if x%2 == 0:
                    ax.set_title('%s Primary' % gene, fontsize = 8)
                else:
                    ax.set_title('%s Met.' % gene, fontsize = 8)
            
    fig.tight_layout()
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Probe level CNA analysis/%s.pdf' % pt)