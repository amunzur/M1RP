# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:56:48 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Import data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/june2021_wes_mutations_betastasis_melted.tsv', sep = '\t')
cn_t = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
cn_w = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_ID8-segments.txt', sep = '\t')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)
tc = tc[tc['Cohort'] == 'M1RP']
tc = tc[tc['Final tNGS_TC'] > 0]

# =============================================================================
# Keep WES samples
# =============================================================================

samples = cn_w.columns.tolist()[3:]
tc = tc[tc['Sample ID'].isin(samples)]
muts = muts[muts['Sample ID'].isin(samples)]
cn_t = cn_t[cn_t['Sample ID'].isin(samples)]

# =============================================================================
# Mut counts
# =============================================================================

muts_count = muts.groupby('Sample ID').count()
muts_count = muts_count[['CHROM']]
muts_count.columns = ['Count']
muts_count = muts_count.reset_index()

# =============================================================================
# Relate to log ratio
# =============================================================================

cn_t = cn_t[cn_t['GENE'].isin(['MSH2','MSH6'])]
cn_t = cn_t.merge(tc[['Sample ID','Final tNGS_TC']], on = 'Sample ID')

cn_t['Est. copies'] = (2**(cn_t['Log_ratio']+1)+2*cn_t['Final tNGS_TC']-2)/cn_t['Final tNGS_TC']

cn_t = cn_t.merge(muts_count, on = 'Sample ID', how = 'left')

# =============================================================================
# Set colors
# =============================================================================

cn_t['color'] = cn_t['GENE'].map({'MSH2':'orange','MSH6':'blue'})

cn_t['Sample cat'] = 'Primary'
cn_t.loc[cn_t['Sample ID'].str.contains('MLN'), 'Sample cat'] = 'Metastatic'

# =============================================================================
# 
# =============================================================================

fig,ax1 = plt.subplots()

cn_t_pri = cn_t[cn_t['Sample cat'] == 'Primary']
cn_t_met = cn_t[cn_t['Sample cat'] == 'Metastatic']

ax1.scatter(cn_t_pri['Count'],cn_t_pri['Est. copies'], color = cn_t_pri['color'], marker = 'o', lw = 0)
ax1.scatter(cn_t_met['Count'],cn_t_met['Est. copies'], color = cn_t_met['color'], marker = '*', lw = 0)

ax1.set_xlabel('Number of mutations')
ax1.set_ylabel('Estimated copies')

from matplotlib.lines import Line2D

handles = [Line2D([],[], color = 'orange', lw = 0, marker = 'o', markeredgewidth=0), Line2D([],[], color = 'blue', lw = 0, marker = 'o', markeredgewidth=0),
           Line2D([],[], color = 'black', lw = 0, marker = 'o', markeredgewidth=0), Line2D([],[], color = 'black', lw = 0, marker = '*', markeredgewidth=0)]
labels = ['MSH2','MSH6', 'Primary','Metastatic']

ax1.legend(handles, labels)