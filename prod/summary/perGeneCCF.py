# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 13:45:46 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

sCat_dict = {'cfDNA':'cfDNA', 'MB':'Metastatic', 'PB':'Primary', 'RP':'Primary', 'MLN':'Metastatic', 'MT':'Metastatic', 'UCC':'UCC'}

# =============================================================================
# Import data
# =============================================================================


ccf = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/ccf_estimations_for_targeted_mutations.tsv', sep = '\t')

ccf['Sample cat'] = ccf['Sample ID'].str.split('_').str[2].str.strip(string.digits)
ccf['Sample cat'] = ccf['Sample cat'].map(sCat_dict)

ccf = ccf[ccf['Patient ID'] != 'ID8']
ccf = keepCodingMutations(ccf)

# =============================================================================
# get boxplots
# =============================================================================

pos = []
x = 0

boxes = []
genes = []

for gene in ccf['GENE'].unique().tolist():
    for cat in ['Primary','Metastatic']:
        pos.append(x)
        tmp = ccf[(ccf['Sample cat'] == cat) & (ccf['GENE'] == gene)]
        boxes.append(tmp['ccf_adjusted'])
        x = x + 1
        genes.append('_'.join([gene,cat[0]]))
    x = x+1
        
        
# =============================================================================
# Figure
# =============================================================================
        
fig,ax = plt.subplots(figsize = (10,4))

ax.boxplot(boxes, positions = pos, widths = 0.3, labels = genes)

ax.set_xticklabels(genes, rotation = 90, fontsize = 6)

plt.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/perGeneCCF.pdf')