# -*- coding: utf-8 -*-
"""
Created on Tue May 18 14:23:55 2021

@author: amurtha
"""

import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =============================================================================
# TC by num bone mets
# =============================================================================

tc = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/cfDNA timepoints.xlsx')
clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx')
tc = tc[tc['Sample ID'] != 'M1RP_ID8_cfDNA_2017Jan30']

# =============================================================================
# Get patient order
# =============================================================================

clin = clin.sort_values(['latitude_raw'], ascending = [False])

clin.loc[~clin['soft_tissue_site'].isna(), 'visceral'] = 1

# =============================================================================
# Sort tc 
# =============================================================================

tc = tc.merge(clin[['patient_id','bone_mets_num','visceral','soft_tissue_site']], left_on = 'Patient ID', right_on = 'patient_id')

tc = tc[tc['Status'].isin(['Tx naive','Post NA ADT'])]

tc = tc.sort_values('Date collected', ascending = True).drop_duplicates('Patient ID')

tc['visceral'] = tc['visceral'].map({1:'blue',0:'grey'})

# =============================================================================
# 
# =============================================================================

fig,ax = plt.subplots(figsize = (2.5,2.5))

ax.scatter(tc['bone_mets_num'],tc['Final tNGS_TC'], c = tc['visceral'], alpha = 0.6, lw = 0, s = 15)

ax.set_xlabel('Num. bone mets', fontsize = 6)
ax.set_ylabel('Baseline ctDNA fraction', fontsize = 6)
ax.tick_params(labelsize = 6)

handles = [Line2D([],[],markerfacecolor='blue',lw=0,markeredgewidth=0,markersize = 5, marker = 'o'),
           Line2D([],[],markerfacecolor='grey',lw=0,markeredgewidth=0,markersize = 5, marker = 'o')]
labels = ['Visceral mets','No visceral mets']

ax.legend(handles,labels, fontsize = 6)

fig.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/ctDNA fraction/preOpRxNaiveCtDNAFractionByNumBoneMets.pdf')

