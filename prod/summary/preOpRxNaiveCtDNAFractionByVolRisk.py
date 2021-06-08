# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:46:49 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

mapping = {'High-volume/High-risk':4,
           'High-volume/Low-risk':3,
           'Low-volume/High-risk':2,
           'Low-volume/Low-risk':1}

labels = ['Low-volume\nLow-risk',
          'Low-volume\nHigh-risk',
          'High-volume\nLow-risk',
          'High-volume\nHigh-risk']

# =============================================================================
# Pre-op Rx naive ctDNA fraction by 
# =============================================================================

tc = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/cfDNA timepoints.xlsx')
clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx')
tc = tc[tc['Sample ID'] != 'M1RP_ID8_cfDNA_2017Jan30']

# =============================================================================
# Get patient order
# =============================================================================

clin = clin.sort_values(['latitude_raw'], ascending = [False])

# =============================================================================
# Sort tc 
# =============================================================================

tc = tc.merge(clin[['patient_id','composite_volrisk']], left_on = 'Patient ID', right_on = 'patient_id')

tc = tc[tc['Status'].isin(['Tx naive','Post NA ADT'])]

tc = tc.sort_values('Date collected', ascending = True).drop_duplicates('Patient ID')

# =============================================================================
# Get x coordinate
# =============================================================================

tc['group'] = tc['composite_volrisk'].map(mapping)
tc['x'] = tc['group'].apply(lambda x: x+np.random.uniform(-0.15,0.15))

boxes = [tc[tc['group'] == i+1]['Final tNGS_TC'].tolist() for i in np.arange(4)]

# =============================================================================
# Plot
# =============================================================================

fig,ax = plt.subplots(figsize = (2.75,2.5))

ax.boxplot(boxes, showfliers = False)
ax.scatter(tc['x'],tc['Final tNGS_TC'], c = 'k', alpha = 0.7, s = 10, lw = 0)

ax.set_xticklabels(['%s\n(%i)'% (labels[i],len(boxes[i])) for i in np.arange(4)], fontsize = 6, ha = 'center')
ax.set_xlabel("Composite volume/risk", fontsize = 6)

ax.tick_params(labelsize = 6)
ax.set_ylabel('Baseline ctDNA fraction', fontsize = 6)
ax.set_ylim(top=1)

fig.tight_layout()
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/ctDNA fraction/preOpRxNaiveCtDNAFractionByVolRisk.pdf')
