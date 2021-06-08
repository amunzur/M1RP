# -*- coding: utf-8 -*-
"""
Created on Tue May 18 09:56:05 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

tc = tc.merge(clin[['patient_id','latitude_raw']], left_on = 'Patient ID', right_on = 'patient_id')

tc = tc[tc['Status'].isin(['Tx naive','Post NA ADT'])]

tc = tc.sort_values('Date collected', ascending = True).drop_duplicates('Patient ID')

# =============================================================================
# Get x coordinate
# =============================================================================

tc['x'] = tc['latitude_raw'].apply(lambda x: x+1+np.random.uniform(-0.15,0.15))

boxes = [tc[tc['latitude_raw'] == i]['Final tNGS_TC'].tolist() for i in np.arange(4)]

# =============================================================================
# Plot
# =============================================================================

fig,ax = plt.subplots(figsize = (2,2.5))

ax.boxplot(boxes, showfliers = False)
ax.scatter(tc['x'],tc['Final tNGS_TC'], c = 'k', alpha = 0.7, s = 10, lw = 0)

ax.set_xticklabels(['%i\n(%i)'% (i,len(boxes[i])) for i in np.arange(4)], fontsize = 6, ha = 'center')
ax.set_xlabel("Num. LATITUDE criteria fulfilled", fontsize = 6)

ax.tick_params(labelsize = 6)
ax.set_ylabel('Baseline ctDNA fraction', fontsize = 6)
ax.set_ylim(top=1)

fig.tight_layout()
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/ctDNA fraction/preOpRxNaiveCtDNAFractionByLatCriteria.pdf')
