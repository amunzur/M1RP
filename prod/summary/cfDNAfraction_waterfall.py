# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:05:43 2021

@author: amurtha
"""

import pandas as pd
import scipy.stats as stats
import numpy as np
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

clin['adt_to_crpc'] = (clin['crpc_date'] - clin['mHSPC_adt_start']) / np.timedelta64(1, 'M')

# =============================================================================
# Sort tc 
# =============================================================================

tc = tc.merge(clin[['patient_id','adt_to_crpc','crpc_status']], left_on = 'Patient ID', right_on = 'patient_id')

tc = tc[tc['Status'].isin(['Tx naive','Post NA ADT'])]

tc = tc.sort_values('Date collected', ascending = True).drop_duplicates('Patient ID')

# =============================================================================
# 
# =============================================================================

tc = tc.sort_values(['crpc_status','adt_to_crpc'], ascending = False)
tc['color'] = 'red'
tc.loc[tc['Final tNGS_TC'] == 0, 'color'] = 'grey'

tc['y'] = np.arange(len(tc))

# =============================================================================
# Fig, ax = 
# =============================================================================

fig,[ax0,ax] = plt.subplots(ncols = 2, figsize = (5.5,3.5), gridspec_kw={'width_ratios':[0.1,1]}, sharey = True)

for index, row in tc.iterrows():
    crpc_time = row['adt_to_crpc']
    y = row['y']
    ax.plot([0,crpc_time],[y+0.4,y+0.4], color = row['color'], lw = 1.5, marker = None, solid_capstyle = 'butt')
    
cspc = tc[tc['crpc_status'] == 0].copy()
crpc = tc[tc['crpc_status'] == 1].copy()

ax.scatter(cspc['adt_to_crpc'],cspc['y']+0.4, c = cspc['color'], marker = '>', lw = 0)
ax.scatter(crpc['adt_to_crpc'],crpc['y']+0.4, c = crpc['color'], marker = '|', lw = 0.8)

ax.set_yticks(np.arange(0.4,len(tc),1))
ax.set_yticklabels(tc['Patient ID'])
ax.set_ylim(-0.1,23.45)

ax.set_xlim(-1,65)
ax.set_xticks(np.arange(0,63,12))

order = dict(zip(tc['Patient ID'],tc['y']))

ax.set_xlabel('Time from ADT to CRPC', fontsize = 6)

handles = [Line2D([],[],lw = 1.5, color = 'red', marker = None),
           Line2D([],[],lw = 1.5, color = 'grey', marker = None),
           Line2D([],[],lw = 0.0, color = 'grey', marker = '|'),
           Line2D([],[],lw = 0.0, color = 'grey', marker = '>')]

labels = ['Baseline ctDNA+','Baseline ctDNA-', 'CRPC progression','Last follow up']

ax.legend(handles, labels, fontsize = 6)

# =============================================================================
# Set up treatment regimine heatmap
# =============================================================================

clin['hspc_arpi'] = 0
for index, row in clin.iterrows():
    if row['crpc_status'] == 1 and row['mCRPC_L1_arpi_administered'] == 1:
        if row['mCRPC_L1_arpi_start'] < row['crpc_date']:
            clin.at[index, 'hspc_arpi'] = 1
    elif row['mCRPC_L1_arpi_administered'] == 1:
        clin.at[index, 'hspc_arpi'] = 1

ti = clin.copy()[['patient_id','neoadj_tx','hspc_arpi','mHSPC_chemo_administered']]
ti['neoadj_tx'] = ti['neoadj_tx'].map({0:'grey',1:'black'})
ti['hspc_arpi'] = ti['hspc_arpi'].map({0:'grey',1:'black'})
ti['mHSPC_chemo_administered'] = ti['mHSPC_chemo_administered'].map({0:'grey',1:'black'})

ti = ti.merge(tc[['Patient ID','y']], left_on = 'patient_id', right_on = 'Patient ID')

ax0.bar([0]*len(ti),0.8, bottom = ti['y'], color = ti['neoadj_tx'])
ax0.bar([1]*len(ti),0.8, bottom = ti['y'], color = ti['mHSPC_chemo_administered'])
ax0.bar([2]*len(ti),0.8, bottom = ti['y'], color = ti['hspc_arpi'])

# ti = pd.DataFrame(columns = ['Neo-adjuvant ADT','mCSPC ARPI','mCSPC chemo.'], index = tc['Patient ID'])

ax0.set_xticks([0,1,2])
ax0.set_xticklabels(['Neo-adj.','mHSPC chemo.','mHSPC ARPI'], fontsize = 8, rotation = 90)
ax0.tick_params(labelsize = 6)
ax.tick_params(labelsize = 6)




# =============================================================================
# Save fig
# =============================================================================

fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/ctDNA fraction/waterfallByCtDNAStatus.pdf')