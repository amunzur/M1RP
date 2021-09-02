# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 09:35:07 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# 
# =============================================================================

altered = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_percentAltered.tsv')

altered['Patient ID'] = altered['Sample ID'].str.split('_').str[1]
altered = altered[~altered['Patient ID'].isin(['ID20','ID9'])]

# =============================================================================
# Import clinical data and TP53 mutation status
# =============================================================================

pts = altered['Patient ID'].unique().tolist()

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')

muts = keepCodingMutations(muts)
muts = muts[muts['GENE'] == 'TP53']
muts = muts[muts['Patient ID'].isin(list(pts))]

cn = cn[cn['GENE'] == 'TP53']
cn = cn[cn['Patient ID'].isin(list(pts))]
cn = cn[cn['Copy_num'] == -2]

clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx')
clin = clin[clin['latitude'] == 'High-risk']

# Set up the heatmap df

hm = pd.DataFrame({'Patient':list(pts)})
hm['TP53 status'] = 0
hm.loc[(hm['Patient'].isin(muts['Patient ID'].tolist()))|(hm['Patient'].isin(cn['Patient ID'].tolist())), 'TP53 status'] = 1

hm['High-risk'] = 0
hm.loc[hm['Patient'].isin(clin['patient_id'].tolist()), 'High-risk'] = 1

hm['Secondary tumor'] = 1
hm.loc[hm['Patient'].isin(['ID8','ID30','ID19']), 'Secondary tumor'] = 0

# =============================================================================
# Get order
# =============================================================================

pt_median = altered.groupby('Patient ID').median().sort_values('Percent altered', ascending = False)[['Percent altered']]
hm = hm.merge(pt_median, left_on = 'Patient', right_index = True)
hm = hm.sort_values(['Secondary tumor','TP53 status','Percent altered'], ascending = False)

x_dict = dict(zip(hm['Patient'].tolist(), np.arange(len(hm))))

altered['x'] = altered['Patient ID'].map(x_dict)
hm['x'] = hm['Patient'].map(x_dict)

color_dict = {0:'grey', 1:'black'}
hm['TP53 status'] = hm['TP53 status'].map(color_dict)
hm['High-risk'] = hm['High-risk'].map(color_dict)

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC','M1RP_ID30_RP3']

altered = altered[~altered['Patient ID'].isin(['ID9','ID20'])]

altered['color'] = 'k'
altered.loc[altered['Patient ID'] == 'ID8', 'color'] = 'orange'
altered.loc[altered['Sample ID'].isin(secondary_tumors), 'color'] = 'blue'

# =============================================================================
# plot scatter
# =============================================================================

fig,[ax1,ax2] = plt.subplots(nrows = 2, gridspec_kw = {'height_ratios':[10,1]}, sharex = True)

ax1.scatter(altered['x'],altered['Percent altered'], alpha = 0.6, lw = 0, s = 11, color = altered['color'])

ax1.set_ylabel('Percent Altered')
ax1.tick_params(bottom = False)

ax1.set_ylim(ax1.get_ylim()[0], 1)

# =============================================================================
# Plot heatmap
# =============================================================================

ax2.bar(hm['x'], [0.8]*len(hm),  color = hm['TP53 status'].astype(str))
ax2.bar(hm['x'], [0.8]*len(hm), bottom = 1, color = hm['High-risk'].astype(str))

ax2.set_xticks(np.arange(len(pt_median)))
ax2.set_xticklabels(list(x_dict.keys()), rotation = 90)

ax2.set_yticks([0.4,1.4])
ax2.set_yticklabels(['TP53 status','High-risk'])

ax2.set_xlim(-1,len(x_dict))

ax2.tick_params(bottom = False, left = False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

fig.tight_layout()
fig.subplots_adjust(hspace = 0.01)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/perPatientPercent altered.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/perPatientPercent altered.png')

# =============================================================================
# 
# =============================================================================

import scipy.stats as stats


sim_no2 = altered.copy()
sim_no2 = sim_no2[~sim_no2['Sample ID'].isin(secondary_tumors)]
sim_no2 = sim_no2[sim_no2['Patient ID'] != 'ID8']

tp53_mut = sim_no2.copy()
tp53_mut = tp53_mut[tp53_mut['Patient ID'].isin(hm[hm['TP53 status'] == 'black']['Patient'].tolist())]
tp53_mut = tp53_mut.groupby('Patient ID').median()['Percent altered']

tp53_intact = sim_no2.copy()
tp53_intact = tp53_intact[tp53_intact['Patient ID'].isin(hm[hm['TP53 status'] != 'black']['Patient'].tolist())]
tp53_intact = tp53_intact.groupby('Patient ID').median()['Percent altered']

print('TP53',stats.mannwhitneyu(tp53_intact, tp53_mut, alternative = 'less'))

high_risk = sim_no2.copy()
high_risk = high_risk[high_risk['Patient ID'].isin(hm[hm['High-risk'] == 'black']['Patient'].tolist())]
high_risk = high_risk.groupby('Patient ID').median()['Percent altered']

low_risk = sim_no2.copy()
low_risk = low_risk[low_risk['Patient ID'].isin(hm[hm['High-risk'] != 'black']['Patient'].tolist())]
low_risk = low_risk.groupby('Patient ID').median()['Percent altered']

print(stats.mannwhitneyu(high_risk, low_risk, alternative = 'greater'))