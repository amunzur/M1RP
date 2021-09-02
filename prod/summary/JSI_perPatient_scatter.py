# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 17:11:04 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string
import scipy.stats as stats

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# 
# =============================================================================

matrix = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t', index_col = 'Unnamed: 0')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

matrix = matrix.where(np.triu(np.ones(matrix.shape)).astype(bool))
matrix = matrix.reset_index().rename(columns = {'index':'s1'})

sim = matrix.melt(id_vars = 's1', var_name = 's2', value_name = 'JSI')

sim = sim[~sim['JSI'].isna()]
sim = sim[sim['s1'] != sim['s2']]

sim['p1'] = sim['s1'].str.split('_').str[1]
sim['p2'] = sim['s2'].str.split('_').str[1]

sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's1', right_on = 'Sample ID')
sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's2', right_on = 'Sample ID')

sim = sim[(sim['Final tNGS_TC_x'] > 0.1)&(sim['Final tNGS_TC_y'] > 0.1)]

# =============================================================================
# Get boxplots for same patient vs other patient
# =============================================================================

sim = sim[sim['p1'] == sim['p2']].copy()

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC','M1RP_ID30_RP3']

sim = sim[~sim['p1'].isin(['ID9','ID20'])]

sim['color'] = 'k'
sim.loc[sim['p1'] == 'ID8', 'color'] = 'orange'
sim.loc[(sim['s1'].isin(secondary_tumors))|(sim['s2'].isin(secondary_tumors)), 'color'] = 'blue'

# =============================================================================
# Import clinical data and TP53 mutation status
# =============================================================================

pts = sim['p1'].unique().tolist()

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
hm.loc[hm['Patient'].isin(['ID8','ID30','ID19','ID15']), 'Secondary tumor'] = 0

# =============================================================================
# Get order
# =============================================================================

pt_median = sim.groupby('p1').median().sort_values('JSI', ascending = False)[['JSI']]
hm = hm.merge(pt_median, left_on = 'Patient', right_index = True)
hm = hm.sort_values(['Secondary tumor','TP53 status','JSI'], ascending = False)

x_dict = dict(zip(hm['Patient'].tolist(), np.arange(len(hm))))

sim['x'] = sim['p1'].map(x_dict)
hm['x'] = hm['Patient'].map(x_dict)

color_dict = {0:'grey', 1:'black'}
hm['TP53 status'] = hm['TP53 status'].map(color_dict)
hm['High-risk'] = hm['High-risk'].map(color_dict)

# =============================================================================
# plot scatter
# =============================================================================

fig,[ax1,ax2] = plt.subplots(nrows = 2, gridspec_kw = {'height_ratios':[10,1]}, sharex = True)

ax1.scatter(sim['x'],sim['JSI'], alpha = 0.6, lw = 0, s = 11, color = sim['color'])

ax1.set_ylabel('JSI')
ax1.tick_params(bottom = False)

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

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/perPatientJSI.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/perPatientJSI.png')

# =============================================================================
# 
# =============================================================================

sim_no2 = sim.copy()
sim_no2 = sim_no2[(~sim_no2['s1'].isin(secondary_tumors))&(~sim_no2['s2'].isin(secondary_tumors))]
sim_no2 = sim_no2[sim_no2['p1'] != 'ID8']

tp53_mut = sim_no2.copy()
tp53_mut = tp53_mut[tp53_mut['p1'].isin(hm[hm['TP53 status'] == 'black']['Patient'].tolist())]
tp53_mut = tp53_mut.groupby('p1').std()['JSI']

tp53_intact = sim_no2.copy()
tp53_intact = tp53_intact[tp53_intact['p1'].isin(hm[hm['TP53 status'] != 'black']['Patient'].tolist())]
tp53_intact = tp53_intact.groupby('p1').std()['JSI']

print('TP53',stats.mannwhitneyu(tp53_intact, tp53_mut, alternative = 'greater'))

high_risk = sim_no2.copy()
high_risk = high_risk[high_risk['p1'].isin(hm[hm['High-risk'] == 'black']['Patient'].tolist())]
high_risk = high_risk.groupby('p1').std()['JSI']

low_risk = sim_no2.copy()
low_risk = low_risk[low_risk['p1'].isin(hm[hm['High-risk'] != 'black']['Patient'].tolist())]
low_risk = low_risk.groupby('p1').std()['JSI']

print(stats.mannwhitneyu(high_risk, low_risk))