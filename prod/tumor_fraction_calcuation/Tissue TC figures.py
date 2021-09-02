# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 13:52:30 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string
from matplotlib.patches import Patch

# =============================================================================
# TC
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc['mut_TC'] = tc['mut_TC'].astype(float)
tc['snp_TC'] = tc['snp_TC'].astype(float)

# =============================================================================
# Remove cfDNA samples and sort by TC
# =============================================================================

tc = tc[~tc['Sample ID'].str.lower().str.contains('cfdna')]
tc = tc[~tc['Sample ID'].str.lower().str.contains('ucc')]

tc = tc.sort_values(['Final tNGS_TC','mut_TC','snp_TC'], ascending = False)

# =============================================================================
# Plot TC methods
# =============================================================================

fig,[ax1,ax2,ax3] = plt.subplots(nrows = 3, sharex = True, figsize = (4,3.5))

# =============================================================================
# Plot final tNGS_TC
# =============================================================================

ax1.scatter(tc['Sample ID'],tc['Final tNGS_TC'], lw = 0, c = 'k', s = 6)

ax1.set_xlim(-4, len(tc)+2)
ax1.tick_params(bottom = False, labelbottom = False, labelsize = 6)

ax1.set_yticks(np.arange(0,1.01,0.25))
ax1.set_ylabel('Final tumor content', fontsize = 6)

# =============================================================================
# Plot mutation tNGS_TC
# =============================================================================

ax2.scatter(tc['Sample ID'],tc['mut_TC'], lw = 0, c = 'k', s = 6)

ax2.set_xlim(-4, len(tc)+2)
ax2.tick_params(bottom = False, labelbottom = False, labelsize = 6)

ax2.set_yticks(np.arange(0,1.01,0.25))
ax2.set_ylabel('Mutation tumor content', fontsize = 6)

# =============================================================================
# Plot SNP  tNGS_TC
# =============================================================================

ax3.scatter(tc['Sample ID'],tc['snp_TC'], lw = 0, c = 'k', s = 6)

ax3.set_xlim(-4, len(tc)+2)
ax3.tick_params(bottom = False, labelbottom = False, labelsize = 6)

ax3.set_yticks(np.arange(0,1.01,0.25))
ax3.set_ylabel('SNP tumor content', fontsize = 6)

# =============================================================================
# Save figure
# =============================================================================

plt.tight_layout()
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_allTissue.pdf')

# =============================================================================
# Plot by tissue location
# =============================================================================

tc['Sample cat'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
tc['Sample cat'] = tc['Sample cat'].replace({'MB':'MT'})

xlabels = []
boxes = []
s_type_n = {}

tc['x'] = 0

for x,cat in enumerate(tc['Sample cat'].unique().tolist()):
    boxes.append(tc[tc['Sample cat'] == cat]['Final tNGS_TC'].tolist())
    xlabels.append('%s\nn=%i'%(cat,len(boxes[-1])))
    tc.loc[tc['Sample cat'] == cat, 'x'] = x+1
    s_type_n[cat] = len(boxes[-1])
    
tc['x'] = tc['x'].apply(lambda x: x + np.random.uniform(-0.15,0.15))

fig,ax = plt.subplots(figsize = (2.25,2))

ax.boxplot(boxes, showfliers = False)
ax.scatter(tc['x'],tc['Final tNGS_TC'], lw = 0, s = 4, c = 'k')

ax.set_xticklabels(xlabels)
ax.set_ylabel('Tumor content', fontsize = 6)
ax.tick_params(labelsize = 6)

plt.tight_layout()
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_bySampleType.pdf')

# =============================================================================
# Plot by sample type propoption negative low and high
# =============================================================================

tc['TC cat'] = 0
tc.loc[(tc['Final tNGS_TC']>0)&(tc['Final tNGS_TC']<0.2), 'TC cat'] = 1
tc.loc[(tc['Final tNGS_TC']>=0.2)&(tc['Final tNGS_TC']<0.4), 'TC cat'] = 2
tc.loc[(tc['Final tNGS_TC']>=0.4), 'TC cat'] = 3

tc = tc.sort_values('TC cat')

tc_cat_counts = tc.groupby('Sample cat').count()
tc_cat_counts = tc_cat_counts[['Cohort']]
tc_cat_counts.columns = ['Total']
tc_cat_counts = tc_cat_counts.reset_index()

tc_grp = tc.groupby(["Sample cat",'TC cat']).count()

tc_grp = tc_grp[['Cohort']]
tc_grp.columns = ['Count']
tc_grp = tc_grp.reset_index()

tc_grp = tc_grp.merge(tc_cat_counts, on = 'Sample cat')
tc_grp['Proportion'] = tc_grp['Count'] / tc_grp['Total']

tc_grp = tc_grp[['Sample cat','TC cat','Proportion']]
tc_grp = tc_grp.sort_values('TC cat')

fig,ax = plt.subplots(figsize = (1.75,2))

bottom = dict(zip(['PB','RP','MLN','MT'],[0,0,0,0]))
color_dict = dict(zip([0,1,2,3],['grey','#fee0d2','#fc9272','#de2d26']))
x_dict = dict(zip(['PB','RP','MLN','MT'],np.arange(4)))

for index, row in tc_grp.iterrows():
    c = color_dict.get(row['TC cat'])
    x = x_dict.get(row['Sample cat'])
    ax.bar(x, row['Proportion'], color = c, bottom = bottom.get(row['Sample cat']))
    bottom[row['Sample cat']] = bottom.get(row['Sample cat']) + row['Proportion']
    
ax.set_xlabel('Sample type', fontsize = 6)
ax.set_xticks(np.arange(4))
ax.set_xticklabels(['PB','RP','MLN','MT'])

ax.set_ylabel("Percent of samples", fontsize = 6)
ax.tick_params(labelsize = 6)

handles = [Patch(color = c) for c in ['grey','#fee0d2','#fc9272','#de2d26']]
labels = ['TC negative','0 - 20%','20 - 40%','> 40%']

# fig.legend(handles, labels, loc = 'upper right', fontsize = 6, handlelength = 0.8)

fig.tight_layout()
# fig.subplots_adjust(right = 0.75)
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_propotion by sample type.pdf')

# =============================================================================
# Plot by tissue location
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

samples = samples[['Sample ID','GGG']]

tc_ggg = tc.merge(samples, on = 'Sample ID', how = 'left')
tc_ggg = tc_ggg[tc_ggg['Sample cat'].isin(['PB','RP'])]

tc_ggg = tc_ggg[(~tc_ggg['GGG'].isna())&(tc_ggg['GGG'] != '?')]

xlabels = []
boxes = []
gleason_n = {}

tc_ggg['x'] = tc_ggg['GGG'].astype(int)
tc_ggg = tc_ggg.sort_values('x')

for x,ggg in enumerate(tc_ggg['GGG'].unique().tolist()):
    boxes.append(tc_ggg[tc_ggg['GGG'] == ggg]['Final tNGS_TC'].tolist())
    xlabels.append('%s\nn=%i'%(ggg,len(boxes[-1])))
    gleason_n[ggg] = len(boxes[-1])
    
tc_ggg['x'] = tc_ggg['x'].apply(lambda x: x + np.random.uniform(-0.15,0.15))

fig,ax = plt.subplots(figsize = (2.5,2))

ax.boxplot(boxes, showfliers = False)
ax.scatter(tc_ggg['x'],tc_ggg['Final tNGS_TC'], lw = 0, s = 4, c = 'k')

ax.set_xticklabels(xlabels)
ax.set_ylabel('Tumor content', fontsize = 6)
ax.tick_params(labelsize = 6)

plt.tight_layout()
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_byGGG.pdf')

# =============================================================================
# 
# =============================================================================

tc_ggg['TC cat'] = 0
tc_ggg.loc[(tc_ggg['Final tNGS_TC']>0)&(tc_ggg['Final tNGS_TC']<0.2), 'TC cat'] = 1
tc_ggg.loc[(tc_ggg['Final tNGS_TC']>=0.2)&(tc_ggg['Final tNGS_TC']<0.4), 'TC cat'] = 2
tc_ggg.loc[(tc_ggg['Final tNGS_TC']>=0.4), 'TC cat'] = 3

tc_ggg = tc_ggg.sort_values('TC cat')

tc_cat_counts = tc_ggg.groupby('GGG').count()
tc_cat_counts = tc_cat_counts[['Cohort']]
tc_cat_counts.columns = ['Total']
tc_cat_counts = tc_cat_counts.reset_index()

tc_grp = tc_ggg.groupby(["GGG",'TC cat']).count()

tc_grp = tc_grp[['Cohort']]
tc_grp.columns = ['Count']
tc_grp = tc_grp.reset_index()

tc_grp = tc_grp.merge(tc_cat_counts, on = 'GGG')
tc_grp['Proportion'] = tc_grp['Count'] / tc_grp['Total']

tc_grp = tc_grp[['GGG','TC cat','Proportion']]
tc_grp = tc_grp.sort_values('TC cat')

fig,ax = plt.subplots(figsize = (2.25,2))

bottom = dict(zip(['1','2','3','4','5'],[0,0,0,0,0]))
color_dict = dict(zip([0,1,2,3],['grey','#fee0d2','#fc9272','#de2d26']))
x_dict = dict(zip(['1','2','3','4','5'],np.arange(5)))

for index, row in tc_grp.iterrows():
    c = color_dict.get(row['TC cat'])
    x = x_dict.get(row['GGG'])
    ax.bar(x, row['Proportion'], color = c, bottom = bottom.get(row['GGG']))
    bottom[row['GGG']] = bottom.get(row['GGG']) + row['Proportion']
    
ax.set_xlabel('GGG', fontsize = 6)
ax.set_xticks(np.arange(5))
ax.set_xticklabels(['1','2','3','4','5'])

ax.set_ylabel("Percent of samples", fontsize = 6)
ax.tick_params(labelsize = 6)

handles = [Patch(color = c) for c in ['grey','#fee0d2','#fc9272','#de2d26']]
labels = ['TC negative','0 - 20%','20 - 40%','> 40%']

fig.legend(handles, labels, loc = 'upper right', fontsize = 6, handlelength = 0.8)

fig.tight_layout()
fig.subplots_adjust(right = 0.75)
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_propotion by GGG.pdf')

del tc_ggg

# =============================================================================
# 
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

samples = samples[['Sample ID','PathTC']]

tc_path = tc.merge(samples, on = 'Sample ID', how = 'left')
tc_path = tc_path[tc_path['Sample cat'].isin(['PB','RP'])]
tc_path = tc_path[~tc_path['PathTC'].isna()]

tc_path['pathTC_int'] = tc_path['PathTC'].str.split('%').str[0].astype(int)
tc_path.loc[tc_path['pathTC_int'] <= 60, 'PathTC'] = '60% or less'

tc_path['TC cat'] = 0
tc_path.loc[(tc_path['Final tNGS_TC']>0)&(tc_path['Final tNGS_TC']<0.2), 'TC cat'] = 1
tc_path.loc[(tc_path['Final tNGS_TC']>=0.2)&(tc_path['Final tNGS_TC']<0.4), 'TC cat'] = 2
tc_path.loc[(tc_path['Final tNGS_TC']>=0.4), 'TC cat'] = 3

tc_path = tc_path.sort_values('TC cat')

tc_cat_counts = tc_path.groupby('PathTC').count()
tc_cat_counts = tc_cat_counts[['Cohort']]
tc_cat_counts.columns = ['Total']
tc_cat_counts = tc_cat_counts.reset_index()

tc_grp = tc_path.groupby(["PathTC",'TC cat']).count()

tc_grp = tc_grp[['Cohort']]
tc_grp.columns = ['Count']
tc_grp = tc_grp.reset_index()

tc_grp = tc_grp.merge(tc_cat_counts, on = 'PathTC')
tc_grp['Proportion'] = tc_grp['Count'] / tc_grp['Total']

tc_grp = tc_grp[['PathTC','TC cat','Proportion']]
tc_grp = tc_grp.sort_values('TC cat')

fig,ax = plt.subplots(figsize = (2.25,2))

bottom = dict(zip(tc_cat_counts['PathTC'].tolist(),list(np.zeros(len(tc_cat_counts)))))
color_dict = dict(zip([0,1,2,3],['grey','#fee0d2','#fc9272','#de2d26']))
x_dict = dict(zip(tc_cat_counts['PathTC'].tolist(),np.arange(len(tc_cat_counts))))

for index, row in tc_grp.iterrows():
    c = color_dict.get(row['TC cat'])
    x = x_dict.get(row['PathTC'])
    ax.bar(x, row['Proportion'], color = c, bottom = bottom.get(row['PathTC']))
    bottom[row['PathTC']] = bottom.get(row['PathTC']) + row['Proportion']
    
ax.set_xlabel('PathTC', fontsize = 6)
ax.set_xticks(np.arange(len(tc_cat_counts)))
ax.set_xticklabels(tc_cat_counts['PathTC'].tolist(),rotation = 90)

ax.set_ylabel("Percent of samples", fontsize = 6)
ax.tick_params(labelsize = 6)

fig.tight_layout()
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_propotion by PathTC.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Tumor content/TC_propotion by PathTC.png', dpi = 500)