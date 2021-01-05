# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:54:50 2020

@author: Sarah
"""


import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats



plt.rcParams.update({'font.size': 12})

# =============================================================================
# TC from google sheets
# =============================================================================
sheet = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022")
sheet.columns = sheet.iloc[0]
sheet = sheet.drop(0)

sheet['Final tNGS_TC'] = sheet['Final tNGS_TC'].str[:-1].astype(np.float32)
sheet['Sequenza WES_TC'] = sheet['Sequenza WES_TC'].str[:-1].astype(np.float32)


# =============================================================================
# WES vs Targeted Scatterplot
# =============================================================================
tc = sheet.copy()

tc = tc[['Final tNGS_TC', 'Sequenza WES_TC']]
tc = tc.dropna()

#Make figure
fig = plt.figure(figsize=(6,5))
ax1 = fig.add_subplot()

ax1.set_title('Targeted and WES TC estimates \n M1RP')

ax1.plot(tc["Final tNGS_TC"],tc["Sequenza WES_TC"],'o', markersize = 2)

ax1.set_xlabel('Targeted Seq TC Estimate (%)')
ax1.set_ylabel('WES TC Estimate (%)')
ax1.set_ylim(0,100)
ax1.set_xlim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)



fig.savefig('C:\\Users\\Sarah\\Desktop\\Targeted_WES_compare_TC_scatter.png', dpi=600, bbox_inches='tight')



# =============================================================================
# Separate by cohort and sample type
# =============================================================================
tc = sheet.copy()
tc = tc[['Cohort', 'Sample ID', 'Final tNGS_TC', 'Sequenza WES_TC']]

tc['type'] = ''
tc = tc.set_index('Sample ID', drop = False)

for index, row in tc.iterrows():
        #Sample type matrix row
        sample = row['Sample ID']
        if '_MLN' in sample:
            tc.at[sample, 'type'] = 'MLN'
        elif '_MB' in sample:
            tc.at[sample, 'type'] = 'MB'     
        elif '_RP' in sample:
            tc.at[sample, 'type'] = 'RP'
        elif 'PB' in sample:
            tc.at[sample, 'type'] = 'PB'            
        elif 'cfDNA' in sample:            
            tc.at[sample, 'type'] = 'cfDNA'
            
m1rp = tc[tc['Cohort'] == 'M1RP']            
m1b = tc[tc['Cohort'] == 'M1B']            
      
# M1RP plot
fig = plt.figure(figsize=(5.5,5))
ax1 = fig.add_subplot()

ax1.set_title('M1RP TC estimates\nTargeted Sequencing')

sns.stripplot(x="type",y="Final tNGS_TC",data=m1rp, palette=['#A6CEE3', '#B2DF8A', '#33A02C', '#FD5E53', '#1F78B4'],ax=ax1,zorder=0,dodge=True, size = 3)
sns.boxplot(x="type",y="Final tNGS_TC",data=m1rp, whis=np.inf,ax=ax1,color='grey',boxprops=dict(alpha=.1),zorder=10,dodge=True, width=0.7, linewidth = 1)

ax1.set_ylabel('Tumour Content %')
ax1.set_xlabel('Sample Type')
ax1.set_ylim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# M1B plot
fig = plt.figure(figsize=(4,5))
ax1 = fig.add_subplot()

ax1.set_title('M1B TC estimates\nTargeted Sequencing')

sns.stripplot(x="type",y="Final tNGS_TC",data=m1b, palette=['#B2DF8A', '#FD5E53'],ax=ax1,zorder=0,dodge=True, size = 3)
sns.boxplot(x="type",y="Final tNGS_TC",data=m1b, showfliers = False, whis=np.inf,ax=ax1,color='grey',boxprops=dict(alpha=.1),zorder=10,dodge=True, width=0.5, linewidth = 1)

ax1.set_ylabel('Tumour Content %')
ax1.set_xlabel('Sample Type')
ax1.set_ylim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Paired t test
m1b = tc[tc['Cohort'] == 'M1B']['Final tNGS_TC']
m1rp = tc[tc['Cohort'] == 'M1RP']['Final tNGS_TC']
tstat, pvalue = stats.ttest_ind(m1rp, m1b)


# MannwhitneyU
stat, p = stats.mannwhitneyu(m1b, m1rp)

fig.savefig('C:\\Users\\Sarah\\Desktop\\M1B_TC_Boxplot.png', bbox_inches='tight', dpi = 600)

#fig.savefig('C:\\Users\\Sarah\\Desktop\\M1RPvsM1B_TC_boxplot.png', bbox_inches='tight', dpi = 600)





######################################################################################################
#ARCHIVED

# =============================================================================
# WES vs Targeted Boxplot
# =============================================================================
tc = sheet.copy()

tc = tc[['Final tNGS_TC', 'Sequenza WES_TC']]
tc = tc.dropna()

target = tc[['Final tNGS_TC']]
wes = tc[['Sequenza WES_TC']]

target.columns = ['TC']
wes.columns = ['TC']

target['Type'] = 'Targeted'
wes['Type'] = 'WES'

tstat, pvalue = stats.ttest_rel(target['TC'], wes['TC'])
stat, p = stats.mannwhitneyu(target['TC'], wes['TC'])

tc = target.append(wes)

#Make figure
fig = plt.figure(figsize=(4,5))
ax1 = fig.add_subplot()

ax1.set_title('Targeted and WES TC estimates \n M1RP')

sns.stripplot(x="Type",y="TC",data=tc,palette="Set2",ax=ax1,zorder=0,dodge=True, size = 3)
sns.boxplot(x="Type",y="TC",data=tc, whis=np.inf,ax=ax1,color='grey',boxprops=dict(alpha=.1),zorder=10,dodge=True, width=0.5, linewidth = 1)

ax1.set_ylabel('Tumour Content %')
ax1.set_ylim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)



fig.savefig('C:\\Users\\Sarah\\Desktop\\Targeted_WES_compare_TC_boxplot.png', dpi=600, bbox_inches='tight')



# =============================================================================
# Separate by cohort and sample type
# =============================================================================
tc = sheet.copy()
tc = tc[['Cohort', 'Final tNGS_TC', 'Sequenza WES_TC']]


#Make figure
fig = plt.figure(figsize=(4,5))
ax1 = fig.add_subplot()

ax1.set_title('M1RP and M1B TC estimates \n Final tNGS_TC')

sns.stripplot(x="Cohort",y="Final tNGS_TC",data=tc,palette="Set2",ax=ax1,zorder=0,dodge=True, size = 3)
sns.boxplot(x="Cohort",y="Final tNGS_TC",data=tc, whis=np.inf,ax=ax1,color='grey',boxprops=dict(alpha=.1),zorder=10,dodge=True, width=0.5, linewidth = 1)

ax1.set_ylabel('Tumour Content %')
ax1.set_ylim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Paired t test
m1b = tc[tc['Cohort'] == 'M1B']['Final tNGS_TC']
m1rp = tc[tc['Cohort'] == 'M1RP']['Final tNGS_TC']
tstat, pvalue = stats.ttest_ind(m1rp, m1b)


# MannwhitneyU
stat, p = stats.mannwhitneyu(m1b, m1rp)


fig.savefig('C:\\Users\\Sarah\\Desktop\\M1RPvsM1B_TC_boxplot.png', bbox_inches='tight', dpi = 600)


