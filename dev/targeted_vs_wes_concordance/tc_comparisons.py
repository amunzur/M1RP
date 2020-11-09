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
# M1B vs M1RP
# =============================================================================
tc = sheet.copy()
tc = tc[['Cohort', 'Final tNGS_TC']]


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

# =============================================================================
# WES vs Targeted 
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




