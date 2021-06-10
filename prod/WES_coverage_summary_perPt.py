# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:53:25 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
from matplotlib.patches import Patch


def p(c):
    return Patch(color = c)

c_dict = dict(zip( ['1','2','3','4'],['green','yellow','red','black']))

# =============================================================================
# Import wes
# =============================================================================

wes = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/exome_sequencing_metrics_June2021.tsv', sep = '\t')

wes['Patient ID'] = wes['SAMPLE'].str.split('_').str[1]

wes['rating'] = 1
wes.loc[~wes['Additional depth needed (90% 8 mutant reads)'].isna(), 'rating'] = 2
wes.loc[~wes['Additional depth needed (median 8 mutant reads)'].isna(), 'rating'] = 3
wes.loc[wes['Median mutant reads'] == 0, 'rating'] = 4

# =============================================================================
# 
# =============================================================================

df = pd.DataFrame(index = wes['Patient ID'].unique().tolist())

wes = wes.groupby(['Patient ID','rating']).count().reset_index()
wes = wes[['Patient ID','rating','SAMPLE']]
wes.columns = ['Patient ID','rating','Count']

for i in [1,2,3,4]:
    tmp = wes[wes['rating'] == i].copy().set_index('Patient ID')
    tmp = tmp[['Count']]
    tmp.columns = [str(i)]
    df = df.merge(tmp, left_index = True, right_index = True, how = 'left').fillna(0)
    
# =============================================================================
# 
# =============================================================================

df = df.sort_values(['1','2','3','4'], ascending = False)

fig,ax = plt.subplots()

df['bottom'] = 0

for i in ['1','2','3','4']:
    ax.bar(df.index, df[i], bottom = df['bottom'], color = c_dict.get(i))
    df['bottom'] = df['bottom'] + df[i]
    
ax.set_xticklabels(df.index, rotation = 90, fontsize = 8)

ax.set_ylabel('Number of samples')

handles = [p(c) for c in ['green','yellow','red','black']]
labels = ['>50% mutation called', '<50% mutation called', '<10% mutations called', 'TC-']

ax.legend(handles, labels)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Extraction and Sequencing/Figures/WES_coverage_summary_perPt.pdf')