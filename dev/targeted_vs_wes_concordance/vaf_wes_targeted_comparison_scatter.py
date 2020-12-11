# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:32:56 2020

@author: Sarah
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy import stats



sheet = pd.read_csv("https://docs.google.com/spreadsheets/d/1zi7UiPteDA3VU4jtcp1dxET4UGCMKvtLqab-7G13_oU/export?format=csv&gid=265363046")
sheet = sheet[['t_vaf', 'wes_vaf', 'Call status', 'wes_depth', 'Note']]
sheet['t_vaf'] = sheet['t_vaf']*100
sheet['wes_vaf'] = sheet['wes_vaf']*100

sheet = sheet[sheet['wes_vaf'].notna()]
sheet = sheet[sheet['t_vaf'].notna()]
#sheet=sheet.fillna(0)
sheet = sheet[sheet['wes_depth'] > 29]

#Scatterplot
cat = ['WES', 'Both', 'Targeted']        
cat_map = {'WES':'#ffae00', 'Targeted':'#ad1100', 'Both':'#0d6eff'}

fig, ax = plt.subplots()

#ax.scatter(sheet['Targeted VAF'], sheet['WES VAF'], c=sheet['Independent mutation call'].apply(lambda x: cat_map[x]), marker='o', s=2.5)
groups = sheet.groupby('Call status')
for name, group in groups:
    group.plot(ax=ax, kind='scatter', x ='t_vaf', y='wes_vaf', marker='o', s=2.5, label=name, color = cat_map[name])
ax.set_xlabel('Targeted VAF (%)')
ax.set_ylabel('WES VAF (%)')

ax.set_ylim(-2,100)
ax.set_xlim(-2,100)
ax.set_xticks([0,20,40,60,80,100])
ax.set_yticks([0,20,40,60,80,100])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(sheet['t_vaf'], sheet['wes_vaf'])

textstr = " p = 4.566e-108 \n r = 0.965" 
props = dict(boxstyle='square', facecolor='white', alpha=0)
ax.text(0.8, 0.16, textstr, transform=ax.transAxes, fontsize=8, style = 'italic',
        verticalalignment='top', bbox=props)

#Legend
n = sheet['Call status'].value_counts().to_dict()   

#ax.legend(loc="upper left")
labels = ['Independently called in:', 'WES  n=' + str(n.get('WES')), 'Targeted  n=' + str(n.get('Targeted')), 'Both  n=' + str(n.get('Both'))]
handles = [Patch(color = 'none', linewidth = 0), Patch(color = '#ffae00', linewidth = 0),
           Patch(color = '#ad1100',linewidth = 0), Patch(color = '#0d6eff', linewidth = 0)]
ax.legend(handles, labels, loc = 'upper left', prop={'size': 7}, handlelength = 0.8, frameon = False)


fig.savefig('C:\\Users\\Sarah\\Desktop\\VAF_WES_vs_Targeted_read_depth_min_30.png', bbox_inches='tight', dpi = 600)

