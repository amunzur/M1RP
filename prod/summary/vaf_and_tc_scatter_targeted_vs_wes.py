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

#---------------------------------------------------------------------------------------
#Scatterplot black and white, where wes depth is minimum of 30 (Elie's cutoff for calling a mutation in wes)

fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot()

sheet = sheet[sheet['wes_depth'] > 29]

ax.scatter(x =sheet['t_vaf'], y=sheet['wes_vaf'], marker='o', s=1.5, zorder=10, color='k')
ax.set_xlabel('Mutant Allele (%)\n(Targeted Gene Panel)')
ax.set_ylabel('Mutant Allele (%)\n(Whole Exome Panel)')

ax.set_ylim(0,100)
ax.set_xlim(0,100)
ax.set_xticks([0,20,40,60,80,100])
ax.set_yticks([0,20,40,60,80,100])
ax.grid(color='#bdbdbd',ls='--', lw=0.5, alpha=0.3, dashes=(5, 4), zorder=0)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

x = np.linspace(*ax.get_xlim())
ax.plot(x, x, color='red', lw=1.2, ls='--',dashes=(7.5, 5), zorder=1)

#Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(sheet['t_vaf'], sheet['wes_vaf'])

textstr = "Pearson\nr = " + str(r_value)[:4] 
props = dict(boxstyle='square', facecolor='white', alpha=0)
ax.text(0.10, 0.86, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)



fig.savefig('C:\\Users\\Sarah\\Desktop\\VAF_WES_vs_Targeted_wes_depth_min_30.pdf', bbox_inches='tight', dpi = 600)


#---------------------------------------------------------------------------------------
#Scatterplot coloured WES vs targeted separated by colour, if mutation detected in one or both seq
cat = ['WES', 'Both', 'Targeted']        
cat_map = {'WES':'#ffae00', 'Targeted':'#ad1100', 'Both':'#0d6eff'}

fig, ax = plt.subplots()

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



################################################################################################
#############################################################################################

#Plot scatter for TC 
tc_est = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022")
tc_est.columns = tc_est.iloc[0]
tc_est = tc_est.drop(tc_est.index[0]) 

tc_est = tc_est[tc_est['Sequenza WES_TC'].notna()]      
tc_est['Sequenza WES_TC'] = tc_est['Sequenza WES_TC'].str[:-1].astype(np.float32)
tc_est['Final tNGS_TC'] = tc_est['Final tNGS_TC'].str[:-1].astype(np.float32)
        

#Plot scatter
fig, ax = plt.subplots()
ax.scatter(x=tc_est['Final tNGS_TC'], y=tc_est['Sequenza WES_TC'], marker='o', s=2.5, color = '#0d6eff')

#Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(tc_est['Final tNGS_TC'], tc_est['Sequenza WES_TC'])

#P and r value text
textstr = " n=" + str(len(tc_est)) + "\n p = 4.84e-16 \n r = 0.739"
props = dict(boxstyle='square', facecolor='white', alpha=0)
ax.text(0.8, 0.17, textstr, transform=ax.transAxes, fontsize=8, style = 'italic',
        verticalalignment='top', bbox=props)

#Formatting
ax.set_xlabel('Targeted TC estimate (%)')
ax.set_ylabel('WES TC estimate (%)')

ax.set_ylim(-1,100)
ax.set_xlim(-1,100)
ax.set_xticks([0,20,40,60,80,100])
ax.set_yticks([0,20,40,60,80,100])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig('C:\\Users\\Sarah\\Desktop\\tcscatter.pdf', bbox_extra_artists=(), bbox_inches='tight')


################################################################################################
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        