# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:49:46 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
import scipy.stats as stats
import math
import numpy as np

# =============================================================================
# Pairwise comparison ID10
# =============================================================================


tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc.set_index('Sample ID')

for pt in np.arange(43):
    pt_id = 'ID%i' % (pt+1)
    pt = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_%s-segments.txt' % pt_id, sep = '\t')


    pt = pt[~pt['CHR'].isin(['chrY','chrM'])]
    pt = pt[pt['END']-pt['START'] >= 100000]
    # pt = pt[(pt['START'] != 125121500)&(pt['START'] != 143271500)]
    
    s_list = pt.columns.tolist()[3:]
    s_comb = it.combinations(s_list, 2)
    
    if len(s_list) <= 1:
        continue;
    
    fig,axs = plt.subplots(nrows = len(s_list)-1,ncols = len(s_list)-1)
    
    for i,s in enumerate(s_list):
        pt[s] = pt[s].str.split(':').str[0].astype(float)
        s_tc = tc.at[s,'Final tNGS_TC']
    
    for (s1, s2) in s_comb:
        x = s_list.index(s1)
        y = s_list.index(s2) - 1
        
        s1_tc = tc.at[s1,'Final tNGS_TC']
        s2_tc = tc.at[s2,'Final tNGS_TC']
        
        s1_dd = 2*(math.log2(2-2*s1_tc)-1)
        s2_dd = 2*(math.log2(2-2*s2_tc)-1)
        
        pt_tmp = pt[(pt[s1] > s1_dd) & (pt[s2] > s2_dd)].copy()
        
        if len(s_list) > 2 and y > x:
            axs[y,x].set_visible(False)
        if len(s_list) > 2:
            ax = axs[x,y]
        else:
            ax = axs 
        ax.scatter(pt_tmp[s2],pt_tmp[s1], lw = 0, color = 'k', s = 10)
        lin = stats.linregress(pt_tmp[s1],pt_tmp[s2])
        if x == 0:
            s2 = s2.split('_')[2]
            ax.set_title(s2, fontsize = 6)
        if y == len(s_list) - 2:
            s1 = s1.split('_')[2]
            ax.set_ylabel(s1, fontsize = 6)
            ax.yaxis.set_label_position('right')
        ax.tick_params(labelsize = 6)
        ax.set_xlabel('r=%.2f' % lin[2], fontsize = 6)
        fig.suptitle(pt_id, fontsize=6)
        fig.tight_layout()
    

    