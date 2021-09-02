# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 14:09:35 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# 
# =============================================================================

segs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', sep = '\t')

# =============================================================================
# 
# =============================================================================

segs['Patient ID'] = segs['Sample ID'].str.split('_').str[1]
segs = segs[~segs['CHR'].isin(['chrY','chrM'])]

# =============================================================================
# 
# =============================================================================

segs = segs.drop_duplicates(['CHR','START','END','Patient ID'])

n_segs = segs.groupby(['CHR','Patient ID']).count()
n_segs = n_segs.reset_index()

n_segs_pCHR = n_segs.groupby(['CHR']).median().reset_index()
n_segs_pCHR = n_segs_pCHR.sort_values(['START'], ascending = False)

x_dict = dict(zip(n_segs_pCHR['CHR'], np.arange(len(n_segs_pCHR))))
n_segs['x'] = n_segs['CHR'].apply(lambda x: x_dict.get(x)+ np.random.uniform(-0.1,0.1))

# =============================================================================
# 
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(n_segs['x'],n_segs['START'], lw = 0, alpha = 0.5)

ax.set_xticks(list(x_dict.values()))
ax.set_xticklabels(list(x_dict.keys()), rotation = 90)
