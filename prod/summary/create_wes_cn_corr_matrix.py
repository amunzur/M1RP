# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 13:54:34 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools as it
import multiprocessing as mp
import math
from joblib import Parallel, delayed

# =============================================================================
# 
# =============================================================================

def get_corr(s1, s2, cn, tc):
    cn1 = cn[cn['Sample ID'] == s1]
    cn2 = cn[cn['Sample ID'] == s2]
    segs = seg_dict[s1.split('_')[1]][s2.split('_')[1]].copy()
    for index, row in segs.iterrows():
        c1_tmp = cn1[(cn1['CHR'] == row['CHR'])&(cn1['START']<=row['START'])&(cn1['END']>=row['END'])].reset_index(drop = True)
        segs.at[index, s1] = c1_tmp.at[0,'Log_ratio']
        c2_tmp = cn2[(cn2['CHR'] == row['CHR'])&(cn2['START']<=row['START'])&(cn2['END']>=row['END'])].reset_index(drop = True)
        segs.at[index, s2] = c2_tmp.at[0,'Log_ratio']
    
    segs = segs[~segs['CHR'].isin(['chrY','chrM'])]
    segs = segs[segs['END']-segs['START'] >= 100000]
    
    s1_tc = tc.at[s1,'Final tNGS_TC']
    s2_tc = tc.at[s2,'Final tNGS_TC']
    
    s1_dd = 2*(math.log2(2-2*s1_tc)-1)
    s2_dd = 2*(math.log2(2-2*s2_tc)-1)
    
    segs = segs[(segs[s1] > s1_dd) & (segs[s2] > s2_dd)].copy()
    
    segs = segs.set_index(['CHR','START','END'])
    corr_matrix = segs.corr('pearson')
    
    return corr_matrix.at[s1,s2]


# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Keep WES samples
# =============================================================================

tc = tc[tc['Sample ID'].isin(cn['Sample ID'].tolist())]
tc = tc[tc['Final tNGS_TC'] > 0.1]

cn = cn[cn['Sample ID'].isin(tc['Sample ID'].tolist())]

samples = cn['Sample ID'].unique().tolist()
pts = tc['Patient ID'].unique().tolist()

# =============================================================================
# Get segment dictionary
# =============================================================================

seg_dict = {}

for c in it.combinations_with_replacement(pts, 2):
    p1 = c[0]
    p2 = c[1]
    if p1 not in list(seg_dict.keys()):
        seg_dict[p1] = {}
    segs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome segmentation intersections/%s_%s.bed' % (p1, p2), sep = '\t', header = None, names = ['CHR','START','END'])
    segs['CHR'] = 'chr'+segs['CHR'].astype(str)
    tmp = seg_dict.get(p1)
    tmp[p2] = segs
    seg_dict[p1] = tmp

# =============================================================================
# CN samples
# =============================================================================

tc = tc.set_index('Sample ID')
matrix = pd.DataFrame(index = samples, columns = samples)

def func(s1, s2):
    return get_corr(s1, s2, cn, tc)

# results = [func(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2)]

n_cores = mp.cpu_count() 
results = Parallel(n_jobs=n_cores)(delayed(func)(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2))

for (s1,s2),r in zip(it.combinations_with_replacement(samples, 2), results):
    matrix.at[s1,s2] = r
    matrix.at[s2,s1] = r

matrix.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_wesCNA_pearson_corr.tsv', sep = '\t')
