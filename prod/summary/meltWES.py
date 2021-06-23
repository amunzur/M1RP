# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 13:36:06 2021

@author: amurtha
"""
import pandas as pd
import numpy as np

# =============================================================================
# 
# =============================================================================

df = pd.DataFrame(columns = ['CHR','START','END','Sample ID','Log_ratio','snp_dev'])

for i in np.arange(1,44,1):
    pt = 'ID%i' % i
    cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_%s-segments.txt' % pt, sep = '\t')
    tmp = cn.melt(id_vars = ['CHR','START','END'], var_name = 'Sample ID', value_name = 'val')
    tmp['Log_ratio'] = tmp['val'].str.split(':').str[0].astype(float)
    tmp['snp_dev'] = tmp['val'].str.split(':').str[1].astype(float)
    tmp = tmp.drop('val', axis = 1)
    df = df.append(tmp, ignore_index = False)
    
df.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', index = None, sep = '\t')