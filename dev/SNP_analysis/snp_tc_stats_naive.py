# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:23:11 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import math 

mut_TC = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

mut_TC.columns = mut_TC.iloc[0]
mut_TC = mut_TC.drop(mut_TC.index[0])
mut_TC['mut_TC'] = mut_TC['mut_TC'].str.split('%').str[0].astype(float) / 100
mut_TC['snpTC (new)'] = mut_TC['snpTC (new)'].str.split('%').str[0].astype(float) / 100

mut_TC = mut_TC[~mut_TC['snpTC (new)'].isnull()]
mut_TC = mut_TC[mut_TC['Cohort'] == 'M1RP']
mut_TC = mut_TC[(mut_TC['mut_TC'] != 0)&(mut_TC['snpTC (new)'] != 0)]

fig,ax = plt.subplots()

ax.scatter(mut_TC['mut_TC'], mut_TC['snpTC (new)'])
ax.plot([0,1],[0,1])

ax.set_ylabel('snp_TC')
ax.set_ylim(0,1)
ax.set_xlabel('mut_TC')
ax.set_xlim(0,1)

mut_TC['residual'] = mut_TC['snpTC (new)'] - mut_TC['mut_TC']

print(stats.linregress(mut_TC['mut_TC'], mut_TC['snpTC (new)']))
print(math.sqrt(mut_TC['residual'].apply(lambda x: x**2).sum() / len(mut_TC['residual'])))