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
snp_tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/ew_snp_tc.tsv', sep = '\t')

mut_TC.columns = mut_TC.iloc[0]
mut_TC = mut_TC.drop(mut_TC.index[0])
mut_TC['mut_TC'] = mut_TC['mut_TC'].str.split('%').str[0].astype(float) / 100

mut_TC = mut_TC[mut_TC.columns[:6].tolist()]
snp_tc = snp_tc[['Sample ID','snp_TC']]
mut_TC = mut_TC.merge(snp_tc, on = 'Sample ID')

mut_TC = mut_TC[~mut_TC['snp_TC'].isnull()]
mut_TC = mut_TC[mut_TC['Cohort'] == 'M1RP']

fig,ax = plt.subplots()

ax.scatter(mut_TC['mut_TC'],mut_TC['snp_TC'], s = 10, c = 'k', alpha = 0.5)
ax.plot([0,1],[0,1])
ax.set_ylabel('snp_TC')
ax.set_ylim(0,1)
ax.set_xlabel('mut_TC')
ax.set_xlim(0,1)


pos_tc = mut_TC[(mut_TC['mut_TC'] > 0)&(mut_TC['snp_TC'] > 0)]
pos_tc['residual'] = pos_tc['snp_TC'] - pos_tc['mut_TC']

linregress = stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC'])
slope = str(round(linregress[0],2))
r = str(round(linregress[2],2))
std_xy = round(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum()/len(pos_tc['residual'])),2)
std_xy = str(std_xy)

ax.plot([0,1],[linregress[1],linregress[0]+linregress[1]], linestyle = 'dashed')

ax.text(x = 0.01, y = 0.95, s = 'm: %s, r: %s\nstd from y=x: %s\nn=%i' % (slope, r, std_xy, len(mut_TC)))


print(stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC']))
print(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum() / len(pos_tc['residual'])))

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/ew_tc_snpmut_scatter.pdf')
