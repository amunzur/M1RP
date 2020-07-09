# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 13:27:07 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import math

cohort = 'M1RP'

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Variant allele frequency'] = tc['Variant allele frequency'].astype(float)
tc['mut_TC'] = tc['mut_TC'].str.split('%').str[0].astype(float) / 100
tc['snp_TC'] = tc['snp_TC'].str.split('%').str[0].astype(float) / 100

tc = tc[tc['Cohort'] == 'M1RP']
tc = tc[(tc['snp_TC'] != 0)&(tc['mut_TC'] != 0)]
tc['residual'] = tc['snp_TC'] - tc['mut_TC']

pos_tc = tc.copy()

subclonal = tc.copy()
subclonal = subclonal[subclonal['mut_TC'] < .20]
subclonal = subclonal[subclonal['mut_TC'] > 0]
subclonal = subclonal[subclonal['snp_TC'] > .40]

tc['Color'] = 'k'
tc.loc[tc['Sample ID'].isin(subclonal['Sample ID'].tolist()), 'Color'] = 'b'

tc = tc[~tc['snp_TC'].isnull()]
pos_tc = pos_tc[~pos_tc['snp_TC'].isnull()]
pos_tc = pos_tc[~pos_tc['Sample ID'].isin(subclonal['Sample ID'].tolist())]

fig,ax = plt.subplots()

clonal = tc[tc['Color'] == 'k']
subclonal = tc[tc['Color'] == 'b']
ax.scatter(clonal['mut_TC'],clonal['snp_TC'], s = 10, c = clonal['Color'], alpha = 0.5, label = 'Clonal')
ax.scatter(subclonal['mut_TC'],subclonal['snp_TC'], s = 10, c = subclonal['Color'], alpha = 0.5, label = 'Subclonal')
ax.plot([0,1],[0,1])
ax.set_ylabel('snp_TC')
ax.set_ylim(0,1)
ax.set_xlabel('mut_TC')
ax.set_xlim(0,1)

ax.set_title('M1RP')

linregress = stats.linregress(clonal['mut_TC'], clonal['snp_TC'])
slope = str(round(linregress[0],2))
r = str(round(linregress[2],2))
std_xy = round(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum()/len(tc['residual'])),2)
std_xy = str(std_xy)

ax.plot([0,1],[linregress[1],linregress[0]+linregress[1]], linestyle = 'dashed')

ax.text(x = 0.01, y = 0.95, s = 'm: %s, r: %s\nstd from y=x: %s\nn=%i' % (slope, r, std_xy, len(tc)))
ax.legend()


print(stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC']))
print(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum() / len(pos_tc['residual'])))

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_mutsnpTC_scatter_woSubclonal.pdf' % cohort)
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_mutsnpTC_scatter_woSubclonal.png' % cohort, dpi = 500)
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/%s_mutsnpTC_scatter_woSubclonal.pdf' % cohort)