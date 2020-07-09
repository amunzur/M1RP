# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:15:34 2020

@author: amurtha
"""

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
seq_stats = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=878831703')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Variant allele frequency'] = tc['Variant allele frequency'].astype(float)
tc['mut_TC'] = tc['mut_TC'].str.split('%').str[0].astype(float) / 100
tc['snp_TC'] = tc['snp_TC'].str.split('%').str[0].astype(float) / 100

tc = tc[tc['Cohort'] == 'M1RP']

seq_stats = seq_stats[seq_stats['Cohort'] == 'M1RP']
seq_stats = seq_stats[seq_stats['Targeted or WES'] == 'Targeted (73 gene)']
seq_stats = seq_stats[~seq_stats['Sample ID'].str.contains('cfDNA')]
seq_stats = seq_stats[~seq_stats['Sample ID'].str.contains('gDNA')]
seq_stats = seq_stats[~seq_stats['Sample ID'].str.contains('WBC')]

pos_tc = tc.copy()
pos_tc = pos_tc[(pos_tc['snp_TC'] != 0)&(pos_tc['mut_TC'] != 0)]

pos_tc['residual'] = pos_tc['snp_TC'] - pos_tc['mut_TC']
pos_tc = pos_tc[~pos_tc['snp_TC'].isnull()]

low_depth = pos_tc.copy()
low_depth = low_depth.merge(seq_stats[['Sample ID','Median coverage']], on = 'Sample ID')
low_depth = low_depth[low_depth['Median coverage'] < low_depth['Median coverage'].quantile(.5)]

high_depth = pos_tc.copy()
high_depth = high_depth.merge(seq_stats[['Sample ID','Median coverage']], on = 'Sample ID')
high_depth = high_depth[high_depth['Median coverage'] >= high_depth['Median coverage'].quantile(.5)]

stat, p = stats.mannwhitneyu(low_depth['residual'],high_depth['residual'],alternative='two-sided')


fig,ax = plt.subplots()

ax.scatter(low_depth['mut_TC'],low_depth['snp_TC'], s = 10, c = 'b', alpha = 0.5, label = 'Depth < median')
ax.scatter(high_depth['mut_TC'],high_depth['snp_TC'], s = 10, c = 'k', alpha = 0.5, label = 'Depth > median')

ax.plot([0,1],[0,1])
ax.set_ylabel('snp_TC')
ax.set_ylim(0,1)
ax.set_xlabel('mut_TC')
ax.set_xlim(0,1)

ax.set_title('M1RP')

linregress = stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC'])
all_r = str(round(linregress[2],2))

linregress = stats.linregress(low_depth['mut_TC'], low_depth['snp_TC'])
low_r = str(round(linregress[2],2))

linregress = stats.linregress(high_depth['mut_TC'], high_depth['snp_TC'])
high_r = str(round(linregress[2],2))
p = str(round(p, 3))

ax.plot([0,1],[linregress[1],linregress[0]+linregress[1]], linestyle = 'dashed')

ax.text(x = 0.1, y = 1, s = 'All: r = %s\nHigh: r = %s\nLow: r = %s\nMWU p = %s' % (all_r, high_r, low_r,p), va = 'top')
ax.legend()

print(stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC']))
print(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum() / len(pos_tc['residual'])))

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_mutsnpTC_scatter_wolow_depth.pdf' % cohort)
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_mutsnpTC_scatter_wolow_depth.png' % cohort, dpi = 500)
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/%s_mutsnpTC_scatter_wolow_depth.pdf' % cohort)