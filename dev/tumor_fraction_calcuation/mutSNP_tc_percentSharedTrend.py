# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 10:44:36 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# =============================================================================
# Import SNP and mut TC data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['mut_TC'] = tc['mut_TC'].str.split('%').str[0].astype(float) / 100
tc['snp_TC'] = tc['snp_TC'].str.split('%').str[0].astype(float) / 100
tc['Position'] = tc['Position'].astype(float)
tc = tc[tc['Cohort'] == 'M1RP']
tc = tc[tc['mut_TC'] > 0]
tc = tc[~tc['snp_TC'].isnull()]

# tc = tc[['Cohort','Sample ID','Patient ID','mut_TC','snp_TC']]

# =============================================================================
# Import mutation and merge sample relatedness
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations.tsv', sep = '\t')

sample_counts = tc[['Sample ID','Patient ID']].drop_duplicates().groupby('Patient ID').count().reset_index()
mutation_counts = muts.groupby(['Patient ID','CHROM','POSITION','GENE','EFFECT']).count()[['Cohort']].reset_index()

mutation_counts = mutation_counts.merge(sample_counts, on = 'Patient ID', how = 'left')

mutation_counts = mutation_counts.rename(columns = {'Cohort':'Mutation count','Sample ID':'Sample count','CHROM':'Chromosome','GENE':'Gene','POSITION':'Position','EFFECT':'Effect'})

mutation_counts['Percent'] = mutation_counts['Mutation count'] / mutation_counts['Sample count']

tc = tc.merge(mutation_counts, on = ['Patient ID','Chromosome','Position','Gene','Effect'], how = 'left' )

# =============================================================================
# Create r value for mutations shared in at least x % of samples
# =============================================================================

rs = []

for x in np.arange(0,1,0.01):
    tmp = tc[tc['Percent'] >= x].copy()
    r = stats.linregress(tmp['mut_TC'],tmp['snp_TC'])[2]
    rs.append(r)
    del tmp
    
# =============================================================================
# Plot r values as % increases
# =============================================================================

fig, ax = plt.subplots()

ax.scatter(np.arange(0,1,0.01), rs)
ax.set_ylabel('r')
ax.set_xlabel('Minimum % of samples with mutation')

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/tumor_fraction_calcuation/mutSNP_tc_percentSharedTrend.pdf')