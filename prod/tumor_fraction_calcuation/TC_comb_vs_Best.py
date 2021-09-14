# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 11:45:37 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Import TC and mutation data
# =============================================================================

tc_comb = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/tumor_fraction/tumor_content.tsv', sep = '\t')
muts_comb = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/mutations/melted.tsv', sep = '\t')
cn_comb = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/copy_num/gene_cna.ffpe.tsv', sep = '\t')

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)


samples['PathTC'] = samples['PathTC'].str.split('%').str[0].astype(float)

# =============================================================================
# 
# =============================================================================

samples = samples[(samples['GGG'] != '?')&(~samples['GGG'].isna())]
samples['GGG'] = samples['GGG'].astype(float)

samples = samples[samples['Sample Category'].isin(['RP','PB'])]
samples = samples[samples['PathTC'] >= 50]

max_ggg = samples[['Patient ID','Sample Category','GGG']].groupby(['Patient ID','Sample Category']).max()
max_ggg = max_ggg.reset_index()[['Patient ID','Sample Category','GGG']]
max_ggg.columns = ['Patient ID','Sample Category','Max GGG']
samples = samples.merge(max_ggg, on = ['Patient ID','Sample Category'])

max_pathTC = samples[['Patient ID','Sample Category','PathTC']].groupby('Patient ID').max()
max_pathTC = max_pathTC.reset_index()[['Patient ID','Sample Category','PathTC']]
max_pathTC.columns = ['Patient ID','Sample Category','Max PathTC']
samples = samples.merge(max_pathTC, on = ['Patient ID','Sample Category'])

samples = samples[(samples['Max PathTC'] == samples['PathTC'])&(samples['GGG'] == samples['Max GGG'])]


