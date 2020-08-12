# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:18:56 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# Constants
# =============================================================================

cohort = 'M1RP'

# =============================================================================
# Import mutation, CN, TC data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_allSamples_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)
tc = tc[tc['Cohort'] == cohort]

# =============================================================================
# Separate low TC
# =============================================================================

lowtc = tc[tc['Final tNGS_TC'] < 0.2].copy()
muts_lowtc = muts[muts['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()
cn_lowtc = cn[cn['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()

# =============================================================================
# Eliminate low TC from main dataframes
# =============================================================================

tc = tc[tc['Final tNGS_TC'] >= 0.2]
muts = muts[muts['Sample ID'].isin(tc['Sample ID'])]
cn = cn[cn['Sample ID'].isin(tc['Sample ID'])]


