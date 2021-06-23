# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:13:27 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Import jsi and nmuts
# =============================================================================

jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t', index_col = 'Unnamed: 0')

jsi = jsi.where(np.triu(np.ones(jsi.shape)).astype(bool))
jsi = jsi.reset_index().rename(columns = {'index':'s1'})

jsi = jsi.melt(id_vars = 's1', var_name = 's2', value_name = 'JSI')

jsi = jsi[jsi['s1'] != jsi['s2']]
jsi = jsi[~jsi['JSI'].isna()]

# =============================================================================
# Import n muts
# =============================================================================

nmuts = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts.tsv', sep = '\t', index_col = 'Unnamed: 0')

nmuts = nmuts.where(np.triu(np.ones(nmuts.shape)).astype(bool))
nmuts = nmuts.reset_index().rename(columns = {'index':'s1'})

nmuts = nmuts.melt(id_vars = 's1', var_name = 's2', value_name = 'nmuts')

nmuts = nmuts[nmuts['s1'] != nmuts['s2']]
nmuts = nmuts[~nmuts['nmuts'].isna()]

# =============================================================================
# Merge 2. Keep same patient
# =============================================================================

jsi = jsi.merge(nmuts, on =['s1','s2'])
jsi['p1'] = jsi['s1'].str.split('_').str[1]
jsi['p2'] = jsi['s2'].str.split('_').str[1]

jsi = jsi[jsi['p1'] == jsi['p2']]

jsi = jsi[jsi['p1'] != 'ID8']

# =============================================================================
# Plot scatter
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(jsi['nmuts'], jsi['JSI'])