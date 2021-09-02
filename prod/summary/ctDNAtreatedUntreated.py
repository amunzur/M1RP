# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 10:58:52 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

# =============================================================================
# Import clinical and TC data
# =============================================================================

clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Get only RP and MLN samples
# =============================================================================

tc['Sample cat'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
tc = tc[tc['Sample cat'].isin(['MLN','RP'])]

tc['color'] = tc['Sample cat'].map({'MLN':'k', 'RP':'b'})

# =============================================================================
# Order patients
# =============================================================================

med_tc = tc.groupby('Patient ID').median()

clin = clin.merge(med_tc[['Final tNGS_TC']], left_on = 'patient_id', right_index = True)

clin = clin.sort_values(['neoadj_tx', 'Final tNGS_TC'], ascending = False)
x_order = dict(zip(clin['patient_id'].tolist(), np.arange(len(clin))))

tc['x'] = tc['Patient ID'].map(x_order)

# =============================================================================
# scatter
# =============================================================================

neo = tc[tc['Patient ID'].isin(clin[clin['neoadj_tx'] == 1]['patient_id'].tolist())]['Final tNGS_TC']
tx_naive = tc[tc['Patient ID'].isin(clin[clin['neoadj_tx'] != 1]['patient_id'].tolist())]['Final tNGS_TC']


fig,ax = plt.subplots(figsize = (1.5,2))

ax.boxplot([neo, tx_naive], showfliers = False)

xs = pd.Series([1]*len(neo)).apply(lambda x: x + np.random.uniform(-0.1,0.1)).tolist() + \
    pd.Series([2]*len(tx_naive)).apply(lambda x: x + np.random.uniform(-0.1,0.1)).tolist()
ax.scatter(xs, neo.tolist() + tx_naive.tolist(), lw = 0, alpha = 0.8, s = 5)

ax.set_xticks([1,2])
ax.set_xticklabels(['Treated','Untreated'], fontsize = 5)
ax.set_ylabel('Tumor content', fontsize = 6)

ax.tick_params(labelsize = 6)

fig.tight_layout()