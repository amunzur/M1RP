# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 14:32:35 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

cohort = 'M1RP'
tc_cutoff = 37.5

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)
tc = tc[tc['Cohort'] == cohort]

tc = tc[~tc['Patient ID'].isin(['ID32','ID30','ID20','ID26'])]

tc['color'] = 'red'
tc.loc[tc['Final tNGS_TC'] < tc_cutoff, 'color'] = 'grey'

grouped = tc.groupby(['Patient ID','color']).count().reset_index()
total = tc.groupby('Patient ID').count().reset_index()[['Patient ID','Sample ID']]
total.columns = ['Patient ID', 'Sample count']

grouped = grouped.merge(total, on = 'Patient ID')
grouped['Sample ID'] = grouped['Sample ID'] / grouped['Sample count'] * 100

all_patients = grouped[['Patient ID']].drop_duplicates()

high = grouped[grouped['color'] == 'red']
high = high.merge(all_patients, how = 'right', on = 'Patient ID')
low = grouped[grouped['color'] == 'grey']
low = low.merge(all_patients, how = 'right', on = 'Patient ID')

# =============================================================================
# Order samples
# =============================================================================

low = low.sort_values('Sample ID',ascending = False)
low['order'] = np.arange(len(low))
high = high.merge(low[['Patient ID','order']], on = 'Patient ID')
tc = tc.merge(low[['Patient ID','order']], on = 'Patient ID')

tc = tc.sort_values('order', ascending = False)
high = high.sort_values('order', ascending = True)

# =============================================================================
# Create plot
# =============================================================================

fig,[ax1,ax2] = plt.subplots(nrows = 2, sharex = True)

# =============================================================================
# plot scatter
# =============================================================================

ax1.scatter(tc['order'], tc['Final tNGS_TC'], c = tc['color'], s = 8)

ax1.tick_params(axis = 'x', bottom = False)
ax1.set_ylabel('Tumor content')

# =============================================================================
# Plot bars
# =============================================================================

ax2.bar(high['order'], high['Sample ID'], color = 'red')
ax2.bar(low['order'], low['Sample ID'], bottom = high['Sample ID'], color = 'grey')

ax2.tick_params(axis = 'x', bottom = False, rotation = 90)
ax2.set_xticks(np.arange(len(high)))
ax2.set_xticklabels(high['Patient ID'])
ax2.set_ylabel('Percent of samples')

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/tumor_fraction_calcuation/low_tumorContent_figure.png')