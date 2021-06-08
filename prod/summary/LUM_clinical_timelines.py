# -*- coding: utf-8 -*-
"""
Created on Tue May 25 13:44:38 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime

dpm = 30.436875

met_colors = {'1':'red',
              '2':'green',
              '3':'darkblue'}

treatment_colors = {'1':'yellow'}

mdt_marker = '1'
mdt_x_offset = 1
y_offset = 0.3


# =============================================================================
# Import clinical data
# =============================================================================

clin = pd.read_csv('https://docs.google.com/spreadsheets/d/1hhevDJUzG7wjsuWkbQc_WIJnhZWKDZSqFManz2By6jE/export?format=csv')

clin.columns = clin.iloc[0]
clin = clin.drop(clin.index[0])

clin = clin[clin['Study ID'] != 'LuM-ID10']
clin['y'] = np.arange(0,len(clin),1)

# =============================================================================
# Create plot
# =============================================================================

fig,ax = plt.subplots(figsize = (6,3.5))

# =============================================================================
# Fill RP and fix datetime types
# =============================================================================

clin.loc[clin['Date RP']=='-', 'Date RP'] = clin['Date primary treatment']


clin['Date RP']= pd.to_datetime(clin['Date RP'], format = '%Y-%m-%d')
clin['Date CR1']= pd.to_datetime(clin['Date CR1'], format = '%Y-%m-%d')
clin['Date CR2']= pd.to_datetime(clin['Date CR2'], format = '%Y-%m-%d')
clin['Date CR3']= pd.to_datetime(clin['Date CR3'], format = '%Y-%m-%d')
clin['Date CR4']= pd.to_datetime(clin['Date CR4'], format = '%Y-%m-%d')
clin['Date CR5']= pd.to_datetime(clin['Date CR5'], format = '%Y-%m-%d')
clin['Date CR6']= pd.to_datetime(clin['Date CR6'], format = '%Y-%m-%d')

clin.loc[clin['Date of last FU'].isna(), 'Date of last FU'] = clin['Date CR2']+datetime.timedelta(days=100)

clin['Date of last FU']= pd.to_datetime(clin['Date of last FU'], format = '%Y-%m-%d')
clin['Date lung lesion']= pd.to_datetime(clin['Date lung lesion'], format = '%Y-%m-%d')


# =============================================================================
# Plot first recurrence
# =============================================================================

cr1 = clin[~clin['Date CR1'].isna()].copy()

cr1['date_cr1'] = (cr1['Date CR1'] - cr1['Date RP']).dt.days / dpm
cr1['cr1_color'] = cr1['Metastatic site_C1_code1'].map(met_colors)
cr1['cr1_color2'] = cr1['Metastatic site_C1_code2'].map(met_colors)

ax.scatter(cr1['date_cr1'],cr1['y'], c = cr1['cr1_color'], lw = 0.8, s = 70, edgecolors = 'k', zorder = 100)

# Plot the metastatic directed therapy

cr1_mdt = cr1[cr1['Coded therapy_CR1'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr1_mdt['date_cr1']-mdt_x_offset,cr1_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

# Plot the second color

cr1_2 = cr1[~cr1['cr1_color2'].isna()].copy()

for index, row in cr1_2.iterrows():
    ax.plot(row['date_cr1'], row['y'], fillstyle = 'right', color = row['cr1_color2'], markerfacecoloralt = 'none', markersize = 7.559, lw = 0, marker = 'o', markeredgewidth = 0, zorder = 1000)
    
# Plot treatment post 1st recurrence

cr1_tx = cr1[~cr1['Tx_cr1'].isna()].copy()

for index, row in cr1_tx.iterrows():
    if not pd.isna(row['Time on Tx_cr1']):
        duration = float(row['Time on Tx_cr1'])
    elif not pd.isna(row['Date CR2']):
        duration = (row['Date CR2'] - row['Date CR1']).days / dpm
    else:
        duration = (row['Date of last FU'] - row['Date CR1']).days / dpm
    if row['cr1_multiple_tx'] == '0':   
        c = treatment_colors.get(row['Tx_cr1'])
        ax.plot([row['date_cr1'],row['date_cr1']+duration],[row['y']]*2, c = c, marker = None, lw = 1, zorder = 10, solid_capstyle = 'butt')
    else:
        for i in [-0.1,0.1]:
            y = row['y']+i
            if i < 0:
                c = treatment_colors.get(row['Tx_cr1'])
            else:
                c = treatment_colors.get(row['Tx_cr1_2'])
            ax.plot([row['date_cr1'],row['date_cr1']+duration],[row['y']]*2, c = c, marker = None, lw = 1, zorder = 10, solid_capstyle = 'butt')
            

# =============================================================================
# Plot second recurrence
# =============================================================================

cr2 = clin[~clin['Date CR2'].isna()].copy()

cr2['date_cr2'] = (cr2['Date CR2'] - cr2['Date RP']).dt.days / dpm
cr2['cr2_color'] = cr2['Metastatic site_C2_code1'].map(met_colors)
cr2['cr2_color2'] = cr2['Metastatic site_C2_code2'].map(met_colors)

ax.scatter(cr2['date_cr2'],cr2['y'], c = cr2['cr2_color'], lw = 0.8, s = 70, edgecolors  = 'k', zorder = 100)

# Plot the metastatic directed therapy

cr2_mdt = cr2[cr2['Coded therapy_CR2'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr2_mdt['date_cr2']-mdt_x_offset,cr2_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

# Plot the second color

cr2 = cr2[~cr2['cr2_color2'].isna()]

for index, row in cr2.iterrows():
    ax.plot(row['date_cr2'], row['y'], fillstyle = 'right', color = row['cr2_color2'], markerfacecoloralt = 'none', markersize = 7.559, lw = 0, marker = 'o', markeredgewidth = 0, zorder = 1000)
    
    
# =============================================================================
# Plot third recurrence
# =============================================================================

cr3 = clin[~clin['Date CR3'].isna()].copy()

cr3['date_cr3'] = (cr3['Date CR3'] - cr3['Date RP']).dt.days / dpm
cr3['cr3_color'] = cr3['Metastatic site_C3_code1'].map(met_colors)
cr3['cr3_color2'] = cr3['Metastatic site_C3_code2'].map(met_colors)

ax.scatter(cr3['date_cr3'],cr3['y'], c = cr3['cr3_color'], lw = 0.8, s = 70, edgecolors  = 'k', zorder = 100)

# Plot the metastatic directed therapy

cr3_mdt = cr3[cr3['Coded therapy_CR3'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr3_mdt['date_cr3']-mdt_x_offset,cr3_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

# Plot the second color

cr3 = cr3[~cr3['cr3_color2'].isna()]

for index, row in cr3.iterrows():
    ax.plot(row['date_cr3'], row['y'], fillstyle = 'right', color = row['cr3_color2'], markerfacecoloralt = 'none', markersize = 7.559, lw = 0, marker = 'o', markeredgewidth = 0, zorder = 1000)
    
# =============================================================================
# Plot fourth recurrence
# =============================================================================

cr4 = clin[~clin['Date CR4'].isna()].copy()

cr4['date_cr4'] = (cr4['Date CR4'] - cr4['Date RP']).dt.days / dpm
cr4['cr4_color'] = cr4['Metastatic site_C4_code1'].map(met_colors)

# Plot the metastatic directed therapy

cr4_mdt = cr4[cr4['Coded therapy_CR4'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr4_mdt['date_cr4']-mdt_x_offset,cr4_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

ax.scatter(cr4['date_cr4'],cr4['y'], c = cr4['cr4_color'], lw = 0.8, s = 70, edgecolors  = 'k', zorder = 100)

# =============================================================================
# Plot fifth recurrence
# =============================================================================

cr5 = clin[~clin['Date CR5'].isna()].copy()

cr5['date_cr5'] = (cr5['Date CR5'] - cr5['Date RP']).dt.days / dpm
cr5['cr5_color'] = cr5['Metastatic site_C5_code1'].map(met_colors)

# Plot the metastatic directed therapy

cr5_mdt = cr5[cr5['Coded therapy_CR5'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr5_mdt['date_cr5']-mdt_x_offset,cr5_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

ax.scatter(cr5['date_cr5'],cr5['y'], c = cr5['cr5_color'], lw = 0.8, s = 70, edgecolors  = 'k', zorder = 100)

# =============================================================================
# Plot sixth recurrence
# =============================================================================

cr6 = clin[~clin['Date CR6'].isna()].copy()

cr6['date_cr6'] = (cr6['Date CR6'] - cr6['Date RP']).dt.days / dpm
cr6['cr6_color'] = cr6['Metastatic site_C6_code1'].map(met_colors)
cr6['cr6_color2'] = cr6['Metastatic site_C6_code2'].map(met_colors)

ax.scatter(cr6['date_cr6'],cr6['y'], c = cr6['cr6_color'], lw = 0.8, s = 70, edgecolors  = 'k', zorder = 100)

# Plot the metastatic directed therapy

cr6_mdt = cr6[cr6['Coded therapy_CR6'].str.lower().str.contains('mdt')].copy()
ax.scatter(cr6_mdt['date_cr6']-mdt_x_offset,cr6_mdt['y']-y_offset, lw = 1.2, marker = mdt_marker, c = 'k', edgecolors = 'k')

# Plot second color

cr6 = cr6[~cr6['cr6_color2'].isna()]

for index, row in cr6.iterrows():
    ax.plot(row['date_cr6'], row['y'], fillstyle = 'right', color = row['cr6_color2'], markerfacecoloralt = 'none', markersize = 7.559, lw = 0, marker = 'o', markeredgewidth = 0, zorder = 1000)
    
# =============================================================================
# Plot time to last followup
# =============================================================================

clin['Time2FU'] = (clin['Date of last FU'] - clin['Date RP']).dt.days / dpm

for index, row in clin.iterrows():
    ax.plot([0,row['Time2FU']],[row['y']]*2, c = 'k', zorder = 0, solid_capstyle = 'butt', marker = '>', markevery = [False, True])
    
    
# =============================================================================
# Plot lung resection
# =============================================================================

clin['Date_LR'] = (clin['Date lung lesion'] - clin['Date RP']).dt.days / dpm

ax.scatter(clin['Date_LR'], clin['y']-y_offset, marker = '^', lw = 0.8, edgecolors  = 'k', c = 'yellow', zorder = 10000)

# =============================================================================
# Aethetics 
# =============================================================================

ax.spines['left'].set_visible(False)
ax.tick_params(axis = 'y',left = False, pad = 0)

# Flip Y axis
ax.set_yticks(clin['y'])
ax.set_yticklabels(clin['Study ID'])
ax.invert_yaxis()

ax.set_xlim(-2,156)
ax.set_xticks(np.arange(0,12*14,12))
ax.set_xticklabels(np.arange(0,14,1))
ax.set_xlabel('Years after Prostatectomy')


plt.tight_layout()
fig.savefig('C:/Users/amurtha/Dropbox/2020 - Ghent lung metastases/Figures/Sandbox/Fig1_clinicalTimelines.pdf')