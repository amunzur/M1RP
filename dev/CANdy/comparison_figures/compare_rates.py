# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 10:12:17 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Import copy number data
# =============================================================================

limit_tf = 'High'

nas = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_M1RP_FFPE_cna.tsv', sep = '\t')
mat = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_FFPE_cna.tsv', sep = '\t')
tumor_frac = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tumor_frac.columns = tumor_frac.iloc[0]
tumor_frac = tumor_frac.drop(tumor_frac.index[0])
tumor_frac['Final tNGS_TC'] = tumor_frac['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100

# =============================================================================
# Keep tumor fraction limit
# =============================================================================

if limit_tf == 'All':
    tumor_frac = tumor_frac
elif limit_tf == 'Low':
    tumor_frac = tumor_frac[tumor_frac['Final tNGS_TC'] < 0.4]
elif limit_tf == 'High':
    tumor_frac = tumor_frac[tumor_frac['Final tNGS_TC'] >= 0.4]


mat = mat.merge(tumor_frac[['Sample ID','Final tNGS_TC']], on = 'Sample ID')
nas = nas.merge(tumor_frac[['Sample ID','Final tNGS_TC']], on = 'Sample ID')

# =============================================================================
# Get order
# =============================================================================

cna = mat[mat['Copy_num'] != 0].copy()
cna = cna.groupby(['GENE']).count()[['START']]
cna = cna.sort_values(['START'], ascending=[False]).drop(['START'], axis = 1)

# =============================================================================
# Merge Mat
# =============================================================================

mat = mat.groupby(['GENE','Copy_num']).count().reset_index()[['GENE','Copy_num','START']]

cna = cna.merge(mat[mat['Copy_num'] == -2][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_deepdel'})
cna = cna.merge(mat[mat['Copy_num'] == -1][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_monodel'})
cna = cna.merge(mat[mat['Copy_num'] == 1][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_gain'})
cna = cna.merge(mat[mat['Copy_num'] == 2][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_amp'})

# =============================================================================
# Merge nas
# =============================================================================

nas = nas.groupby(['GENE','Adjusted_copy_num']).count().reset_index()[['GENE','Adjusted_copy_num','START']]

cna = cna.merge(nas[nas['Adjusted_copy_num'] == 0][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'nas_deepdel'})
cna = cna.merge(nas[nas['Adjusted_copy_num'] == 1][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'nas_monodel'})
cna = cna.merge(nas[nas['Adjusted_copy_num'] == 3][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'nas_gain'})
cna = cna.merge(nas[nas['Adjusted_copy_num'] == 4][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'nas_amp'})

cna = cna.fillna(0)

# =============================================================================
# Graph
# =============================================================================

    # 4:'#EE2D24',
    # 3:'#F59496',
    # 2:'#E6E7E8',
    # 1:'#9CC5E9',
    # 0:'#3F60AC',
    # -1:'#ffffff',
    # }

fig,ax = plt.subplots(figsize = (17,10))
w = 0.35

ax.bar(cna.index, cna['mat_deepdel'], width = w, color = '#3F60AC')
bottom = cna['mat_deepdel']
ax.bar(cna.index, cna['mat_monodel'], width = w, color = '#9CC5E9', bottom = bottom)
bottom = bottom + cna['mat_monodel']
ax.bar(cna.index, cna['mat_gain'], width = w, color = '#F59496', bottom = bottom)
bottom = bottom + cna['mat_gain']
ax.bar(cna.index, cna['mat_amp'], width = w, color = '#EE2D24', bottom = bottom)


ax.bar(cna.index+w, cna['nas_deepdel'], width = w, color = '#3F60AC')
bottom = cna['nas_deepdel']
ax.bar(cna.index+w, cna['nas_monodel'], width = w, color = '#9CC5E9', bottom = bottom)
bottom = bottom + cna['nas_monodel']
ax.bar(cna.index+w, cna['nas_gain'], width = w, color = '#F59496', bottom = bottom)
bottom = bottom + cna['nas_gain']
ax.bar(cna.index+w, cna['nas_amp'], width = w, color = '#EE2D24', bottom = bottom)

ax.set_xticks(np.arange(w/2, len(cna), 1))
ax.set_xticklabels(cna.GENE, rotation = 90)
ax.set_xlim(-0.5, len(cna))

if limit_tf == 'All':
    ax.set_ylim(0, 180)
elif limit_tf == 'Low':
    ax.set_ylim(0, 50)
elif limit_tf == 'High':
    ax.set_ylim(0,140)
ax.set_ylabel('Sample count')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/comparison_figures/rates_comparison%sTC.pdf' % limit_tf)