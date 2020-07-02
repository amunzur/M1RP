# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:46:11 2020

@author: amurtha
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 10:12:17 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# =============================================================================
# Import copy number data
# =============================================================================


if len(sys.argv) > 1:
    alpha = sys.argv[1]
    limit_tf = sys.argv[2]
else:
    alpha = '0.01'
    limit_tf = 'High'

nas = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/abi_enza/Targeted_sequencing_error_correction/CANdy_Abienza_ctDNA_cna.tsv', sep = '\t')

# =============================================================================
# Keep tumor fraction limit and highest sample from each patient
# =============================================================================

nas = nas.sort_values('ctDNA%', ascending = False)
nas = nas.drop_duplicates(['Patient ID','GENE'])

if limit_tf == 'All':
    nas = nas
elif limit_tf == 'Low':
    nas = nas[nas['ctDNA%'] < 0.4]
elif limit_tf == 'High':
    nas = nas[nas['ctDNA%'] >= 0.4]

print(len(nas['Patient ID'].unique()))
# =============================================================================
# eliminate unevaluable samples
# =============================================================================

mat = nas.copy()
nas = nas[nas['Adjusted_copy_num'] != -1]
nas_count_genes = nas.groupby('GENE').count().reset_index()[['GENE','START']]
nas_count_genes.columns = ['GENE','count']

mat = mat[mat['ctDNA%'] >= 0.19]
print(len(mat['Patient ID'].unique()))
mat_count_genes = mat.groupby('GENE').count().reset_index()[['GENE','START']]
mat_count_genes.columns = ['GENE','count']

# =============================================================================
# Get order
# =============================================================================

cna = mat[mat['Copy_num'] != 0].copy()
cna = cna.groupby(['GENE']).count()[['START']]
cna = cna.sort_values(['START'], ascending=[False]).drop(['START'], axis = 1)

# =============================================================================
# Merge nas
# =============================================================================

mat = mat.groupby(['GENE','Copy_num']).count().reset_index()[['GENE','Copy_num','START']]
mat = mat.merge(mat_count_genes, on = 'GENE', how = 'left')

mat['START'] = mat['START']/mat['count']
mat = mat.drop('count', axis = 1)

cna = cna.merge(mat[mat['Copy_num'] == -2][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_deepdel'})
cna = cna.merge(mat[mat['Copy_num'] == -1][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_monodel'})
cna = cna.merge(mat[mat['Copy_num'] == 1][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_gain'})
cna = cna.merge(mat[mat['Copy_num'] == 2][['GENE','START']], on = 'GENE', how = 'left').rename(columns = {'START':'mat_amp'})

# =============================================================================
# Merge nas
# =============================================================================

nas = nas.groupby(['GENE','Adjusted_copy_num']).count().reset_index()[['GENE','Adjusted_copy_num','START']]

nas = nas.merge(nas_count_genes, on = 'GENE', how = 'left')

nas['START'] = nas['START']/nas['count']
nas = nas.drop('count', axis = 1)


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

ax.set_ylim(0, 1)
ax.set_ylabel('CNA %')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/comparison_figures/pt_rates_comparison%sTC_%s.pdf' % (limit_tf, alpha))