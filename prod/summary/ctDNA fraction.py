# -*- coding: utf-8 -*-
"""
Created on Mon May 10 15:30:01 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

mapping = {'Tx naive':1,
           'Progression':3,
           'On treatment':4,
           'Stable disease':5}

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/cfDNA timepoints.xlsx')

preOpRxNaive = tc[(tc['Status'] == 'Tx naive')&(tc['Pre Op'] == True)]['Final tNGS_TC'].tolist()
postOpRxNaive = tc[(tc['Status'] == 'Tx naive')&(tc['Pre Op'] == False)]['Final tNGS_TC'].tolist()

progression = tc[tc['Status'] == 'Progression']['Final tNGS_TC'].tolist()
on_treatment = tc[tc['Status'] == 'On treatment']['Final tNGS_TC'].tolist()
stable_disease = tc[tc['Status'] == 'Stable disease']['Final tNGS_TC'].tolist()

tc['x'] = tc['Status'].map(mapping)
tc.loc[(tc['Status'] == 'Tx naive')&(tc['Pre Op'] == False), 'x'] = 2

tc['x'] = tc['x'].apply(lambda x: x+np.random.uniform(-0.15,0.15))
tc['color'] = tc['Final tNGS_TC'].map({0:'grey'}).fillna('k')

# =============================================================================
# Plot some boxes
# =============================================================================

boxes = [preOpRxNaive,postOpRxNaive,progression,on_treatment,stable_disease]

fig,ax = plt.subplots(figsize = (3.25,1.75))

ax.boxplot(boxes, showfliers = False, zorder = 10)

tc_0 = tc[tc['Final tNGS_TC'] == 0]
tc_1 = tc[tc['Final tNGS_TC'] > 0]

ax.scatter(tc_0['x'],tc_0['Final tNGS_TC'], lw = 0, color = tc_0['color'], s = 10, alpha = 0.6, zorder = 100)
ax.scatter(tc_1['x'],tc_1['Final tNGS_TC'], lw = 0, color = tc_1['color'], s = 10, alpha = 0.9, zorder = 100)


ax.set_ylabel('ctDNA fraction', fontsize = 6)
ax.set_ylim(top = 1)

ax.set_xticklabels(['Pre-op\nRx. Naive\nn=%i'%len(preOpRxNaive),
                    'Post-op\nRx. Naive\nn=%i'%len(postOpRxNaive),
                    'Progression\nn=%i'%len(progression),
                    'On Rx.\nn=%i'%len(on_treatment),
                    'Stable disease\noff treatment\nn=%i'%len(stable_disease)], fontsize = 6)


ax.set_yticks(np.arange(0,1.1,0.25))
ax.tick_params(labelsize = 6)

fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/ctDNA fraction/ctDNAfraction_byTimepoint.pdf')