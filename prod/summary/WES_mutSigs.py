# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:20:28 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

mut_sigs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/all_sbs_weights.tsv', sep = '\t')

mut_sigs['patient'] = mut_sigs['Sample'].str.split('_').str[1]
pt_s_counts = mut_sigs.groupby('patient').count()['Sample']
pts = pt_s_counts.index.tolist()

mut_sigs = mut_sigs.drop('Unnamed: 31', axis = 1)

mut_sigs['MMRd'] = mut_sigs['Signature_6'] + mut_sigs['Signature_15'] + mut_sigs['Signature_20'] + mut_sigs['Signature_26']
mut_sigs['BRCA2'] = mut_sigs['Signature_3']
mut_sigs['APOBEC'] = mut_sigs['Signature_2'] + mut_sigs['Signature_13']
mut_sigs['Ageing'] = mut_sigs['Signature_1']

#Quick exploration of the data before incorporating into main oncoprint
mut_sigs = mut_sigs[mut_sigs['patient'].isin(pts)]
mut_sigs = mut_sigs.set_index(['Sample', 'patient'])
mut_sigs = mut_sigs[['Ageing','MMRd','BRCA2','APOBEC',]]
mut_sigs['Other'] = 1 - mut_sigs.sum(axis=1)
mut_sigs = mut_sigs.reset_index(drop=False)
mut_sigs = mut_sigs.loc[:, (mut_sigs != 0).any(axis=0)]
#Boxplot
def get_cmap(n, name='nipy_spectral'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
cmap = get_cmap(len(mut_sigs))

#Barplot
fig = plt.figure(figsize=(20,4))
gs  = gridspec.GridSpec(2,len(pts), height_ratios = [0.25,1], width_ratios=pt_s_counts.tolist())
ax_l = fig.add_subplot(gs[0, :])
mut_sigs.plot.bar(x='Sample', cmap=cmap, stacked=True, ax=ax_l, width=0)
axs = []
for i, pt in enumerate(pts):
    sig_pt = mut_sigs[mut_sigs['patient'] == pt].copy()
    sig_pt = sig_pt[sig_pt.apply(lambda r: (r != 0).any(), axis=1)]
    if i == 0:
        ax = fig.add_subplot(gs[1,i])
        ax1 = ax
    else:
        ax = fig.add_subplot(gs[1,i],sharey = ax1)
    sig_pt.plot.bar(x='Sample', cmap=cmap, stacked=True, ax=ax, width=0.9)
    ax.get_legend().remove()
    ax.set_xlabel(None)
    axs.append(ax)
for i, ax in enumerate(axs):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    if i == 0:
        ax.tick_params(axis='both', labelsize=6, left=True, labelleft=True, bottom=True, labelbottom=True, length=3, width=0.5)
    else:
        ax.tick_params(axis='both', labelsize=6, left=True, labelleft=False, bottom=True, labelbottom=True, length=3, width=0.5)


ax_l.legend(loc='center', ncol=4, fontsize = 6)
ax_l.spines['left'].set_visible(False)
ax_l.spines['bottom'].set_visible(False)
ax_l.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
# fig.subplots_adjust(left=0.05, top=0.65, right=0.98, bottom=0.35)
plt.tight_layout()
plt.savefig("C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/mut_sigs_wes.pdf")
#plt.clf()
