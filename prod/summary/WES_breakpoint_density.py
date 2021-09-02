# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 15:03:39 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib.patches import Polygon

# =============================================================================
# 
# =============================================================================

chr_size = pd.read_excel('G:/Andy Murtha/panelDesign/chr_size.xlsx')
chr_size = chr_size.set_index('CHR')

acen = pd.read_csv('G:/Andy Murtha/panelDesign/cytoBandhg38.txt', sep = '\t', header = None, names = ['CHR','START','END','BAND','cat'])

all_segs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', sep = '\t')

all_segs['Patient ID'] = all_segs['Sample ID'].str.split('_').str[1]
all_segs = all_segs[all_segs['CHR'] != 'chrM']
all_segs = all_segs.drop('Sample ID', axis = 1)

# =============================================================================
# Drop duplicate segements for same patient samples
# =============================================================================

all_segs = all_segs.drop_duplicates(['CHR','START','END','Patient ID'])

# =============================================================================
# Drop breakpoints within the centromere or 100kbp in either direction
# =============================================================================

acen = acen[acen['cat'] == 'acen']

# for index, row in acen.iterrows():
#     c = row['CHR']
#     s = row['START']
#     e = row['END']
#     all_segws = all_segs[(all_segs['CHR'] != c)|((all_segs['END'] < s-100000)|(all_segs['END'] > s+100000))]

# =============================================================================
# Drop "max" endpoint
# =============================================================================

max_segs = all_segs.groupby(['Patient ID','CHR']).max().reset_index().set_index(['Patient ID','CHR','END'])

all_segs = all_segs[~all_segs.set_index(['Patient ID','CHR','END']).index.isin(max_segs.index)]

# =============================================================================
# Create figure
# =============================================================================

fig,axs = plt.subplots(ncols = 24, gridspec_kw={'width_ratios':chr_size['Total length (bp)']}, figsize = (15,2))

# =============================================================================
# plot breakpoints
# =============================================================================

for ax,c in zip(axs, chr_size.index.tolist()):
    ax.spines['left'].set_visible(False)
    ax.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
    ax.set_xlim(0,chr_size.at[c,'Total length (bp)'])
    ax.set_xlabel(c, fontsize = 6, rotation = 90)
    chr_segs = all_segs[all_segs['CHR'] == c].copy()
    x_kde = stats.gaussian_kde(chr_segs['END'])
    xs = np.linspace(0,chr_size.at[c,'Total length (bp)'], chr_size.at[c,'Total length (bp)']//10000)
    x_kde.covariance_factor = lambda : .15
    x_kde._compute_covariance()
    ax.plot(xs,x_kde(xs), lw = 0.5)
    ax.fill_between(xs, 0, x_kde(xs))
    ax.set_ylim(bottom = 0)
    c_acen = acen[acen['CHR'] == c].copy()
    acen_start = c_acen['START'].min()
    acen_end = c_acen['END'].max()
    y_top = ax.get_ylim()[1]
    coords = np.array([(acen_start,0),(acen_end,0),(acen_end,y_top),(acen_start,y_top)])
    p = Polygon(coords, color = 'k')
    ax.add_patch(p)


plt.tight_layout()
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/breakpoint_distribution.pdf')