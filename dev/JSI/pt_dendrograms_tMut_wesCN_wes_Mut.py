# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 13:33:45 2020

@author: amurtha
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage,fcluster


# =============================================================================
# Import various JSI values
# =============================================================================

targeted = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_mutThenCN.tsv', sep = '\t', index_col = 'Unnamed: 0')
wes_cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_wxs_CN.tsv', sep = '\t', index_col = 'Sample ID correct')
wes_mut = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_wxs.tsv', sep = '\t', index_col = 'Unnamed: 0')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)
tc = tc[tc['Cohort'] == 'M1RP']
tc = tc[~tc['Sample ID'].str.contains('cfDNA')]
tc = tc[tc['Final tNGS_TC'] > 0]

# =============================================================================
# Loop over-patients
# =============================================================================

for pt in tc['Patient ID'].unique().tolist():
# for pt in ['ID1']:
    samples = tc[(tc['Patient ID'] == pt)&(tc['Final tNGS_TC'] >= 37.5)]['Sample ID'].copy()
    wes_samples = tc[(tc['Patient ID'] == pt)&(~tc['Sequenza WES_TC'].isnull())]['Sample ID'].copy()
    pt_targeted = targeted[targeted.index.isin(samples)][samples.tolist()]
    pt_wes_cn = wes_cn[wes_cn.index.isin(wes_samples)][wes_samples.tolist()]
    pt_wes_mut = wes_mut[wes_mut.index.isin(wes_samples)][wes_samples.tolist()]
    pt_targeted_wes = targeted[targeted.index.isin(wes_samples)][wes_samples.tolist()]
    
    pt_targeted = pt_targeted.fillna(1)
    pt_wes_cn = pt_wes_cn.fillna(1)
    pt_wes_mut = pt_wes_mut.fillna(1)
    pt_targeted_wes = pt_targeted_wes.fillna(1)
    pt_targeted_wes = pt_targeted_wes.fillna(1)
    
    if len(wes_samples) == 1: continue;
    
    # Set up plots
    # fig, [[ax1,ax2],[ax3,ax4]] = plt.subplots(ncols = 2, nrows = 2)
    fig = plt.figure()
    gs = gridspec.GridSpec(5,2, height_ratios = [1,0.05,0.7,1,0.05])
    
    ax1 = fig.add_subplot(gs[0,0])
    ax1_c = fig.add_subplot(gs[1,0])
    
    ax3 = fig.add_subplot(gs[3,0])
    ax3_c = fig.add_subplot(gs[4,0])
    
    ax2 = fig.add_subplot(gs[0,1])
    ax2_c = fig.add_subplot(gs[1,1])
    ax4 = fig.add_subplot(gs[3,1])
    ax4_c = fig.add_subplot(gs[4,1])
    
    gs.update(hspace = 0)
    
    fig.suptitle(pt)
    
    # Set up and plot targeted dendrogram
    ln = linkage(pt_targeted, method = 'average', metric = 'euclidean')
    
    dn = dendrogram(ln, ax = ax1, color_threshold = 0, above_threshold_color='k')
    samples = pt_targeted.columns.tolist()
    dn_data = dendrogram(ln,no_plot = True)
    
    order = pd.DataFrame({'Sample ID':samples, 'i':np.arange(len(samples))})
    order['x'] = order['i'].apply(lambda x: dn.get('leaves').index(x))
    order_dict = dict(zip(order['Sample ID'],order['x']))
    
    samples = order.sort_values('x')['Sample ID']
    samples = samples.str.split('_').str[2:].str.join('_').tolist()
    
    ax1.yaxis.set_visible(False)
    ax1.tick_params(left = False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_xticklabels(samples, rotation = 90, fontsize = 8)
    [t.set_color('red') for t in ax1.xaxis.get_ticklabels() if t.get_text() in wes_samples.str.split('_').str[2:].str.join('_').tolist()]
    ax1.set_title('JSI from targeted\n(mutation -> copy number)')
    
    # Plot clustering wes muts
    
    cmap = plt.get_cmap('tab10')
    
    k = len(list(set(dn_data.get('color_list'))))
    T = fcluster(ln, k, 'maxclust')
    
    # calculate labels
    for index, row in order.iterrows():
        order.at[index, 'cluster'] = T[row['i']]
        order.at[index, 'cluster_color'] = matplotlib.colors.to_hex(cmap(T[row['i']]))
        
    ax1_c.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax1_c.set_xlim(-0.5,len(samples)-0.5)
    ax1_c.set_yticks([0.66])
    ax1_c.set_ylim(0,1)
    ax1_c.set_yticklabels([''], fontsize = 6)
    ax1_c.spines['left'].set_visible(False)
    ax1_c.spines['bottom'].set_visible(False)
    ax1_c.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    # Set up and plot wes mutation dendrogram
    ln = linkage(pt_wes_mut, method = 'average', metric = 'euclidean')
    
    dn = dendrogram(ln, ax = ax2, color_threshold = 0, above_threshold_color='k')
    samples = pt_wes_mut.columns.tolist()
    dn_data = dendrogram(ln,no_plot = True)
    
    order = pd.DataFrame({'Sample ID':samples, 'i':np.arange(len(samples))})
    order['x'] = order['i'].apply(lambda x: dn.get('leaves').index(x))
    order_dict = dict(zip(order['Sample ID'],order['x']))
    
    samples = order.sort_values('x')['Sample ID']
    samples = samples.str.split('_').str[2:].str.join('_').tolist()
    
    ax2.yaxis.set_visible(False)
    ax2.tick_params(left = False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.set_xticklabels(samples, rotation = 90, fontsize = 8)
    ax2.set_title('JSI from WES mutations')
    
    # Plot clustering wes muts
    
    cmap = plt.get_cmap('tab10')
    
    k = len(list(set(dn_data.get('color_list'))))
    T = fcluster(ln, k, 'maxclust')
    
    # calculate labels
    for index, row in order.iterrows():
        order.at[index, 'cluster'] = T[row['i']]
        order.at[index, 'cluster_color'] = matplotlib.colors.to_hex(cmap(T[row['i']]))
        
    ax2_c.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax2_c.set_xlim(-0.5,len(samples)-0.5)
    ax2_c.set_yticks([0.66])
    ax2_c.set_ylim(0,1)
    ax2_c.set_yticklabels([''], fontsize = 6)
    ax2_c.spines['left'].set_visible(False)
    ax2_c.spines['bottom'].set_visible(False)
    ax2_c.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    # Plot dendrogram of targeted data with only samples that got WES sequencing
    
    ln = linkage(pt_targeted_wes, method = 'average', metric = 'euclidean')
    
    dn = dendrogram(ln, ax = ax3, color_threshold = 0, above_threshold_color='k')
    samples = pt_targeted.columns.tolist()
    dn_data = dendrogram(ln,no_plot = True)
    
    order = pd.DataFrame({'Sample ID':wes_samples, 'i':np.arange(len(wes_samples))})
    order['x'] = order['i'].apply(lambda x: dn.get('leaves').index(x))
    order_dict = dict(zip(order['Sample ID'],order['x']))
    
    samples = order.sort_values('x')['Sample ID']
    samples = samples.str.split('_').str[2:].str.join('_').tolist()
    
    ax3.yaxis.set_visible(False)
    ax3.tick_params(left = False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.set_xticklabels(samples, rotation = 90, fontsize = 8)
    ax3.set_title('JSI from targeted\n(WES samples only)')
    
    # Plot clustering wes muts
    
    cmap = plt.get_cmap('tab10')
    
    k = len(list(set(dn_data.get('color_list'))))
    T = fcluster(ln, k, 'maxclust')
    
    # calculate labels
    for index, row in order.iterrows():
        order.at[index, 'cluster'] = T[row['i']]
        order.at[index, 'cluster_color'] = matplotlib.colors.to_hex(cmap(T[row['i']]))
        
    ax3_c.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax3_c.set_xlim(-0.5,len(samples)-0.5)
    ax3_c.set_yticks([0.66])
    ax3_c.set_ylim(0,1)
    ax3_c.set_yticklabels([''], fontsize = 6)
    ax3_c.spines['left'].set_visible(False)
    ax3_c.spines['bottom'].set_visible(False)
    ax3_c.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    # Set up and plot wes CN dendrogram
    ln = linkage(pt_wes_cn, method = 'average', metric = 'euclidean')
    
    dn = dendrogram(ln, ax = ax4, color_threshold = 0, above_threshold_color='k')
    samples = pt_wes_cn.columns.tolist()
    dn_data = dendrogram(ln,no_plot = True)
    
    order = pd.DataFrame({'Sample ID':samples, 'i':np.arange(len(samples))})
    order['x'] = order['i'].apply(lambda x: dn.get('leaves').index(x))
    order_dict = dict(zip(order['Sample ID'],order['x']))
    
    samples = order.sort_values('x')['Sample ID']
    samples = samples.str.split('_').str[2:].str.join('_').tolist()
    
    ax4.yaxis.set_visible(False)
    ax4.tick_params(left = False)
    ax4.spines['left'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.set_xticklabels(samples, rotation = 90, fontsize = 8)
    ax4.set_title('JSI from WES copy number')
    
    # Plot clustering wes cn
    
    cmap = plt.get_cmap('tab10')
    
    k = len(list(set(dn_data.get('color_list'))))
    T = fcluster(ln, k, 'maxclust')
    
    # calculate labels
    for index, row in order.iterrows():
        order.at[index, 'cluster'] = T[row['i']]
        order.at[index, 'cluster_color'] = matplotlib.colors.to_hex(cmap(T[row['i']]))
        
    ax4_c.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax4_c.set_xlim(-0.5,len(samples)-0.5)
    ax4_c.set_yticks([0.66])
    ax4_c.set_ylim(0,1)
    ax4_c.set_yticklabels([''], fontsize = 6)
    ax4_c.spines['left'].set_visible(False)
    ax4_c.spines['bottom'].set_visible(False)
    ax4_c.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/JSI/pt_dendrograms_byDataType/%s.pdf' % pt)