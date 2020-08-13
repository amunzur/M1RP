# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:18:56 2020

@author: amurtha
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage,fcluster

# =============================================================================
# Constants
# =============================================================================

matplotlib.rcParams['lines.linewidth'] = 0.5

cohort = 'M1RP'
pt = 'ID1'
cn_color_dict = {-2:'#3F60AC',-1:'#9CC5E9',0:'#E6E7E8',1:'#F59496',2:'#EE2D24'}
mut_color_dict = {'Frameshift': '#FFC907', 'Missense':'#79B443','Stopgain': '#FFC907', 'Splice': '#BD4398', 'Non-frameshift': '#a9a9a9'}

# =============================================================================
# Helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False))]


# =============================================================================
# Import mutation, CN, TC, JSI data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_allSamples_cna.tsv', sep = '\t')
jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/jsi_matrix.tsv', sep = '\t', index_col='Unnamed: 0')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Cohort'] == cohort]

muts = keepCodingMutations(muts)

# =============================================================================
# Separate low TC
# =============================================================================

lowtc = tc[tc['Final tNGS_TC'] < 0.2].copy()
muts_lowtc = muts[muts['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()
cn_lowtc = cn[cn['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()

# =============================================================================
# Eliminate low TC from main dataframes
# =============================================================================

tc = tc[tc['Final tNGS_TC'] >= 0.2]
muts = muts[muts['Sample ID'].isin(tc['Sample ID'])]
cn = cn[cn['Sample ID'].isin(tc['Sample ID'])]
muts = muts[muts['GENE'].isin(cn['GENE'].unique().tolist())]

# =============================================================================
# Assign color to copy number and mutation
# =============================================================================

cn['color'] = cn['Copy_num'].apply(lambda x: cn_color_dict.get(x))
muts['color'] = muts['EFFECT'].apply(lambda x: mut_color_dict.get(x.split(' ')[0]))

for pt in tc['Patient ID'].unique().tolist():
    # =============================================================================
    # Get sample order
    # =============================================================================
    print(pt)
    pt_samples = tc[tc['Patient ID'] == pt]['Sample ID'].tolist()
    matrix_tmp = jsi[jsi.index.isin(pt_samples)][pt_samples].copy().fillna(1)
    
    pt_cn = cn[cn['Sample ID'].isin(pt_samples)].copy()
    pt_muts = muts[muts['Sample ID'].isin(pt_samples)].copy()
    
    # =============================================================================
    # Set up figure
    # =============================================================================
    
    fig = plt.figure(figsize = (8.5,11))
    gs = gridspec.GridSpec(7,1, height_ratios=[1,0.1,10,0.1,0.4,0.01,0.1])
    
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[4,0])
    ax4 = fig.add_subplot(gs[2,0])
    ax5 = fig.add_subplot(gs[6,0])
    
    gs.update(hspace=0.005)
    
    # =============================================================================
    # Plot dendrogram()
    # =============================================================================
    
    ln = linkage(matrix_tmp, method = 'average', metric = 'euclidean')
    
    dn = dendrogram(ln, ax = ax1, color_threshold = 0, above_threshold_color='k')
    samples = matrix_tmp.columns.tolist()
    dn_data = dendrogram(ln,no_plot = True)
    
    order = pd.DataFrame({'Sample ID':samples, 'i':np.arange(len(samples))})
    order['x'] = order['i'].apply(lambda x: dn.get('leaves').index(x))
    order_dict = dict(zip(order['Sample ID'],order['x']))
    
    samples = order.sort_values('x')['Sample ID'].tolist()
    
    ax1.yaxis.set_visible(False)
    ax1.xaxis.set_visible(False)
    ax1.tick_params(left = False, bottom = False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    
    # =============================================================================
    # Plot clusters
    # =============================================================================
    
    cmap = plt.get_cmap('tab10')
    
    k = len(list(set(dn_data.get('color_list'))))
    T = fcluster(ln, k, 'maxclust')
    
    # calculate labels
    for index, row in order.iterrows():
        order.at[index, 'cluster'] = T[row['i']]
        order.at[index, 'cluster_color'] = matplotlib.colors.to_hex(cmap(T[row['i']]))
        
    order = order.merge(tc[['Sample ID', 'Final tNGS_TC']], on = 'Sample ID')
    order['tc_color'] = 'red'
    order.loc[order['Final tNGS_TC'] < 0.4, 'tc_color'] = 'grey'
    
    ax2.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax2.set_xlim(-0.5,len(samples)-0.5)
    ax2.set_yticks([0.66])
    ax2.set_ylim(0,1)
    ax2.set_yticklabels(['Cluster'], fontsize = 8)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    ## Make another 
    ax5.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax5.set_xlim(-0.5,len(samples)-0.5)
    ax5.set_yticks([0.66])
    ax5.set_ylim(0,1)
    ax5.set_yticklabels([''], fontsize = 8)
    ax5.spines['left'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    # =============================================================================
    # Plot ctDNA%
    # =============================================================================
    
    ax3.bar(order['x'], order['Final tNGS_TC'], color = order['tc_color'])
    ax3.plot([-0.5,max(order['x'])+0.5],[0.4,0.4], color = 'k', lw = 0.5, ls = 'dashed')
    
    ax3.set_xlim(-0.5,len(samples)-0.5)
    ax3.set_ylim(0,1)
    
    ax3.tick_params(bottom = False, labelbottom = True)
    
    ax3.set_ylabel('Tumor\ncontent', fontsize = 8)
    ax3.set_yticks([0,0.4,1])
    ax3.set_yticklabels([0,0.4,1], fontsize = 8)
    ax3.set_xticks(np.arange(0,len(samples),1))
    ax3.set_xticklabels(samples, rotation = 90)
    
    # =============================================================================
    # Merge order (x-coord) on mutations and cn.
    # =============================================================================
    
    pt_cn['x'] = pt_cn['Sample ID'].apply(lambda x: order_dict.get(x))
    pt_muts['x'] = pt_muts['Sample ID'].apply(lambda x: order_dict.get(x))
    
    gene_list = cn['GENE'].unique().tolist()
    gene_dict = dict(zip(gene_list, np.arange(len(gene_list))))
    
    pt_cn['y'] = pt_cn['GENE'].apply(lambda x: gene_dict.get(x))
    pt_muts['y'] = pt_muts['GENE'].apply(lambda x: gene_dict.get(x) + 0.4)
    
    # =============================================================================
    # Plot copy number
    # =============================================================================
    
    ax4.bar(pt_cn['x'], 0.8, bottom = pt_cn['y'], color = pt_cn['color'], zorder = 10)
    ax4.scatter(pt_muts['x'], pt_muts['y'], color = pt_muts['color'], marker = 's', s = 6, zorder = 100)
    
    ax4.set_yticks(np.arange(0.4,len(gene_list),1))
    ax4.set_xticks(np.arange(len(samples)))
    ax4.set_yticklabels(gene_list, fontsize = 8)
    
    ax4.set_xlim(-0.5,len(samples)-0.5)
    ax4.set_ylim(0,len(gene_list))
    ax4.invert_yaxis()
    
    ax4.tick_params(left = False, bottom = False, labelbottom = False)
    ax4.spines['left'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    
    plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/oncoprints/%s_%s.pdf' % (cohort, pt),bbox_inches='tight')