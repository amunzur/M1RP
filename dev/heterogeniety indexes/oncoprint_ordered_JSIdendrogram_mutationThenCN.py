# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:27:58 2020

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
matplotlib.rcParams['hatch.linewidth'] = 0.25

cohort = 'M1RP'
tc_cutoff = 0.375

cn_color_dict = {-2:'#3F60AC',-1:'#9CC5E9',0:'#E6E7E8',1:'#F59496',2:'#EE2D24'}
mut_color_dict = {'Frameshift': '#FFC907', 'Missense':'#79B443','Stopgain': '#FFC907', 'Splice': '#BD4398', 'Non-frameshift': '#a9a9a9','Other':'grey'}
coding = ['Missense','Stopgain','Spice','Frameshift','Non-framshift']

# =============================================================================
# Helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False))]

def get_missing_chr(cn):
    missing_chr = pd.DataFrame({'CHROM':np.arange(1,23,1)})
    missing_chr = missing_chr[~missing_chr['CHROM'].isin(cn['CHROM'])]
    chr_sub = pd.DataFrame({'CHROM':np.arange(1,23,1)})
    for index, row in chr_sub.iterrows():
        sub = len(missing_chr[missing_chr['CHROM'] < row['CHROM']])
        chr_sub.at[index, 'sub'] = sub
    return dict(zip(chr_sub['CHROM'], chr_sub['sub']))

def get_ytick_lables(cn):
    cn = cn[['GENE','CHROM']].drop_duplicates().reset_index(drop = True)
    length = len(cn)
    blank_row = pd.DataFrame({'GENE': [''], 'CHROM': [np.nan]}, index = [0])
    c = 0
    for index, row in cn.iterrows():
        if index+1 < length and row['CHROM'] != cn.iloc[index+1+c]['CHROM'] and cn.iloc[index+c]['CHROM'] != '':
            dfs = np.split(cn, [index+c+1])
            cn = pd.concat([dfs[0], blank_row, dfs[1]], ignore_index=True)
            c+=1
    cn['GENE'] = cn['GENE'].apply(lambda x: "$\it{%s}$" % x)
    return cn['GENE'].tolist();

def plot_mutations(pt_muts, ax):
    muts_i = pt_muts[pt_muts['Independent'] == True].copy()
    muts_ni = pt_muts[pt_muts['Independent'] == False].copy()
    if len(muts_i) > 0:
        ax.scatter(muts_i['x'], muts_i['y'], color = muts_i['color'], marker = 's', s = 12, zorder = 100, lw = 0)
    if len(muts_ni) > 0:
        ax.scatter(muts_ni['x'], muts_ni['y'], facecolor = muts_ni['color'], marker = 's', s = 12, zorder = 100, hatch = '/////////////////', lw = 0, edgecolor = 'w')
    

# =============================================================================
# Import mutation, CN, TC, JSI data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_allSamples_cna.tsv', sep = '\t')
jsi_muts = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/jsi_matrix_mutThenCN.tsv', sep = '\t', index_col='Unnamed: 0')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Cohort'] == cohort]


cn=cn.rename(columns = {'CHROMOSOME':'CHROM'})
muts = muts[muts['GENE'].isin(cn['GENE'].unique().tolist())]
muts['CHROM'] = muts['CHROM'].replace('chrX','chr22')
cn['CHROM'] = cn['CHROM'].replace('chrX','chr22')
muts['CHROM'] = muts['CHROM'].str.extract('(\d+)').astype(int)
cn['CHROM'] = cn['CHROM'].str.extract('(\d+)').astype(int)

# =============================================================================
# Assign color to copy number and mutation
# =============================================================================

muts['EFFECT'] = muts['EFFECT'].str.split(' ').str[0]
muts.loc[~muts['EFFECT'].isin(coding), 'EFFECT'] = 'Other'

cn['color'] = cn['Copy_num'].apply(lambda x: cn_color_dict.get(x))
muts['color'] = muts['EFFECT'].apply(lambda x: mut_color_dict.get(x.split(' ')[0]))

tc['tc_color'] = 'red'
tc.loc[tc['Final tNGS_TC'] < tc_cutoff, 'tc_color'] = 'grey'

# =============================================================================
# Separate low TC
# =============================================================================

lowtc = tc[tc['Final tNGS_TC'] < tc_cutoff].copy()
lowtc = lowtc.sort_values('Final tNGS_TC', ascending = False)
muts_lowtc = muts[muts['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()
cn_lowtc = cn[cn['Sample ID'].isin(lowtc['Sample ID'].tolist())].copy()

# =============================================================================
# Eliminate low TC from main dataframes
# =============================================================================

tc = tc[tc['Final tNGS_TC'] >= tc_cutoff]
muts = muts[muts['Sample ID'].isin(tc['Sample ID'])]
cn = cn[cn['Sample ID'].isin(tc['Sample ID'])]

# =============================================================================
# For loop to plot entire figure
# =============================================================================

# for pt in ['ID4']:
for pt in tc['Patient ID'].unique().tolist():
    # =============================================================================
    # Get sample order
    # =============================================================================
    print(pt)
    if pt=='ID20':continue;
    pt_samples = tc[tc['Patient ID'] == pt]['Sample ID'].tolist()
    pts_ltc = lowtc[lowtc['Patient ID'] == pt].copy()
    pt_lowtc_samples = pts_ltc['Sample ID'].tolist()
    matrix_tmp = jsi_muts[jsi_muts.index.isin(pt_samples)][pt_samples].copy().fillna(1)
    
    pt_cn = cn[cn['Sample ID'].isin(pt_samples)].copy()
    pt_muts = muts[muts['Sample ID'].isin(pt_samples)].copy()
    
    pt_ltc_cn = cn_lowtc[cn_lowtc['Sample ID'].isin(pt_lowtc_samples)].copy()
    pt_ltc_muts = muts_lowtc[muts_lowtc['Sample ID'].isin(pt_lowtc_samples)].copy()
    
    # =============================================================================
    # Set up figure
    # =============================================================================
    fig_width = 5/8 + .25 * (len(pt_samples)+len(pt_lowtc_samples)+1)
    fig = plt.figure(figsize = (fig_width,9))
    gs = gridspec.GridSpec(7,4,
                           height_ratios=[1,0.1,10,0.01,0.4,0.01,0.1],
                           width_ratios=[.1,len(pt_samples),0.1,len(pt_lowtc_samples)])
    
    ax0 = fig.add_subplot(gs[2,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,1])
    ax3 = fig.add_subplot(gs[4,1])
    ax4 = fig.add_subplot(gs[2,1])
    ax5 = fig.add_subplot(gs[6,1])
    ax6 = fig.add_subplot(gs[2,3])
    ax7 = fig.add_subplot(gs[1,3])
    ax8 = fig.add_subplot(gs[4,3], sharey = ax3)
    
    gs.update(hspace=0.005, wspace = 0.01)
    
    # =============================================================================
    # Get gene list
    # =============================================================================
    
    genes = pt_cn[pt_cn['Copy_num'] != 0]['GENE']
    genes = genes.append(pt_ltc_cn[pt_ltc_cn['Copy_num'] != 0]['GENE']).unique().tolist()
    genes = genes + pt_muts['GENE'].append(pt_ltc_muts['GENE']).unique().tolist()
    genes = list(set(genes))
    
    pt_cn = pt_cn[pt_cn['GENE'].isin(genes)]
    pt_ltc_cn = pt_ltc_cn[pt_ltc_cn['GENE'].isin(genes)]
    
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
        
    order = order.merge(tc[['Sample ID', 'Final tNGS_TC', 'tc_color']], on = 'Sample ID')

    
    ax2.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax2.set_xlim(-0.5,len(samples)-0.5)
    ax2.set_yticks([0.66])
    ax2.set_ylim(0,1)
    ax2.set_yticklabels([''], fontsize = 6)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    ## Make another 
    ax5.bar(order['x'], 0.67, bottom = 0.33, color = order['cluster_color'])
    
    ax5.set_xlim(-0.5,len(samples)-0.5)
    ax5.set_yticks([0.66])
    ax5.set_ylim(0,1)
    ax5.set_yticklabels([''], fontsize = 6)
    ax5.spines['left'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.tick_params(left = False, bottom = False, labelleft = True, labelbottom = False)
    
    # =============================================================================
    # Plot ctDNA%
    # =============================================================================
    
    ax3.bar(order['x'], order['Final tNGS_TC'], color = order['tc_color'])
    ax3.plot([-0.5,max(order['x'])+0.5],[tc_cutoff,tc_cutoff], color = 'k', lw = 0.5, ls = 'dashed')
    
    ax3.set_xlim(-0.5,len(samples)-0.5)
    ax3.set_ylim(0,1)
    
    ax3.tick_params(bottom = False, labelbottom = True)
    
    ax3.set_ylabel('Tumor\ncontent', fontsize = 6)
    ax3.set_yticks([0,tc_cutoff,1])
    ax3.set_yticklabels([0,tc_cutoff,1], fontsize = 6)
    ax3.set_xticks(np.arange(0,len(samples),1))
    ax3.set_xticklabels(pd.Series(samples).str.split('_').str[2:].str.join('_'), rotation  = 90, fontsize = 6)
    
    # =============================================================================
    # Merge order (x-coord) on mutations and cn.
    # =============================================================================
    
    pt_cn['x'] = pt_cn.apply(lambda x: order_dict.get(x['Sample ID']), axis = 1)
    if len(pt_muts)>0:
        pt_muts['x']=pt_muts.apply(lambda x: order_dict.get(x['Sample ID']), axis = 1)
    else:
        pt_muts['x'] = 0
    
    gene_list = pt_cn.append(pt_ltc_cn).sort_values(['CHROM', 'START'])['GENE'].unique().tolist()
    gene_dict = dict(zip(gene_list, np.arange(len(gene_list))))
    
    missing_chr = get_missing_chr(pt_cn.append(pt_ltc_cn))
    
    pt_cn['y'] = pt_cn.apply(lambda x: gene_dict.get(x['GENE']) + x["CHROM"]-(missing_chr.get(x['CHROM']) + 1), axis=1)
    if len(pt_muts) > 0:
        pt_muts['y'] = pt_muts.apply(lambda x: gene_dict.get(x['GENE']) + 0.4 + x["CHROM"] - (missing_chr.get(x['CHROM']) + 1), axis=1)
    else:
        pt_muts['y'] = 0
    
    # =============================================================================
    # Plot copy number and mutations
    # =============================================================================
    
    ax4.bar(pt_cn['x'], 0.8, bottom = pt_cn['y'], color = pt_cn['color'], zorder = 10)
    plot_mutations(pt_muts, ax4)
    ylables = get_ytick_lables(pt_cn.append(pt_ltc_cn))
    ax4.set_yticks(np.arange(0.4,len(ylables),1))
    ax4.set_xticks(np.arange(len(samples)))
    ax4.set_yticklabels(ylables, fontsize = 6)
    
    ax4.set_xlim(-0.5,len(samples)-0.5)
    ax4.set_ylim(0,len(gene_list)+len(pt_cn.append(pt_ltc_cn)['CHROM'].unique()))
    ax4.invert_yaxis()
    
    ax4.tick_params(left = False, bottom = False, labelbottom = False)
    ax4.spines['left'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    
    # =============================================================================
    #     Plot black bars
    # =============================================================================
    
    chr_bars = pt_cn.groupby('CHROM').max()[['y']]
    chr_bars.columns = ['max_y']
    chr_bars = chr_bars.merge(pt_cn.groupby('CHROM').min()[['y']], left_index = True, right_index = True)
    chr_bars.columns = ['max_y','min_y']
    chr_bars['y_ticks'] = (chr_bars['max_y'] + chr_bars['min_y']) / 2
        
    ax0.bar(0,chr_bars['max_y']-chr_bars['min_y']+0.8, bottom = chr_bars['min_y'], color = 'k')
    
    chr_labels = pt_cn.append(pt_ltc_cn)['CHROM'].astype(str).apply(lambda x: (x if x != '22' else 'X')).unique().tolist()
    ax0.set_yticks(chr_bars['y_ticks']+0.4)
    ax0.set_yticklabels(chr_labels, fontsize = 6, rotation = 90, va = 'center')
    ax0.set_ylim(0,len(gene_list)+len(pt_cn.append(pt_ltc_cn)['CHROM'].unique()))
    ax0.invert_yaxis()
    ax0.tick_params(axis = 'y', pad = 35, left = False)
    
    ax0.xaxis.set_visible(False)
    ax0.yaxis.set_visible(True)
    ax0.spines['left'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    
    # =============================================================================
    #     Plot unplotted samples on right
    # =============================================================================
    
    ltc_order_dict = dict(zip(pt_lowtc_samples, np.arange(len(pt_lowtc_samples))))
    
    pt_ltc_cn['x'] = pt_ltc_cn.apply(lambda x: ltc_order_dict.get(x['Sample ID']), axis = 1)
    if len(pt_muts)>0:
        pt_ltc_muts['x']=pt_ltc_muts.apply(lambda x: ltc_order_dict.get(x['Sample ID']), axis = 1)
    else:
        pt_ltc_muts['x'] = 0
    
    gene_list = pt_cn.append(pt_ltc_cn).sort_values(['CHROM','START'])['GENE'].unique().tolist()
    gene_dict = dict(zip(gene_list, np.arange(len(gene_list))))
    
    pt_ltc_cn['y'] = pt_ltc_cn.apply(lambda x: gene_dict.get(x['GENE']) + x["CHROM"]-(missing_chr.get(x['CHROM']) + 1), axis=1)
    if len(pt_ltc_muts) > 0:
        pt_ltc_muts['y'] = pt_ltc_muts.apply(lambda x: gene_dict.get(x['GENE']) + 0.4 + x["CHROM"] - (missing_chr.get(x['CHROM']) + 1), axis=1)
    else:
        pt_ltc_muts['y'] = 0
    
    # =============================================================================
    # Plot copy number and mutations
    # =============================================================================
    
    ax6.bar(pt_ltc_cn['x'], 0.8, bottom = pt_ltc_cn['y'], color = pt_ltc_cn['color'], zorder = 10)
    plot_mutations(pt_ltc_muts, ax6)
    ylables = get_ytick_lables(cn)
    ax6.set_yticks(np.arange(0.4,len(ylables),1))
    ax6.set_xticks(np.arange(len(samples)))
    ax6.set_yticklabels(ylables, fontsize = 6)
    
    ax6.set_xlim(-0.5,len(pt_lowtc_samples)-0.5)
    ax6.set_ylim(0,len(gene_list)+len(pt_cn.append(pt_ltc_cn)['CHROM'].unique()))
    ax6.invert_yaxis()
    
    ax6.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
    ax6.spines['left'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    
    # =============================================================================
    #     Label clustering unavailable
    # =============================================================================
       
    ax7.barh(0.66, len(pt_lowtc_samples)-0.2, left = -0.4, color = 'k', height = 0.67)
    
    ax7.set_xlim(-0.5,len(pt_lowtc_samples)-0.5)
    ax7.set_yticks([0.66])
    ax7.set_ylim(0,1)
    ax7.set_title('Low tumor content\nprevents accurate clustering', fontsize = 6)
    ax7.spines['left'].set_visible(False)
    ax7.spines['bottom'].set_visible(False)
    ax7.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)
    
    # =============================================================================
    #     Plot tumor content for low tc samples
    # =============================================================================
    
    ax8.bar(pts_ltc['Sample ID'], pts_ltc['Final tNGS_TC'], color = pts_ltc['tc_color'])
    ax8.plot([-0.5,len(pt_lowtc_samples)+0.5],[tc_cutoff,tc_cutoff], color = 'k', lw = 0.5, ls = 'dashed')
    
    ax8.set_xlim(-0.5,len(pt_lowtc_samples)-0.5)
    ax8.set_ylim(0,1)
    
    ax8.tick_params(bottom = False, labelbottom = True, left = False, labelleft = False)
    ax8.set_xticks(np.arange(0,len(pt_lowtc_samples),1))
    ax8.set_xticklabels(pt_lowtc_samples, rotation = 90)
    ax8.set_xticklabels(pd.Series(pt_lowtc_samples).str.split('_').str[2:].str.join('_'), rotation  = 90, fontsize = 6)
    
    # =============================================================================
    # Save figure    
    # =============================================================================
    
    plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/oncoprints/%i/%s_%s_mutsThenCN.png' % (int(tc_cutoff * 100), cohort, pt),bbox_inches='tight', dpi = 500)