# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:04:12 2021

@author: amurtha
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.legend as legend
from datetime import datetime as dt
import natsort as ns

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from os import walk

plt.rcParams.update({'font.size': 7.5})

sampleid_dict = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/wes_validation/wes_genesegment_names_targeted_id_dict.tsv', sep = '\t')
sampleid_dict = sampleid_dict.set_index('wes_id', drop=False)

#Get rid of samples not in sampleid_dict, ie not WES sequenced yet
target = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
target = target[target['Sample ID'].isin(sampleid_dict['targeted_id'].tolist())]

# Get all the segmentation file names
mypath = 'Y:/eritch/projects/ghent_m1rp/wxs/sequenza_segments_with_genes2'

f = []
for (dirpath, dirnames, filenames) in walk(mypath):
    f.extend(filenames)

del filenames, dirnames, mypath, dirpath

# For samples with duplicate files, get rid of old ones, use new ones
old = ['M1RP_ID13_C13_MLN4_WES_segments_withgene.txt',
       'M1RP_ID1_MLN4_segments_withgene.txt',
       'M1RP_ID1_RP1_segments_withgene.txt',
       'M1RP_ID1_RP2_segments_withgene.txt',
       'M1RP_ID7_MLN3_segments_withgene.txt',
       'M1RP_ID7_RP1_segments_withgene.txt']

f = [x for x in f if x not in old]
f = [x for x in f if 'wbc' not in x]

del old
# =============================================================================
# Format
# =============================================================================
inf_to_zero = True

targets = target['GENE'].unique().tolist()

target = target[['Sample ID', 'GENE', 'Copy_num', 'Log_ratio']]

target['cn_wes'] = 0
target['log_ratio_wes'] = 0.000


for sample in f:
    segfile = pd.read_csv('Y:/eritch/projects/ghent_m1rp/wxs/sequenza_segments_with_genes2/' + sample, sep = '\t')
    segfile = segfile[['depth.ratio', 'CNt', 'genes']]
    segfile['CNt'] = segfile['CNt'] - 2
    segfile['genes'] = segfile['genes'].str.split(';')

    for i, row in segfile.iterrows():
        all_genes = row['genes']
        if all_genes == all_genes:
            tar_genes = [x for x in all_genes if x in targets]
        else:
            tar_genes = []
        segfile.at[i, 'genes'] = tar_genes

    for i, row in segfile.iterrows():
        all_genes = row['genes']
        cn = row['CNt']
        lr = row['depth.ratio']
        for gene in all_genes:
            tar_id = sampleid_dict.at[sample, 'targeted_id']
            index = target[(target['Sample ID'] == tar_id) & (target['GENE'] == gene)].index[0]
            target.at[index, 'cn_wes'] = cn
            target.at[index, 'log_ratio_wes'] = lr

sampleid_dict.columns = ['wes_file_id', 'Sample ID']
sampleid_dict = sampleid_dict.reset_index(drop=True)
target = sampleid_dict.merge(target, how='left', on='Sample ID')

#Take log2 of depth ratios
target['log_ratio_wes'] = target['log_ratio_wes'].apply(lambda x: np.log2(x))

if inf_to_zero:
    # Set -inf or inf in log_ratio_wes to 0
    target = target.replace([np.inf, -np.inf], 0)
else:
    # OR delete rows with inf
    target = target.replace([np.inf, -np.inf], np.nan).dropna(subset=["log_ratio_wes"], how="all")



# =============================================================================
# Plot - black and white
# ============================================================================
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot()

ax.scatter(target['Log_ratio'], y=target['log_ratio_wes'] , marker='o',
        s=0.08, color='k', zorder = 100)

ax.set_xlabel('Targeted Log Ratio\nPer Gene')
ax.set_ylabel('WES Log Ratio\nVia Segmentation')

ax.set_ylim(-2,2)
ax.set_xlim(-2,2)
ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
ax.set_yticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

x = np.linspace(*ax.get_xlim())
ax.plot(x, x, color='red', lw=1, ls='--',dashes=(7.5, 5), zorder=1)

ax.axvline(x=0, color='k', lw=0.5)
ax.axhline(y=0, color='k', lw=0.5)
ax.grid(color='#bdbdbd',ls='--', lw=0.5, alpha=0.3, dashes=(5, 4), zorder=0)


# Statistical tests
slope, intercept, r_value, p_value, std_err = stats.linregress(target['Log_ratio'], target['log_ratio_wes'])

textstr = "r = " + str(r_value)[:4] + '\np = ' + str(p_value)[:4]
props = dict(boxstyle='square', facecolor='white', alpha=0)
ax.text(0.10, 0.86, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)


fig.savefig('C:\\Users\\Sarah\\Desktop\\M1_Figures\\log_ratio_wes_targeted_keep_weslr_0_black.pdf', bbox_extra_artists=(), bbox_inches='tight')



# =============================================================================
# Plot per gene - black and white
# ============================================================================
# WES vs targeted log ratio for each gene
genes = ['TP53','PTEN','RB1','BRCA2','ATM', 'AKT1', 'MYC']

tp53 = target[target['GENE'] == 'TP53']
pten = target[target['GENE'] == 'PTEN']
rb1 = target[target['GENE'] == 'RB1']
brca2 = target[target['GENE'] == 'BRCA2']
atm = target[target['GENE'] == 'ATM']
akt1 = target[target['GENE'] == 'AKT1']
myc = target[target['GENE'] == 'MYC']


# Set up plots and subplots for each gene
fig = plt.figure(figsize = (33,4))
gs = gridspec.GridSpec(ncols=7, nrows = 1)
gs.update(right = 0.97, top = 0.97, bottom = 0.25, left = 0.17)

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[0,2])
ax4 = fig.add_subplot(gs[0,3])
ax5 = fig.add_subplot(gs[0,4])
ax6 = fig.add_subplot(gs[0,5])
ax7 = fig.add_subplot(gs[0,6])


# To be used in loop
ax_to_gene = {ax1:'TP53', ax2:'PTEN', ax3:'RB1', ax4:'BRCA2', ax5:'ATM', ax6:'AKT1', ax7:'MYC'}
ax_to_df = {ax1:tp53, ax2:pten, ax3:rb1, ax4:brca2, ax5:atm, ax6:akt1, ax7:myc}

# Make plots
for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6, ax7]):

    ax.scatter(ax_to_df.get(ax)['Log_ratio'], y=ax_to_df.get(ax)['log_ratio_wes'] , marker='o',
        s=0.8, color='k', zorder = 100)

    ax.set_xlabel('Targeted Log Ratio\nPer Gene')
    ax.set_ylabel('WES Log Ratio\nVia Segmentation')
    ax.set_title((ax_to_gene.get(ax)))


    ax.set_ylim(-2,2)
    ax.set_xlim(-2,2)
    ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
    ax.set_yticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, color='red', lw=1, ls='--',dashes=(7.5, 5), zorder=1)

    ax.axvline(x=0, color='k', lw=0.5)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.grid(color='#bdbdbd',ls='--', lw=0.5, alpha=0.3, dashes=(5, 4), zorder=0)

    # Statistical tests
    slope, intercept, r_value, p_value, std_err = stats.linregress(ax_to_df.get(ax)['Log_ratio'], ax_to_df.get(ax)['log_ratio_wes'])

    textstr = "r = " + str(r_value)[:4] + '\np = ' + str(p_value)[:4] + str(p_value)[-4:]
    props = dict(boxstyle='square', facecolor='white', alpha=0)
    ax.text(0.10, 0.86, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)


fig.savefig('C:\\Users\\Sarah\\Desktop\\M1_Figures\\log_ratio_wes_vs_targeted_per_gene.pdf', bbox_extra_artists=(), bbox_inches='tight')



# =============================================================================
# Plot all with colour, separating deletions and gains
# =============================================================================

#final['status'] = ''
#
#for i, row in final.iterrows():
#    tar_cn = row['Copy_num']
#    wes_cn = row['cn_wes']
#    if tar_cn<0 and wes_cn<0:
#        final.at[i, 'status'] = 'both_deletion'
#    elif tar_cn<0 or wes_cn<0:
#        final.at[i, 'status'] = 'one_deletion'
#    elif tar_cn==0 and wes_cn==0:
#        final.at[i, 'status'] = 'neutral'
#    elif tar_cn>0 and wes_cn>0:
#        final.at[i, 'status'] = 'both_gain'
#    elif tar_cn>0 or wes_cn>0:
#        final.at[i, 'status'] = 'one_gain'
#
#nan_for_wes = final[final['status'] == '']
#
##Scatterplot coloured WES vs targeted separated by colour, if mutation detected in one or both seq
#cat = ['both_deletion', 'one_deletion', 'both_gain', 'one_gain', 'neutral']
#cat_map = {'both_deletion':'#000ceb', 'one_deletion':'#7077ff', 'both_gain':'red', 'one_gain':'#ff8585', 'neutral':'#cfcfcf'}
#zorder_map = {'both_deletion':2, 'one_deletion':1, 'both_gain':2, 'one_gain':1, 'neutral':10}
#
#
#fig = plt.figure(figsize=(4,4))
#ax = fig.add_subplot()
#
#groups = final.groupby('status')
#for name, group in groups:
#    group.plot(ax=ax, kind='scatter', x ='Log_ratio', y='log_ratio_wes', marker='o',
#               s=0.0001, label=name, color = cat_map[name], zorder=zorder_map[name])
#ax.set_xlabel('Targeted Log Ratio\nPer Gene')
#ax.set_ylabel('WES Log Ratio\nVia Segmentation')
#
#
#ax.set_ylim(-2,2)
#ax.set_xlim(-2,2)
#ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
#ax.set_yticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
#
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#
#x = np.linspace(*ax.get_xlim())
#ax.plot(x, x, color='red', lw=1, ls='--',dashes=(7.5, 5), zorder=1)
#
#ax.axvline(x=0, color='k', lw=0.5)
#ax.axhline(y=0, color='k', lw=0.5)
#ax.grid(color='#bdbdbd',ls='--', lw=0.5, alpha=0.3, dashes=(5, 4), zorder=0)
#
#
## Statistical tests
#slope, intercept, r_value, p_value, std_err = stats.linregress(final['Log_ratio'], final['log_ratio_wes'])
#
#textstr = "Pearson\nr = " + str(r_value)[:4]
#props = dict(boxstyle='square', facecolor='white', alpha=0)
#ax.text(0.10, 0.86, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#
##Legend
#labels = ['Gain in both', 'Gain in one', 'Loss in both', 'Loss in one', 'Copy neutral']
#handles = [mlines.Line2D([], [], color='red', markeredgecolor='red', marker='o', lw=0, markersize=7),
#           mlines.Line2D([], [], color='#ff8585', markeredgecolor='#ff8585', marker='o', lw=0, markersize=7),
#           mlines.Line2D([], [], color='#000ceb', markeredgecolor='#000ceb', marker='o', lw=0, markersize=7),
#           mlines.Line2D([], [], color='#7077ff', markeredgecolor='#7077ff', marker='o', lw=0, markersize=7),
#           mlines.Line2D([], [], color='#cfcfcf', markeredgecolor='#cfcfcf', marker='o', lw=0, markersize=7)]
#ax.legend(handles, labels, loc = 'lower right', prop={'size': 7}, handlelength = 1.6, frameon = False)
#
#
#fig.savefig('C:\\Users\\Sarah\\Desktop\\log_ratio_wes_targeted_keep_weslr_0.pdf', bbox_extra_artists=(), bbox_inches='tight')
#
#
## See the calls with neg log ratio but still called as gain
#test = final[final['status'] == 'one_gain']
#test = test[test['log_ratio_wes'] < 0]
#test = test[test['Log_ratio'] < 0]
#
#







