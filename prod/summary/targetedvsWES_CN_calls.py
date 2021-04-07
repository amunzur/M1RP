# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 10:30:56 2020

@author: Sarah
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from os import walk

plt.rcParams.update({'font.size': 7.5})

target = pd.read_csv('C:/Users/Sarah/Desktop/dropbox data m1rp Jan 28/M1RP_allSamples_cna.tsv', sep = '\t')

sampleid_dict = pd.read_csv('C:/Users/Sarah/Desktop/wes_genesegment_names_targeted_id_dict.tsv', sep = '\t')
sampleid_dict = sampleid_dict.set_index('wes_id', drop=False)

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
    final = target.replace([np.inf, -np.inf], 0)
else:
    # OR delete rows with inf
    final = target.replace([np.inf, -np.inf], np.nan).dropna(subset=["log_ratio_wes"], how="all")

        
    
# =============================================================================
# Plot - black and white
# ============================================================================
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot()

ax.scatter(final['Log_ratio'], y=final['log_ratio_wes'] , marker='o',
        s=0.001, color='k', zorder = 100)

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
slope, intercept, r_value, p_value, std_err = stats.linregress(final['Log_ratio'], final['log_ratio_wes'])

textstr = "Pearson\nr = " + str(r_value)[:4] 
props = dict(boxstyle='square', facecolor='white', alpha=0)
ax.text(0.10, 0.86, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    
fig.savefig('C:\\Users\\Sarah\\Desktop\\log_ratio_wes_targeted_keep_weslr_0_black.pdf', bbox_extra_artists=(), bbox_inches='tight')
    
    
    
    

# =============================================================================
# Plot with colour, separating deletions and gains
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














