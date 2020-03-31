# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:19:44 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import multiprocessing
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

# =============================================================================
# Constants
# =============================================================================

color_dict = {0:'blue',1:'green',2:'purple',3:'red',4:'orange'}

bin_size = 1000 ## Units: kb
# sample = sys.argv[1]
sample = 'M1RP_ID5_RP1'
if '.pdf' in sample:
	sample = sample.split('.pdf')[0]

# =============================================================================
# Helpers
# =============================================================================

def normalize_chrom_positions(lrs):
    #find length of chromosome
    max_pos = lrs.groupby('CHROM').max()[['END']].rename(columns = {'END':'CHR_LEN'})
    lrs = lrs.merge(max_pos, right_index = True, left_on = 'CHROM', how = 'left')
    #divide positions by length
    lrs['rel_pos'] = lrs['START'] / lrs['CHR_LEN']
    return lrs;

def chrom_to_bins(chr_all, bs, col, new_col):
    chr_len = chr_all['CHR_LEN'][0]
    start = chr_all['START'][0]
    end = start + bs
    rel_start = start / chr_len
    rel_end = end / chr_len
    lr_bins = np.array([])
    mean_pos = np.array([])
    mean_rel_pos = np.array([])
    tmp = np.array([])
    for index, row in chr_all.iterrows():
        if row['END'] < end:
            tmp = np.append(tmp,row[col])
        else:
            if tmp.size >= 5:
                lr_bins = np.append(lr_bins, np.median(tmp))
                mean_pos = np.append(mean_pos, (start+end)//2)
                mean_rel_pos = np.append(mean_rel_pos, (rel_start+rel_end) / 2)
            start = row['START']
            end = start + bs
            rel_start = start / chr_len
            rel_end = end / chr_len
            tmp = np.array(row[col])
    return pd.DataFrame({'mean_pos':mean_pos, 'mean_rel_pos':mean_rel_pos, new_col: lr_bins})

def get_cluster_dataframe(chrom):
    chrom = chrom[chrom['cluster'] != -1]
    chrom_gb = chrom.groupby('cluster')
    starts = chrom_gb.min()['mean_rel_pos']
    ends = chrom_gb.max()['mean_rel_pos']
    median_lrs = chrom_gb.median()['median_lr']
    return pd.DataFrame({'start':starts,'end':ends,'median_lr':median_lrs});

def update_clusters(chrom, clusters, cluster_num, break_pt):
    chrom.loc[chrom['cluster'] > cluster_num, 'cluster'] = chrom['cluster'] + 1
    chrom.loc[(chrom['cluster'] == cluster_num) & (chrom['mean_rel_pos'] > break_pt), 'cluster'] = chrom['cluster'] + 1
    clusters = get_cluster_dataframe(chrom)
    
    return clusters, chrom    

def check_clusters(chrom, clusters, min_samples):
    for cluster_num, cluster in clusters.iterrows():
        cluster_values = chrom[chrom['cluster'] == cluster_num].sort_values('mean_rel_pos')
        cluster_values.index = np.arange(0,len(cluster_values),1)
        min_p = 2.
        break_pt = -1
        n = 0
        for i in range(len(cluster_values)):
            if i < 25 or len(cluster_values) - i < 25: continue
            n += 1
            group1 = cluster_values[cluster_values.index <= i]
            group2 = cluster_values[cluster_values.index > i]
            stat, p = stats.ttest_ind(group1['median_lr'], group2['median_lr'])
            if p < min_p:
                min_p = p
                break_pt = chrom.at[i, 'mean_rel_pos']
        if n > 0 and min_p < 0.05 / (n):
            # print(min_p, break_pt)
            (clusters, chrom) = update_clusters(chrom, clusters, cluster_num, break_pt)
    return clusters
        
# =============================================================================
# =============================================================================
# =============================================================================
    
''' 
MAIN

This script creates segmented WES data using the algorithm DB scan. WES
is seperated by chromosome, then binned. Epoch and n_sample parameteres are 
determined based on the number of bins per chromosome.
'''

# =============================================================================
# Import data
# =============================================================================

lrs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks/'+sample+'_logratio.igv', sep = '\t').rename(columns = {sample:'log_ratio'})

snps = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks/'+sample+'_hetz_snp.igv', sep = '\t').rename(columns = {sample:'hetz_snp'})

tumor_frac = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

# =============================================================================
# Change chromosomes to numeric data
# =============================================================================

lrs = lrs.replace({'chrX':'chr23', 'chrY':'chr24'})
lrs['CHROM'] = lrs['CHROM'].str.extract('(\d+)').astype(int)

snps = snps.replace({'chrX':'chr23', 'chrY':'chr24'})
snps['CHROM'] = snps['CHROM'].str.extract('(\d+)').astype(int)

#change from bases to kb
lrs['START'] = lrs['START'] / 1000
lrs['END'] = lrs['END'] / 1000

snps['START'] = snps['START'] / 1000
snps['END'] = snps['END'] / 1000

snps = snps[~snps['hetz_snp'].isin([0,1])]
snps['hetz_snp'] = (snps['hetz_snp'] - 0.5).abs()

# Add column for relative position

lrs = normalize_chrom_positions(lrs)

snps = snps.merge(lrs[['CHROM', 'CHR_LEN']].drop_duplicates(), on = 'CHROM', how = 'left')

# =============================================================================
# # Get mutation for TF check
# =============================================================================
tumor_frac = tumor_frac.replace({'chrX':'chr23', 'chrY':'chr24'})
tumor_frac['Chromosome'] = tumor_frac['Chromosome'].fillna('0').str.extract('(\d+)').astype(int)

tumor_frac = tumor_frac.loc[tumor_frac['Sample M1RP_ID'] == sample].squeeze()

tumor_frac['Position'] = tumor_frac['Position'] / 1000
tumor_frac['Rel_pos'] = tumor_frac['Position'] / lrs[lrs['CHROM'] == tumor_frac['Chromosome']]['CHR_LEN'].unique()[0]

# =============================================================================
# Create chromosome dictionary on binned chromosome
# =============================================================================

chr_dict = dict(zip(range(1,25), Parallel(n_jobs = 4)(delayed(chrom_to_bins)(lrs[lrs['CHROM'] == i].copy().reset_index(), bin_size, 'log_ratio', 'median_lr') for i in range(1,25))))

chr_snp_dict = dict(zip(range(1,25), Parallel(n_jobs = 4)(delayed(chrom_to_bins)(snps[snps['CHROM']==i].copy().reset_index(), bin_size, 'hetz_snp', 'median_vaf') for i in range(1,25))))


# =============================================================================
# Segement bins by chromosome
# =============================================================================

cluster_dict = {}

for i in range(1,25):
    chrom = chr_dict.get(i)
    epoch = float(len(chrom) / (1000 / (1 + (i+1) / 10)))
    min_samples = max(int(len(chrom) / 20), 3)
    # print(epoch, min_samples)
    outlier_detection = DBSCAN(eps = epoch, metric = 'euclidean', min_samples = min_samples, n_jobs = 4)
    clusters = outlier_detection.fit_predict(chrom[['mean_rel_pos','median_lr']])
    chrom['cluster'] = clusters
    # print(i, chrom['cluster'].unique())
    print(i)
    df_clusters = get_cluster_dataframe(chrom)
    df_clusters = check_clusters(chrom, df_clusters, min_samples)
    cluster_dict[i] = df_clusters

# =============================================================================
# Plot the dots
# =============================================================================

#Scale chromosome sizes to y chromosome
sizes = [12.75, 12.39, 10.13, 9.73, 9.27, 8.72, 8.14, 7.42, 7.08, 6.84, 6.87, 6.82, 4.89, 4.47, 4.17, 4.61, 4.25, 4.10, 2.99, 3.29, 1.93, 1.75, 7.97, 1]


fig,axs = plt.subplots(ncols = len(sizes), nrows = 2, figsize = (15, 9), sharey = 'row', sharex = 'col', gridspec_kw = {'width_ratios':sizes})

for i,ax in enumerate(axs[0]):
    chrom = chr_dict.get(i+1)
    clusters = cluster_dict.get(i+1)
    ax.scatter(chrom['mean_rel_pos'], chrom['median_lr'], s = 1, color = 'k', zorder = 100, marker = 'o')
    for index, row in clusters.iterrows():
        lw = 1.5
        if i+1 == tumor_frac['Chromosome'] and row['start'] < tumor_frac['Rel_pos'] and row['end'] > tumor_frac['Rel_pos']:
            lw = 2
        ax.plot([row['start'],row['end']],[row['median_lr'],row['median_lr']], color = color_dict.get(index), linewidth = lw, zorder = 1000)
    if i+1 == tumor_frac['Chromosome']:
        ax.plot([tumor_frac['Rel_pos'], tumor_frac['Rel_pos']], [-2,2])
    c_str = str(i+1)
    ax.set_xlabel(c_str)
    ax.tick_params(bottom = False, labelbottom = False)
    if i != 0:
        ax.spines['left'].set_linestyle((0,(4,4)))
        ax.tick_params(left = False)
    ax.set_ylim(-2,2)
    ax.set_yticks(np.arange(-2, 2.5, 0.5))
    ax.grid(axis = 'y', lw = 0.5, color = '0.65', linestyle = 'dashed', zorder = 0)
    
for i, ax in enumerate(axs[1]):
    chrom = chr_snp_dict.get(i+1)
    ax.scatter(chrom['mean_rel_pos'], chrom['median_vaf'], s = 1, color = 'k', zorder = 100, marker = 'o')
    c_str = str(i+1)
    ax.set_xlabel(c_str)
    ax.tick_params(bottom = False, labelbottom = False)
    if i != 0:
        ax.spines['left'].set_linestyle((0,(4,4)))
        ax.tick_params(left = False)
    ax.set_ylim(0,1)
    ax.grid(axis = 'y', lw = 0.5, color = '0.65', linestyle = 'dashed', zorder = 0)
    

fig.tight_layout()    
fig.subplots_adjust(wspace = 0)

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Figures/WES_segmentation/'+sample+'.pdf')