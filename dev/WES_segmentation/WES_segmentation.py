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


# =============================================================================
# Constants
# =============================================================================

color_dict = {0:'blue',1:'yellow',2:'purple',3:'red',4:'green'}

bin_size = 1000
sample = 'Ghent_ID5_RP1'

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

def chrom_to_bins(chr_all, bs):
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
            tmp = np.append(tmp,row['log_ratio'])
        else:
            lr_bins = np.append(lr_bins, np.median(tmp))
            mean_pos = np.append(mean_pos, (start+end)//2)
            mean_rel_pos = np.append(mean_rel_pos, (rel_start+rel_end) / 2)
            start = row['START']
            end = start + bs
            rel_start = start / chr_len
            rel_end = end / chr_len
            tmp = np.array(row['log_ratio'])
    return pd.DataFrame({'mean_pos':mean_pos, 'mean_rel_pos':mean_rel_pos, 'median_lr': lr_bins})

def get_cluster_dataframe(chrom):
    chrom = chrom[chrom['cluster'] != -1]
    chrom_gb = chrom.groupby('cluster')
    starts = chrom_gb.min()['mean_rel_pos']
    ends = chrom_gb.max()['mean_rel_pos']
    median_lrs = chrom_gb.median()['median_lr']
    return pd.DataFrame({'start':starts,'end':ends,'median_lr':median_lrs});


def check_clusters(chrom, clusters, min_samples):
    for index, cluster in clusters.iterrows():
        chrom = chrom[chrom['cluster'] == index].sort_values('mean_rel_pos')
        chrom.index = np.arange(0,len(chrom),1)
        min_p = 2.
        break_pt = -1
        n = 0
        for i in range(len(chrom)):
            if i < 15 or len(chrom) - i < 15: continue
            n += 1
            group1 = chrom[chrom.index <= i]
            group2 = chrom[chrom.index > i]
            stat, p = stats.ttest_ind(group1['median_lr'], group2['median_lr'])
            if p < min_p:
                min_p = p
                break_pt = chrom.at[i, 'mean_rel_pos']
        if n > 0 and min_p < 0.05 / (len(chrom) / n):
            print(min_p, break_pt)
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

# =============================================================================
# Change chromosomes to numeric data
# =============================================================================

lrs = lrs.replace({'chrX':'chr23', 'chrY':'chr24'})
lrs['CHROM'] = lrs['CHROM'].str.extract('(\d+)').astype(int)

#change from bases to kb
lrs['START'] = lrs['START'] / 1000
lrs['END'] = lrs['END'] / 1000

# Add column for relative position

lrs = normalize_chrom_positions(lrs)

# =============================================================================
# Create chromosome dictionary on binned chromosome
# =============================================================================

chr_dict = dict(zip(range(1,25), Parallel(n_jobs = 4)(delayed(chrom_to_bins)(lrs[lrs['CHROM'] == i].copy().reset_index(), bin_size) for i in range(1,25))))

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
# Check clusters
# =============================================================================
'''

'''



# =============================================================================
# Plot the dots
# =============================================================================

#Scale chromosome sizes to y chromosome
sizes = [12.75, 12.39, 10.13, 9.73, 9.27, 8.72, 8.14, 7.42, 7.08, 6.84, 6.87, 6.82, 4.89, 4.47, 4.17, 4.61, 4.25, 4.10, 2.99, 3.29, 1.93, 1.75, 7.97, 1]


fig,axs = plt.subplots(ncols = len(sizes), figsize = (15, 5), sharey = True, gridspec_kw = {'width_ratios':sizes})

for i,ax in enumerate(axs):
    chrom = chr_dict.get(i+1)
    clusters = cluster_dict.get(i+1)
    ax.scatter(chrom['mean_rel_pos'], chrom['median_lr'], s = 1, color = 'k')
    for index, row in clusters.iterrows():
        ax.plot([row['start'],row['end']],[row['median_lr'],row['median_lr']], color = color_dict.get(index), linewidth = 1.4)
    c_str = str(i+1)
    ax.set_xlabel(c_str)
    ax.tick_params(bottom = False, labelbottom = False)
    if i != 0:
        ax.spines['left'].set_linestyle((0,(4,4)))
        ax.tick_params(left = False)
    ax.set_ylim(-2,2)
    ax.set_yticks(np.arange(-2, 2.5, 0.5))

fig.tight_layout()    
fig.subplots_adjust(wspace = 0)
