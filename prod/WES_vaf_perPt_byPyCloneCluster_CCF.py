# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 12:12:58 2021

@author: amurtha
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

cmap=plt.get_cmap('tab20')
# cluster2color = dict(zip(np.arange(cmap.N),[matplotlib.colors.rgb2hex(cmap(i)) for i in np.arange(cmap.N)]))

# =============================================================================
# WES mutations
# =============================================================================

for pt in ['ID%i' % (i+1) for i in np.arange(43)]:
# for pt in ['ID13']:
    if pt == 'ID13':
        continue;
    ccf = pd.read_csv('G:\Evan MSc\Ghent M1\Fall2020_Updated_Analysis\Clean_data\ccf_estimations_for_wes_mutations.tsv', sep = '\t')
    ccf = ccf[ccf['Patient ID'] == pt]
    # =============================================================================
    # Keep columns that contain ID33
    # =============================================================================
    
    index_cols = ['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES']
    
    ccf['mutation_id'] = ccf['CHROM']+'_'+ccf['GENE']+'_'+ccf['START'].astype(str)+'_'+ccf['REF']+'_'+ccf['ALT']+'_'+ccf['EFFECT']
    
    ccf = ccf.pivot(index = 'mutation_id', columns = 'Sample ID', values = 'adjusted_ccf')
    
    # =============================================================================
    # Get allele frequency
    # =============================================================================
    
    n_samples = len(ccf.columns.tolist())
    
    p_samples = [col for col in ccf.columns.tolist() if col not in index_cols and not ('MLN' in col or 'MT' in col or 'cfDNA' in col or 'ML' in col)]
    m_samples = [col for col in ccf.columns.tolist() if col not in index_cols and ('MLN' in col or 'MT' in col or 'ML' in col)]
    c_samples = [col for col in ccf.columns.tolist() if col not in index_cols and 'cfDNA' in col]
    
    samples = c_samples+m_samples+p_samples
    
    ccf = ccf.reset_index()
    
    if len(samples) <= 1: 
        continue;
        
    # =============================================================================
    # Check if pyclone file exists. If so color by clone ID
    # =============================================================================
    fp = 'X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/loci_filtered.tsv' % pt
    if os.path.exists(fp):
        loci = pd.read_csv('X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/loci_filtered.tsv' % pt, sep = '\t')
        cluster = pd.read_csv('X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/cluster.tsv' % pt, sep = '\t')
        ccf = ccf.merge(loci[['mutation_id','cluster_id']].drop_duplicates(), on = 'mutation_id')
        cluster_vals = pd.Series(ccf['cluster_id'].unique())
        n_clusters = len(cluster_vals)
        cluster2color = dict(zip(cluster_vals,\
                                 [matplotlib.colors.rgb2hex(cmap(i*cmap.N//n_clusters)) for i in cluster_vals]))
        ccf['color'] = ccf['cluster_id'].map(cluster2color)
        ccf = ccf.sort_values(['cluster_id']+samples, ascending = [True]+[False]*len(samples))
    else:
        print(pt)
        continue;
        
    # =============================================================================
    #     
    # =============================================================================

    ccf['x'] = np.arange(len(ccf))
    if pt == 'ID8':
        fig,axs = plt.subplots(nrows = len(samples), figsize = (20,20))
    else:
        fig,axs = plt.subplots(nrows = len(samples), figsize = (3.5,2.75))
    
    for s, ax in zip(samples, axs):
        if 'color' in ccf.columns.tolist(): 
            ax.bar(ccf['x'], ccf[s], width = 1, color = ccf['color'])
        else:
            ax.bar(ccf['x'], ccf[s], width = 1, color = ccf[s+'_color'])
        ax.tick_params(bottom = False, labelbottom = False, labelsize = 6)
        ax.set_ylabel('VAF', fontsize = 6)
        ax.set_xlim(-1,len(ccf)+1)
        ax.set_xlabel(s, fontsize = 6)
    
    # fig.suptitle(pt, fontsize=6)
    plt.tight_layout()
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Patient specific/WES VAF by cluster/CCF/%s.CCF.pdf' % pt)