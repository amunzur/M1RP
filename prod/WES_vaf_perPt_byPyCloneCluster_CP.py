# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 11:46:16 2021

@author: amurtha
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

cmap=plt.get_cmap('tab20')

# =============================================================================
# WES mutations
# =============================================================================

# for pt in ['ID%i' % (i+1) for i in np.arange(43)]:
for pt in ['ID6']:
    muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/betastasis/june2021_wes_mutations_betastasis.tsv', sep = '\t')
    
    # =============================================================================
    # Keep columns that contain ID33
    # =============================================================================
    
    index_cols = ['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES']
    
    muts = muts[index_cols + [col for col in muts.columns.tolist() if '_%s_' % pt in col and not '_NL' in col and not '_WBC' in col and not '_gDNA' in col]]
    muts['EFFECT'] = muts['EFFECT'].str.split(' ').str[0]
    muts = muts[muts.apply(lambda r: r.str.contains('*', regex = False).any(), axis=1)]
    
    # =============================================================================
    # Get allele frequency
    # =============================================================================
    
    n_samples = len([col for col in muts.columns.tolist() if col not in index_cols])
    
    p_samples = [col for col in muts.columns.tolist() if col not in index_cols and not ('MLN' in col or 'MT' in col or 'cfDNA' in col or 'ML' in col)]
    m_samples = [col for col in muts.columns.tolist() if col not in index_cols and ('MLN' in col or 'MT' in col or 'ML' in col)]
    c_samples = [col for col in muts.columns.tolist() if col not in index_cols and 'cfDNA' in col]
    
    samples = c_samples+m_samples+p_samples
    
    # =============================================================================
    # Check if pyclone file exists. If so color by clone ID
    # =============================================================================
    fp = 'X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/loci_filtered.tsv' % pt
    if os.path.exists(fp):
        loci = pd.read_csv('X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/loci_filtered.tsv' % pt, sep = '\t')
        cluster = pd.read_csv('X:/users/ewarner/ghent_m1/pyclone_wes/%s/pyclone_output/tables/cluster.tsv' % pt, sep = '\t')
        cluster_vals = pd.Series(loci['cluster_id'].unique())
        n_clusters = len(cluster_vals)
        cluster2color = dict(zip(cluster_vals,\
                                 [matplotlib.colors.rgb2hex(cmap(i*cmap.N//n_clusters)) for i in cluster_vals]))
        
        # =============================================================================
        # Pivot loci file to show ccf and map colors onto copy of it        
        # =============================================================================
        
        ccf = loci.pivot(index = 'mutation_id',columns = 'sample_id', values = 'cellular_prevalence')
        ccf = ccf.merge(loci[['mutation_id','cluster_id']].drop_duplicates(), left_index = True, right_on = 'mutation_id')
        ccf['color'] = ccf['cluster_id'].map(cluster2color)
            
        # =============================================================================
        # Set up x values by sorting by cluster then ccf by sample        
        # =============================================================================
        
        ccf = ccf.sort_values(['cluster_id']+samples, ascending = [True]+([False]*len(samples)))
        ccf['x'] = np.arange(len(ccf))
        
        # =============================================================================
        # Create axis. plot bars via x and ccf. Get color from color dataframe
        # =============================================================================
        
        if pt == 'ID8':
            fig,axs = plt.subplots(nrows = len(samples), figsize = (20,20))
        else:
            fig,axs = plt.subplots(nrows = len(samples), figsize = (3.5,2.75))
        for s, ax in zip(samples, axs):
            ax.bar(ccf['x'], ccf[s], width = 1, color = ccf['color'])
            ax.tick_params(bottom = False, labelbottom = False, labelsize = 6)
            ax.set_ylabel('VAF', fontsize = 6)
            ax.set_xlim(-1,len(muts)+1)
            ax.set_xlabel(s, fontsize = 6)
        plt.tight_layout()
        plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Patient specific/WES VAF by cluster/CP/%s.CP.pdf' % pt)
        
'''        
    # =============================================================================
    #     
    # =============================================================================

    if pt == 'ID8':
        fig,axs = plt.subplots(nrows = len(samples), figsize = (20,20))
    else:
        fig,axs = plt.subplots(nrows = len(samples), figsize = (3.5,2.75))
    
    for s, ax in zip(samples, axs):
        if 'color' in muts.columns.tolist(): 
            ax.bar(muts['x'], muts[s+'_vaf'], width = 1, color = muts['color'])
        else:
            ax.bar(muts['x'], muts[s+'_vaf'], width = 1, color = muts[s+'_color'])
        ax.tick_params(bottom = False, labelbottom = False, labelsize = 6)
        ax.set_ylabel('VAF', fontsize = 6)
        ax.set_xlim(-1,len(muts)+1)
        ax.set_xlabel(s, fontsize = 6)
    
    # fig.suptitle(pt, fontsize=6)
    plt.tight_layout()
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Patient specific/WES VAF by cluster/%s.CCF.pdf' % pt)
    
'''