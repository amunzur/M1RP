# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:41:00 2021

@author: amurtha
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

cmap=plt.get_cmap('tab20')


# cluster2color = dict(zip(np.arange(cmap.N),[matplotlib.colors.rgb2hex(cmap(i)) for i in np.arange(cmap.N)]))

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] > 0]

tc_pos = tc['Sample ID'].tolist()
tc = tc.set_index('Sample ID')

rd = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/exome_sequencing_metrics_June2021.tsv', sep = '\t')
rd = rd.drop_duplicates('SAMPLE').set_index('SAMPLE')

# =============================================================================
# WES mutations
# =============================================================================

for pt in ['ID%i' % (i+1) for i in np.arange(43)]:
# for pt in ['ID40']:
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
    
    p_samples = [col for col in muts.columns.tolist() if col not in index_cols and not ('MLN' in col or 'MT' in col or 'cfDNA' in col or 'ML' in col)]
    m_samples = [col for col in muts.columns.tolist() if col not in index_cols and ('MLN' in col or 'MT' in col or 'ML' in col)]
    c_samples = [col for col in muts.columns.tolist() if col not in index_cols and 'cfDNA' in col]
    
    samples = c_samples+m_samples+p_samples
    samples = [s for s in samples if s in tc_pos]
    
    if len(samples) <= 1:
        print('1 sample')
        continue;
    
    for s in samples:
        muts[s+'_vaf'] = muts[s].str.split('%').str[0].astype(float) / 100
        muts[s+'_depth'] = muts[s].str.split('(').str[1].str.split(')').str[0].astype(int)
        muts[s+'_called'] = 0
        muts.loc[muts[s].str.contains('*', regex = False), s+'_called'] = 1
        muts[s+'_color'] = muts[s+'_called'].map({0:'grey',1:'black'})
        
    # muts = muts.sort_values([s+'_called' for s in samples]+[s+'_vaf' for s in samples], ascending = False)
    muts = muts.sort_values([s+'_vaf' for s in samples], ascending = False)
    muts['x'] = np.arange(len(muts))
    
    # =============================================================================
    # Check if pyclone file exists. If so color by clone ID
    # =============================================================================
    fp = 'X:/users/amurtha/Ghent/pyclone_wes/%s/pyclone_output2/tables/loci_filtered.tsv' % pt
    if os.path.exists(fp):
        loci = pd.read_csv('X:/users/amurtha/Ghent/pyclone_wes/%s/pyclone_output2/tables/loci_filtered.tsv' % pt, sep = '\t')
        cluster = pd.read_csv('X:/users/amurtha/Ghent/pyclone_wes/%s/pyclone_output2/tables/cluster.tsv' % pt, sep = '\t')
        muts['mutation_id'] = muts['CHROM']+'_'+muts['GENE']+'_'+muts['POSITION'].astype(str)+'_'+muts['REF']+'_'+muts['ALT']+'_'+muts['EFFECT']
        muts = muts.merge(loci[['mutation_id','cluster_id']].drop_duplicates(), on = 'mutation_id')
        n_clusters = muts['cluster_id'].max()+1
        cluster2color = dict(zip(np.arange(n_clusters),[matplotlib.colors.rgb2hex(cmap(i*cmap.N//n_clusters)) for i in np.arange(n_clusters)]))
        muts['color'] = muts['cluster_id'].map(cluster2color)
        
        cluster_hierarchy = cluster.copy()
        cluster_hierarchy = cluster_hierarchy[cluster_hierarchy['mean'] > 0.05]
        cluster_hierarchy = cluster_hierarchy.groupby('cluster_id').count()[['sample_id']]
        cluster_hierarchy.columns = ['truncal_score']
        muts = muts.merge(cluster_hierarchy, left_on = 'cluster_id', right_index = True, how = 'left')
        muts['truncal_score'] = muts['truncal_score'].fillna(0)
        
        
        muts = muts.sort_values(['truncal_score','cluster_id']+[s+'_vaf' for s in samples], ascending = [False, True]+[False]*len(samples))
        muts['x'] = np.arange(len(muts))
    else:
        print(pt)
        
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
        ax.set_xlabel(s + ', TC: %.2f, Median depth: %i' % (tc.at[s, 'Final tNGS_TC'], rd.at[s,'COVERAGE (MEDIAN)']), fontsize = 6)
    
    # fig.suptitle(pt, fontsize=6)
    plt.tight_layout()
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Patient specific/WES VAF by cluster/pyclone2/%s.pdf' % pt)