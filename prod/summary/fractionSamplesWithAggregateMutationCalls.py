# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 10:14:58 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]


# =============================================================================
# 
# =============================================================================

fig, axs = plt.subplots(ncols = 2, figsize = (8,4), sharey = True, gridspec_kw={'width_ratios':[44,113]})

# =============================================================================
# 
# =============================================================================

for comb_type, ax in zip(['pb_comb','rp_comb'],axs):
    comb_tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/tumor_fraction/tumor_content.tsv', sep = '\t')
    comb_muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/mutations/melted.tsv', sep = '\t')
    comb_cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/copy_num/gene_cna.ffpe.tsv', sep = '\t')
    comb_samples = [s for s in comb_tc['Sample ID'].unique().tolist() if comb_type in s]
    
    comb_tc = comb_tc[comb_tc['Sample ID'].isin(comb_samples)]
    comb_muts = comb_muts[comb_muts['Sample ID'].isin(comb_samples)]
    comb_cn = comb_cn[comb_cn['Sample ID'].isin(comb_samples)]
    
    pts = comb_cn['Patient ID'].unique().tolist()
    
    # =============================================================================
    # 
    # =============================================================================
    
    muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations.tsv', sep = '\t')
    tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
    
    tc.columns = tc.iloc[0]
    tc = tc.drop(tc.index[0])
    tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
    
    tc = tc[~tc['Sample ID'].str.contains('cfDNA|MLN|MB|MT')]
    tc = tc[tc['Patient ID'].isin(pts)]
    
    muts = keepCodingMutations(muts)
    comb_muts = keepCodingMutations(comb_muts)
    
    muts = muts[muts['GENE'].isin(comb_cn['GENE'].unique().tolist())]
    comb_muts = comb_muts[comb_muts['GENE'].isin(comb_cn['GENE'].unique().tolist())]
    
    # =============================================================================
    # Keep only samples > 15% TC and mutations in those samples only
    # =============================================================================
    
    tc = tc[tc['Final tNGS_TC'] >= 0.15]
    muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
    
    # =============================================================================
    # Get index of all mutations
    # =============================================================================
    
    muts['mut_ID'] = muts.apply(lambda row: '%s_%s_%i_%s_%s_%s_%s' %\
                               (row['Patient ID'],row['CHROM'],row['POSITION'],row['REF'],row['ALT'],row['GENE'],row['EFFECT']), axis = 1)
        
    comb_muts['mut_ID'] = comb_muts.apply(lambda row: '%s_%s_%i_%s_%s_%s_%s' %\
                               (row['Patient ID'],row['CHROM'],row['POSITION'],row['REF'],row['ALT'],row['GENE'],row['EFFECT']), axis = 1)
    
    # =============================================================================
    # Fill in True if mutation called in combined samples    
    # =============================================================================
    
    comb_muts['Present_in_comb'] = True
    all_muts = pd.DataFrame({'mut_ID':muts['mut_ID']})
    
    all_muts = all_muts.merge(comb_muts[['mut_ID','Present_in_comb']], on = 'mut_ID', how = 'left')
    all_muts['Present_in_comb'] = all_muts['Present_in_comb'].fillna(False)
    
    # =============================================================================
    # 
    # =============================================================================
    
    pt_s_counts = tc.groupby('Patient ID').count()
    pt_s_counts = pt_s_counts [['Sample ID']]
    pt_s_counts.columns = ['Sample count']
    
    mut_count = muts.groupby(['mut_ID','Patient ID']).count()
    mut_count = mut_count[['Sample ID']]
    mut_count.columns = ['Mutation count']
    
    # =============================================================================
    # 
    # =============================================================================
    
    mut_count = mut_count.reset_index()
    pt_s_counts = pt_s_counts.reset_index()
    
    mut_count = mut_count.merge(pt_s_counts, on = 'Patient ID', how = 'left')
    
    mut_count['Fraction present'] = mut_count['Mutation count'] / mut_count['Sample count']
    
    all_muts = all_muts.merge(mut_count[['mut_ID','Mutation count','Sample count','Fraction present']], on = 'mut_ID', how = 'left')
    
    # =============================================================================
    # Sort all muts and prep for plotting
    # =============================================================================
    
    all_muts = all_muts.sort_values(['Present_in_comb','Fraction present'], ascending = [False, False])
    all_muts['Present_in_comb'] = all_muts['Present_in_comb'].map({True:'k',False:'grey'})
    all_muts = all_muts[all_muts['Sample count'] > 1]
    
    all_muts = all_muts.drop_duplicates()
    all_muts['x'] = np.arange(len(all_muts))
    
    # =============================================================================
    # 
    # =============================================================================
    
    ax.bar(all_muts['x'],all_muts['Fraction present'], color = all_muts['Present_in_comb'])
    print(len(all_muts))
    
    ax.set_xlim(-1,len(all_muts))
    ax.set_xticks([all_muts[all_muts['Present_in_comb'] == 'k']['x'].mean(),all_muts[all_muts['Present_in_comb'] == 'grey']['x'].mean()])
    ax.set_xticklabels(['Called\nin combined','Not detected\nin combined'])
    ax.set_title('%s\nn pt=%i' %(comb_type, len(pts)))
    
axs[0].set_ylabel('Primary samples with mutation present (%)')
fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/combinedSamplePercentConcordance.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/combinedSamplePercentConcordance.png')
