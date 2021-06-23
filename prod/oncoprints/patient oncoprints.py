# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:56:24 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]


gene_list = ['TP53','RB1','PTEN','BRCA2','CDK12','SPOP','AR']

sample_type_dict = {'PB':0,'RP':1,'MLN':2,'MB':2,'MT':2,'cfDNA':3}
cn_color_dict = {-2:'#3F60AC',-1:'#9CC5E9', 0:'#E6E7E8', 1:'#F59496', 2:'#EE2D24'}
mut_effect_dict = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907'}

# =============================================================================
# Import TC, mutation, and CN data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

muts = keepCodingMutations(muts)
muts = muts[muts['GENE'].isin(cn['GENE'].unique().tolist())]

gl['NOTES'] = gl['NOTES'].fillna('')
gl = gl[gl['NOTES'].str.contains('ClinVar:Pathogenic.')]

# =============================================================================
# Loop over patients and create an oncoprint
# =============================================================================

pts = ['ID%i' % (i+1) for i in np.arange(43)]

for pt in pts:
    pt_tc = tc[tc['Patient ID'] == pt].copy()
    pt_muts = muts[muts['Patient ID'] == pt].copy()
    pt_cn = cn[cn['Patient ID'] == pt].copy()
    pt_gl = gl[gl['Patient ID'] == pt]
    
    # =============================================================================
    # Label y coordinates    
    # =============================================================================
    
    pt_genes = gene_list + [gene for gene in pt_muts['GENE'].unique().tolist() if gene not in gene_list]
    pt_muts['y'] = pt_muts['GENE'].map(dict(zip(pt_genes, np.arange(len(pt_genes)))))
    pt_cn['y'] = pt_cn['GENE'].map(dict(zip(pt_genes, np.arange(len(pt_genes)))))
    
    pt_cn = pt_cn[pt_cn['GENE'].isin(pt_genes)]
    pt_muts = pt_muts[pt_muts['GENE'].isin(pt_genes)]
    pt_gl = pt_gl[pt_gl['GENE'].isin(pt_genes)]
    
    # =============================================================================
    # Add pathogenic germline variant to the mutation file    
    # =============================================================================
    
    # pt_samples = pt_tc['Sample ID'].tolist()
    # for index, row in pt_gl.iterrows():
    #     tmp_df = pd.DataFrame({'Sample ID': pt_samples, }))
    
    # =============================================================================
    # Add mutation color    
    # =============================================================================
    
    pt_muts['Color'] = pt_muts['EFFECT'].str.split(' ').str[0].map(mut_effect_dict)
    
    # =============================================================================
    # Label x coordinates     
    # =============================================================================
    
    pt_tc['Sample type'] = pt_tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
    pt_tc['Sample type'] = pt_tc['Sample type'].map(sample_type_dict)
    
    pt_tc = pt_tc.sort_values(['Sample type','Final tNGS_TC'], ascending = [True, False])
    sample_order_dict = dict(zip(pt_tc['Sample ID'],np.arange(len(pt_tc))))
    pt_tc['x'] = pt_tc['Sample ID'].map(sample_order_dict)
    pt_muts['x'] = pt_muts['Sample ID'].map(sample_order_dict)
    pt_cn['x'] = pt_cn['Sample ID'].map(sample_order_dict)
    
    # =============================================================================
    # Create the figure
    # =============================================================================
    
    fig,[ax1,ax2] = plt.subplots(nrows=2, sharex = True, gridspec_kw={'height_ratios':{1,2}}, figsize = (2.5,2))
    
    # =============================================================================
    # Plot tumor content    
    # =============================================================================
    
    ax1.bar(pt_tc['x'],pt_tc['Final tNGS_TC'])
    ax1.set_ylim(0,1)
    ax1.set_yticks(np.arange(0,1.01,0.25))
    ax1.set_ylabel('Tumor\ncontent', fontsize = 6)
    ax1.tick_params(labelsize = 6, bottom = False)
    
    ax1.set_xlim(-0.55, len(sample_order_dict)-0.4)
    
    # =============================================================================
    # Plot oncoprint    
    # =============================================================================
    
    ax2.bar(pt_cn['x'],0.8,bottom = pt_cn['y'], color = pt_cn['Copy_num'].map(cn_color_dict), zorder = 10)
    ax2.scatter(pt_muts['x'],pt_muts['y']+0.4, s = 6, lw = 0, c = pt_muts['Color'], zorder = 100, marker = 's')
    
    # =============================================================================
    # Labels    
    # =============================================================================
    
    ax2.set_yticks(np.arange(0.4,len(pt_genes),1))
    ax2.set_yticklabels(pt_genes, fontsize = 6)
    
    ax2.set_xticks(np.arange(0,len(sample_order_dict),1))
    ax2.set_xticklabels(pt_tc['Sample ID'].str.split('_').str[2], rotation = 90, fontsize = 6)
    
    ax2.tick_params(labelsize = 6, bottom = False, left = False, pad = 0)
    
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    plt.tight_layout()
    fig.subplots_adjust(hspace = 0.05)
    
    
    
    fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/Patient oncoprints/%s.pdf' % pt)