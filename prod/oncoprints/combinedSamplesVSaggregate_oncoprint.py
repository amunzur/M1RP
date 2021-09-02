# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 14:09:50 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

additional_genes = ['AR','TP53','PTEN','RB1','CHD1','SPOP','FOXA1', 'MSH2','MSH6']

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

muts_color = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907'}



def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# 
# =============================================================================
for comb_type in ['rp_comb','pb_comb']:    
    comb_tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/tumor_fraction/tumor_content.tsv', sep = '\t')
    comb_muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/mutations/melted.tsv', sep = '\t')
    comb_cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/combined_primary/copy_num/gene_cna.ffpe.tsv', sep = '\t')
    
    comb_cn = comb_cn[comb_cn['GENE']!='FOXA1']
    comb_samples = [s for s in comb_tc['Sample ID'].unique().tolist() if comb_type in s]
    
    comb_tc = comb_tc[comb_tc['Sample ID'].isin(comb_samples)]
    comb_muts = comb_muts[comb_muts['Sample ID'].isin(comb_samples)]
    comb_cn = comb_cn[comb_cn['Sample ID'].isin(comb_samples)]
    
    
    pts = comb_cn['Patient ID'].unique().tolist()
    
    # =============================================================================
    # 
    # =============================================================================
    
    muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
    cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
    tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
    
    tc.columns = tc.iloc[0]
    tc = tc.drop(tc.index[0])
    tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
    
    tc = tc[~tc['Sample ID'].str.contains('cfDNA|MLN|MB|MT')]
    tc = tc[tc['Patient ID'].isin(pts)]
    
    tc = tc[tc['Final tNGS_TC'] > 0]
    tc_high = tc[tc['Final tNGS_TC'] > 0.2].copy()
    
    muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
    cn = cn[cn['Sample ID'].isin(tc_high['Sample ID'].tolist())]
    
    # =============================================================================
    # 
    # =============================================================================
    
    comb_muts = keepCodingMutations(comb_muts)
    muts = keepCodingMutations(muts)
    
    muts = muts[muts['GENE'].isin(additional_genes + cn['GENE'].unique().tolist())]
    comb_muts = comb_muts[comb_muts['GENE'].isin(additional_genes + cn['GENE'].unique().tolist())]
    
    genes = additional_genes + [g for g in muts['GENE'].unique().tolist() if g not in additional_genes]
    genes = genes + [g for g in comb_muts['GENE'].unique().tolist() if g not in genes]
    
    if comb_type == 'rp_comb':
        genes = additional_genes
    
    y_dict = dict(zip(genes, np.arange(len(genes))))
    
    # =============================================================================
    # format and plot combined sampled
    # =============================================================================
    
    comb_cn = comb_cn.merge(comb_tc[['Sample ID','Final tNGS_TC']], on = 'Sample ID', how = 'left')
    comb_cn = comb_cn.sort_values('Final tNGS_TC', ascending = False)
    
    x_dict = dict(zip(comb_cn['Patient ID'].unique().tolist(),\
                      np.arange(0, len(comb_cn['Patient ID'].unique())*3, 3)))
    
    comb_cn['x'] = comb_cn['Patient ID'].map(x_dict)
    comb_muts['x'] = comb_muts['Patient ID'].map(x_dict)
    
    comb_cn = comb_cn[comb_cn['GENE'].isin(genes)]
    comb_muts = comb_muts[comb_muts['GENE'].isin(genes)]
    
    comb_cn['y'] = comb_cn['GENE'].map(y_dict)
    comb_muts['y'] = comb_muts['GENE'].apply(lambda x: y_dict.get(x)+0.4)
    
    # =============================================================================
    # Assign colors to mutations and copy number
    # =============================================================================
    
    comb_cn['color'] = comb_cn['Copy_num'].map(cn_color)
    comb_muts['color'] = comb_muts['EFFECT'].str.split(' ').str[0].map(muts_color)
    
    # =============================================================================
    # Set up figure
    # =============================================================================
    
    size = 18 if comb_type == 'pb_comb' else 6
    fs = (7,6) if comb_type == 'pb_comb' else (7,2)
    
    fig,ax = plt.subplots(figsize = fs)
    
    ax.bar(comb_cn['x'], 0.8, bottom = comb_cn['y'], color = comb_cn['color'], zorder = 10)
    ax.scatter(comb_muts['x'], comb_muts['y'], color = comb_muts['color'], lw = 0, marker = 's', s = size, zorder = 1000)
    
    # =============================================================================
    # Plot aggregate of samples
    # =============================================================================
    
    # =============================================================================
    # Keep mutations and CNA in genes
    # =============================================================================
    
    muts = keepCodingMutations(muts)
    
    muts = muts[muts['GENE'].isin(genes)]
    cn = cn[cn['GENE'].isin(genes)]
    
    # =============================================================================
    # Get sample counts per patient for mutations
    # =============================================================================
    
    s_counts = tc.groupby('Patient ID').count()[['Sample ID']]
    s_counts.columns = ['pt_sample_count']
    
    muts = muts.merge(s_counts, left_on = 'Patient ID', right_index = True, how = 'left')
    
    # =============================================================================
    # Get sample counts per patient for CNA
    # =============================================================================
    
    high_tc_s_counts = tc_high.groupby('Patient ID').count()[['Sample ID']]
    high_tc_s_counts.columns = ['pt_sample_count']
    
    cn = cn.merge(high_tc_s_counts, left_on = 'Patient ID', right_index = True, how = 'left')
    
    # =============================================================================
    # Add mutation counts. Keep mutations where count > 3 or = num_samples
    # =============================================================================
    
    m_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
    m_counts.columns = ['Mutation count']
    m_counts = m_counts.reset_index()
    
    muts = muts.merge(m_counts, on = ['Patient ID','GENE','POSITION','EFFECT'], how = 'left')
    
    muts = muts[(muts['Mutation count'] >= 3)|((muts['pt_sample_count'] < 3)&(muts['pt_sample_count'] == muts['Mutation count']))]
    
    # =============================================================================
    # Add CNA counts. Keep CNA where count > 3. If not enough samples are present,
    # remove it. Including copy neutral
    # =============================================================================
    
    cn_counts = cn.groupby(['Patient ID','GENE','Copy_num']).count()[['Log_ratio']]
    cn_counts.columns = ['cn_count']
    cn_counts = cn_counts.reset_index()
    
    cn = cn.merge(cn_counts, on = ['Patient ID','GENE','Copy_num'])
    
    ## If using 3 as limit for AMP/DeepDel, then use this code so no shallow dels are lost
    '''
    for index, row in cn.iterrows():
        pt = row['Patient ID']
        g = row['GENE']
        if row['Copy_num'] == -2 and row['cn_count'] < 3:
            cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == -1), 'cn_count'] = cn['cn_count']+row['cn_count']
        elif row['Copy_num'] == 2 and row['cn_count'] < 3:
            cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == 1), 'cn_count'] = cn['cn_count'] + row['cn_count']
    cn = cn[cn['cn_count'] >= 3]
    '''
    # Keep amps and deep dels and gains/shallow losses if count >= 3
    cn = cn[(cn['Copy_num'].isin([-2,2]))|(cn['cn_count'] >= 3)]
    
    
    # =============================================================================
    # Assign x and y coordinates to mutations and CN events 
    # =============================================================================
    
    cn['y'] = cn['GENE'].map(y_dict)
    cn['x'] = cn['Patient ID'].apply(lambda x: x_dict.get(x)+1)
    
    muts['y'] = muts['GENE'].map(y_dict)
    muts['x'] = muts['Patient ID'].apply(lambda x: x_dict.get(x)+1)
    
    # =============================================================================
    # Drop dupbliate mutations
    # =============================================================================
    
    muts = muts.drop_duplicates(['Patient ID','GENE','EFFECT'])
    
    # =============================================================================
    # Find multiple mutations
    # =============================================================================
    
    m_count_pt = muts.groupby(['Patient ID','GENE']).count()[['Sample ID']]
    m_count_pt.columns = ['gene_mut_count']
    m_count_pt = m_count_pt.reset_index()
    
    muts = muts.merge(m_count_pt, on = ['Patient ID','GENE'], how = 'left')
    
    offset = 0.15
    
    for index, row in m_count_pt.iterrows():
        if row['gene_mut_count'] > 1:
            gene = row['GENE']
            pt = row['Patient ID']
            effects = muts[(muts['GENE'] == gene)&(muts['Patient ID'] == pt)]['EFFECT'].drop_duplicates().to_list()
            for i, effect in enumerate(effects):
                if i == 0:
                    # muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']+offset
                    muts.loc[(muts['GENE'] == gene) & (muts['EFFECT'] == effect) & (muts['Patient ID'] == pt), 'y'] = muts['y']+offset
                else:
                    # muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']-offset
                    muts.loc[(muts['GENE'] == gene) & (muts['EFFECT'] == effect) & (muts['Patient ID'] == pt), 'y'] = muts['y']-offset
    
    # =============================================================================
    # Add color
    # =============================================================================
    
    cn['color'] = cn['Copy_num'].map(cn_color)
    muts['color'] = muts['EFFECT'].str.split(' ').str[0].map(muts_color)
    
    # =============================================================================
    # CN 
    # =============================================================================
    
    cn_0 = cn[cn['Copy_num'] == 0].copy()
    ax.bar(cn_0['x'], 0.8, bottom = cn_0['y'], color = cn_0['color'], zorder = 10)
    
    cn_1 = cn[cn['Copy_num'].isin([-1,1])].copy()
    
    ax.bar(cn_1['x'], 0.8, bottom = cn_1['y'], color = cn_1['color'], zorder = 50)
    
    cn_2 = cn[cn['Copy_num'].isin([-2,2])]
    
    ax.bar(cn_2['x'], 0.8, bottom = cn_2['y'], color = cn_2['color'], zorder = 100)
    
    
    # =============================================================================
    # Muts
    # =============================================================================
    
    ax.scatter(muts['x'], muts['y']+0.4, c = muts['color'], marker = 's', lw = 0, zorder = 1000, s = size)
    
    # =============================================================================
    # Aethetics
    # =============================================================================
    
    ax.set_yticks(np.array(list(y_dict.values()))+0.4)
    ax.set_yticklabels(list(y_dict.keys()), fontsize = 6)
    
    ax.set_ylim(-0.1, len(genes))
    
    ax.set_xticks(np.array(list(x_dict.values()))+0.5)
    ax.set_xticklabels(list(x_dict.keys()), fontsize = 6, rotation = 90)
    ax.set_xlim(-0.65, len(pts)*3)
    
    ax.set_title('%s' % comb_type)
    
    plt.tight_layout()
    
    fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/%s_aggregate_oncoprint.pdf' % comb_type)
    fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/%s_aggregate_oncoprint.png' % comb_type)