# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:26:49 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Constants
# =============================================================================

genes = ['TP53','PTEN','RB1','BRCA2','AR','NCOA2','MYC','NKX3-1','CLU','CHD1','CDKN2A']
alpha = 0.001

colors = {
    4:'#EE2D24',
    3:'#F59496',
    2:'#E6E7E8',
    1:'#9CC5E9',
    0:'#3F60AC',
    -1:'#ffffff',
    }

# =============================================================================
# Helpers
# =============================================================================

# Return the copy number for the gene sample pair. 2 if TC > min_tc_1loss and 
# p_val > significance. 0 if p_val is significant and TC < min_tc_1loss and 
# lr < lr_mean. 1 if if p_val is significant and TC > min_tc_1loss and 
# lr < lr_mean. 3 if p_val is significant and TC > min_tc_1gain and 
# lr > lr_mean. 4 if p_val is significant and TC < min_tc_1gain and 
# lr > lr_mean. -1 if TC < min_tc_1loss and p is not significant
def get_cn(cn, gene, sample, alpha):
    print(gene, sample)
    row = cn.loc[(cn['GENE'] == gene)&(cn['Sample ID'] == sample)].reset_index().iloc[0]
    pval = row['p_val']
    if pval < alpha:
        print(True)
        if row['log_ratio'] > row['lr_mean']:
            if row['mut_TC'] > row['min_tc_1gain']:
                return 3;
            else:
                return 4;
        else:
            if row['mut_TC'] > row['min_tc_1loss']:
                return 1;
            else:
                return 0
    else:
        if row['mut_TC'] < row['min_tc_1loss']:
            return -1
        else:
            return 2;
    

def create_matrix(cn, genes, samples, alpha):
    matrix = pd.DataFrame(index = genes, columns = samples)
    for gene in genes:
        for sample in samples:
            matrix.at[gene, sample] = get_cn(cn, gene, sample, alpha)
    return matrix;

def visualize_matrix(matrix, patient):
    fig,ax = plt.subplots()
    for x, sample in enumerate(matrix.columns):
        for y, gene in enumerate(matrix.index):
            ax.bar(x, 0.8, bottom = y, color = colors.get(matrix.at[gene, sample]))
            
    ax.set_title(matrix.columns.tolist()[0].split('_')[1])
    ax.set_xticks(np.arange(0,len(matrix.columns),1))
    ax.set_xticklabels(matrix.columns.str.split('_').str[2], rotation = 90, fontsize = 8)
    ax.set_yticks(np.arange(0.4,len(matrix),1))
    ax.set_yticklabels(matrix.index, fontsize = 8)
    ax.tick_params(bottom = False, left = False)
    fig.tight_layout()
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/CN grouped by patient (alpha=0.001)/%s.pdf' % patient)

# =============================================================================
# Main
# =============================================================================

'''
Create an oncoprint for the specific genes listed above. Copy number changes
are shown based on the alpha threshold above, tumor content, and direction of 
change
'''

# =============================================================================
# Import data
# =============================================================================

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Keep ctDNA and FFPE samples
# =============================================================================

samples = samples[samples['Sample type'].isin(['FFPE tissue'])]
samples = samples[samples['Cohort'] == 'M1RP']
samples = samples[samples['Targeted sequencing'] == 'Completed']

# =============================================================================
# Keep genes designated above in the samples in samples
# =============================================================================

cn = cn[cn['GENE'].isin(genes)]
cn = cn[cn['Sample ID'].isin(samples['Sample ID'].unique().tolist())]

# =============================================================================
# Loop over patients and create a CN matix and oncoprint for each patient
# =============================================================================

for patient in cn['Patient ID'].unique().tolist():
    tmp_cn = cn[cn['Patient ID'] == patient]
    pt_samples = samples[samples['Patient ID'] == patient]['Sample ID']
    matrix = create_matrix(tmp_cn, genes, pt_samples, alpha)
    visualize_matrix(matrix, patient)
