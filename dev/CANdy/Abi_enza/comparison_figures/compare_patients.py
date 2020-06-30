# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:17:59 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# =============================================================================
# Constants
# =============================================================================

genes = ['TP53','PTEN','RB1','BRCA2','AR','NCOA2','MYC','NKX3-1','CLU','CHD1','CDKN2A']

cohort = 'M1RP'

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
def get_cn(cn, gene, sample):
    row = cn.loc[(cn['GENE'] == gene)&(cn['Sample ID'] == sample)].reset_index().iloc[0]
    return row['Adjusted_copy_num']
    

def create_matrix(cn, genes, samples):
    matrix = pd.DataFrame(index = genes, columns = samples)
    for gene in genes:
        for sample in samples:
            matrix.at[gene, sample] = get_cn(cn, gene, sample)
    return matrix;

def visualize_matrix(matrix_nas, matrix_mat, patient):
    fig,[ax1,ax2] = plt.subplots(ncols = 2, sharey = True)
    for x, sample in enumerate(matrix_nas.columns):
        for y, gene in enumerate(matrix_nas.index):
            ax1.bar(x, 0.8, bottom = y, color = colors.get(matrix_nas.at[gene, sample]))
            ax2.bar(x, 0.8, bottom = y, color = colors.get(matrix_mat.at[gene, sample]))
            
    ax1.set_title(matrix_nas.columns.tolist()[0].split('_')[1] + ' NASCNT')
    ax1.set_xticks(np.arange(0,len(matrix_nas.columns),1))
    ax1.set_xticklabels(matrix_nas.columns.str.split('_').str[2], rotation = 90, fontsize = 8)
    ax2.set_title(matrix_mat.columns.tolist()[0].split('_')[1] + ' MATTI')
    ax2.set_xticks(np.arange(0,len(matrix_mat.columns),1))
    ax2.set_xticklabels(matrix_mat.columns.str.split('_').str[2], rotation = 90, fontsize = 8)
    ax1.set_yticks(np.arange(0.4,len(matrix_nas),1))
    ax1.set_yticklabels(matrix_nas.index, fontsize = 8)
    ax1.tick_params(bottom = False, left = False)
    fig.tight_layout()
    plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/comparison_figures/patient_comparisons/%s.pdf' % patient)

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

nas = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_%s_FFPE_cna.tsv' % cohort, sep = '\t')
mat = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_FFPE_cna.tsv' % cohort, sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

mat['Adjusted_copy_num'] = mat['Copy_num'].replace({-2:0,-1:1,0:2,1:3,2:4})

# =============================================================================
# Keep ctDNA and FFPE samples
# =============================================================================

samples = samples[samples['Sample type'].isin(['FFPE tissue'])]
samples = samples[samples['Cohort'] == cohort]
samples = samples[samples['Targeted sequencing'] == 'Completed']

# =============================================================================
# Keep genes designated above in the samples in samples
# =============================================================================

nas = nas[nas['GENE'].isin(genes)]
nas = nas[nas['Sample ID'].isin(samples['Sample ID'].unique().tolist())]

mat = mat[mat['GENE'].isin(genes)]
mat = mat[mat['Sample ID'].isin(samples['Sample ID'].unique().tolist())]

# =============================================================================
# Loop over patients and create a CN matix and oncoprint for each patient
# =============================================================================

for patient in nas['Patient ID'].unique().tolist():
    tmp_cn_nas = nas[nas['Patient ID'] == patient]
    tmp_cn_mat = mat[mat['Patient ID'] == patient]
    pt_samples = samples[samples['Patient ID'] == patient]['Sample ID']
    matrix_nas = create_matrix(tmp_cn_nas, genes, pt_samples)
    matrix_mat = create_matrix(tmp_cn_mat, genes, pt_samples)
    visualize_matrix(matrix_nas, matrix_mat, patient)
