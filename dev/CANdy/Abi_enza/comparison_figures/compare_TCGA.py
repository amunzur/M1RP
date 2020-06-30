# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:07:32 2020

@author: amurtha
"""

import pandas as pd
import random

# =============================================================================
# Helpers
# =============================================================================



def meltCN(df):
    df = pd.melt(df, id_vars=['track_name', 'track_type'])
    df.rename(columns={'value': 'Copy_num'}, inplace=True)
    df.rename(columns={'variable': 'Sample_ID'}, inplace=True)
    df = df.replace({'Amplification':4,'Deep Deletion':0})
    df['Copy_num'] = df['Copy_num'].fillna(0)
    df = df.drop('track_type', axis = 1)
    return df;

# =============================================================================
# Import data
# =============================================================================

tcga = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/comparison_figures/TCGA_cn_cell2015.tsv', sep = '\t')

nas = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_M1RP_FFPE_cna.tsv', sep = '\t')

# =============================================================================
# Melt tcga
# =============================================================================

tcga = meltCN(tcga)

# =============================================================================
# Create patient sample dict
# =============================================================================

pt_sample_dict = {}
for pt in nas['Patient ID'].unique().tolist():
    samples = nas[nas['Patient ID'] == pt]['Sample ID'].unique().tolist()
    pt_sample_dict[pt] = samples
del pt, samples

# =============================================================================
# Get random m1rp samples
# =============================================================================

samples = [random.choice(pt_sample_dict.get(pt)) for pt in pt_sample_dict.keys()]

nas = nas[nas['Sample ID'].isin(samples)] 

nas = nas[nas['Adjusted_copy_num'] != -1]

gene_counts = nas.groupby('GENE').count()
cna_counts = nas.groupby(['GENE','Adjusted_copy_num']).count()