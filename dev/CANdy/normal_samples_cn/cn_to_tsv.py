# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:07:11 2020

@author: amurtha
"""
import pandas as pd

def meltCN(filepath):
    df = pd.read_csv(filepath, delimiter = '\t', index_col=None)
    df = pd.melt(df, id_vars=['GENE', 'CHROMOSOME', 'START', 'END'])
    df.rename(columns={'value': 'Copy_num'}, inplace=True)
    df.rename(columns={'variable': 'Sample ID'}, inplace=True)
    df['Cohort'] = df['Sample ID'].str.split('_').str[0]
    df['Patient ID'] = df['Sample ID'].str.split('_').str[1]
    df['Log_ratio'] = df['Copy_num'].str.split(':').str[1]
    df['Copy_num'] = df['Copy_num'].str.split(':').str[0]
    df[['Copy_num','Log_ratio']] = df[['Copy_num','Log_ratio']].apply(pd.to_numeric)
    df = df[['Cohort','Patient ID','Sample ID', 'GENE', 'Copy_num', 'Log_ratio', 'CHROMOSOME', 'START', 'END']]
    return df;

m1b = meltCN('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1B Copy Number Analysis/Targeted (Ghent normals)/ghentm1b_gene_cna.tsv')
m1rp = meltCN('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Copy number analysis (targeted, using Ghent normals)/gene_cna.ffpe.tsv')

m1rp = m1rp[~m1rp['Sample ID'].str.contains('cfDNA')]
m1b = m1b[~m1b['Sample ID'].str.contains('cfDNA')]

m1rp.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/M1RP_FFPE_cna.tsv', sep = '\t', index = None)
m1b.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/M1B_FFPE_cna.tsv', sep = '\t', index = None)