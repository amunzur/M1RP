# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:31:49 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

cohort = 'M1B'

# =============================================================================
# Import data
# =============================================================================

bet = None
if cohort == 'M1RP':
    bet = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/betastasis_noncoding.xlsx')
elif cohort == 'M1B':
    bet = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/betastasis_M1B_noncoding.xlsx')

# =============================================================================
# Melt betastasis
# =============================================================================

bet = pd.melt(bet, id_vars=['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES'])
bet.rename(columns={'value': 'Allele_frequency'}, inplace=True)
bet.rename(columns={'variable': 'Sample ID'}, inplace=True)
bet['Patient ID'] = bet['Sample ID'].str.split('_').str[1]
bet['Read_depth'] = bet['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[1]
bet['Called'] = False
bet.loc[bet['Read_depth'].str.contains("*", regex=False), 'Called'] = True
bet['Read_depth'] = bet['Read_depth'].replace('\(','', regex=True)
bet['Read_depth'] = bet['Read_depth'].replace('\)','', regex=True)
bet['Read_depth'] = bet['Read_depth'].replace('\*','', regex=True)
bet['Read_depth'] = bet['Read_depth'].replace('\*','', regex=True)
bet['Read_depth'] = bet['Read_depth'].replace("\[[^]]*\]",'', regex=True)
bet['Allele_frequency'] = bet['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[0].astype(float) / 100
bet = bet[['Patient ID','Sample ID', 'CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','NOTES', 'Called']]
bet[['Read_depth','Allele_frequency']] = bet[['Read_depth','Allele_frequency']].apply(pd.to_numeric)


# =============================================================================
# Label oxidatative
# =============================================================================

bet['Damage'] = False
bet.loc[(bet['REF'] == 'C')&(bet['ALT'] == 'A'), 'Damage'] = True
bet.loc[(bet['REF'] == 'G')&(bet['ALT'] == 'T'), 'Damage'] = True

# =============================================================================
# Keep mutations above 5%
# =============================================================================

bet = bet[bet['Allele_frequency'] >= 0.05]
called  = bet[bet['Called'] == True].copy()

# =============================================================================
# Create blacklist. Keep damaged mutation at AF <= 0.15
# =============================================================================

blacklist = called.copy()
blacklist = blacklist[blacklist['Damage'] == True]
blacklist = blacklist[blacklist['Allele_frequency'] <= 0.15]
blacklist = blacklist.set_index(['Patient ID','CHROM','POSITION','EFFECT','GENE'])

# =============================================================================
# Get a count of each unique mutation per patient
# =============================================================================

mut_counts = bet.copy()
mut_counts = mut_counts[mut_counts['Damage'] == True]
mut_counts = mut_counts.groupby(['Patient ID','CHROM','POSITION','EFFECT','GENE']).count()
mut_counts = mut_counts[['Sample ID']].rename(columns = {'Sample ID':'Count'})

# =============================================================================
# Keep mutations in blacklist if they only appear once in mut_counts
# =============================================================================

mut_counts = mut_counts[mut_counts['Count'] == 1]
blacklist = blacklist[blacklist.index.isin(mut_counts.index)]

# =============================================================================
# Check the cosmic score
# =============================================================================

blacklist['NOTES'] = blacklist['NOTES'].str.split('\. ').str[1]
for index, row in blacklist.iterrows():
    if 'COSMIC' in row['NOTES']:
        cosmic = int(row['NOTES'].split('COSMIC:')[1].split('.')[0])
        if cosmic > 10:
            blacklist.drop(index)
            
blacklist = blacklist[~blacklist['Sample ID'].str.contains('cfDNA')]

blacklist.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/blacklists/%s_blacklist_full.tsv' % cohort, sep = '\t', header = None, index = None)

blacklist = blacklist.reset_index()[['CHROM','POSITION','REF','ALT']]

blacklist.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/blacklists/%s_blacklist_upload.tsv' % cohort, sep = '\t', header = None, index = None)