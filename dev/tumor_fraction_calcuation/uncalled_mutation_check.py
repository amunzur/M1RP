# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:07:10 2020

@author: amurtha

"""

import pandas as pd
import sys
import numpy as np

# =============================================================================
# Constants and parameters
# =============================================================================

missing_gDNA = {'M1B':['ID22','ID27','ID1','ID31', 'ID19','ID21'],
                'M1RP':[]}
cohort = sys.argv[1]
min_reads = int(sys.argv[2])
lr_max = float(sys.argv[3])
cna_exclusion = ['FOXA1']
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']

# =============================================================================
# Helpers
# =============================================================================

def meltBet(path):
    df = pd.read_excel(path, index_col=None)
    df = pd.melt(df, id_vars=['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES'])
    df.rename(columns={'value': 'Allele_frequency'}, inplace=True)
    df.rename(columns={'variable': 'Sample ID'}, inplace=True)
    df['Read_depth'] = df['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[1]
    df = df[~df['Read_depth'].str.contains("*", regex=False)]
    df['Read_depth'] = df['Read_depth'].replace('\(','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\)','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\*','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\*','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace("\[[^]]*\]",'', regex=True)
    df['Allele_frequency'] = df['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[0]
    df['Cohort'] = cohort
    df['Patient ID'] = df['Sample ID'].str.split('_').str[1]
    df = df[['Cohort','Patient ID','Sample ID', 'CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','NOTES']]
    df[['Read_depth','Allele_frequency']] = df[['Read_depth', 'Allele_frequency']].apply(pd.to_numeric)
    df['Allele_frequency'] = df['Allele_frequency'] / 100.
    return df;

# =============================================================================
# Import called mutations, tumor fraction estimations, and full mutation table
# =============================================================================
    
muts_uncalled = meltBet('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/%s_betastasis_all.xlsx' % cohort)
called = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/%s_mutations.tsv' % cohort, sep = '\t')
mut_tc = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_allMuts.xlsx' % cohort)
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/%s_cna_melted.tsv' % cohort, sep = '\t')
exclude = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_exclusion.tsv' % cohort, sep = '\t')


# =============================================================================
# Get unique mutations
# =============================================================================

called = called.drop_duplicates(unique_mutation_columns)

# =============================================================================
# Limit to > 1% AF and 3 mutant reads
# =============================================================================

muts_uncalled = muts_uncalled[muts_uncalled['Allele_frequency'] >= 0.01]
muts_uncalled['Mutant reads'] = muts_uncalled['Allele_frequency'] * muts_uncalled['Read_depth']
muts_uncalled = muts_uncalled[muts_uncalled['Mutant reads'] >= 3]

# =============================================================================
# Keep mutation in called. Must be same patient
# =============================================================================

called = called.set_index(unique_mutation_columns)
muts_uncalled = muts_uncalled.set_index(unique_mutation_columns)

muts_uncalled = muts_uncalled[muts_uncalled.index.isin(called.index)]

# =============================================================================
# Keep only uncalled mutations in patients with no tumor content
# =============================================================================

mut_tc = mut_tc[mut_tc['Cohort'] == cohort]
tc_neg = mut_tc[mut_tc['mut_TC'] == 0]
muts_uncalled = muts_uncalled[muts_uncalled['Sample ID'].isin(tc_neg['Sample ID'])]

# =============================================================================
# Add flags
# =============================================================================

muts_uncalled = muts_uncalled.reset_index()

# =============================================================================
# Merge log_ratio onto mutations
# =============================================================================

muts_uncalled = muts_uncalled.merge(cn[['Sample ID','GENE','Log_ratio']], on = ['Sample ID','GENE'], how = 'left')

# =============================================================================
# Flag mutation with < min reads. Depth flag = False if < min_reads
# =============================================================================

muts_uncalled['Depth_flag'] = True
muts_uncalled.loc[muts_uncalled['Read_depth'] < min_reads, 'Depth_flag'] = False

# =============================================================================
# Flag mutations on amplified genes (log_ratio > lr_max). Exclude mutations on 
# genes in cna_exclusion
# Lr_flag is false if mutation is on an included amplified gene. Lr_flag also 
# gets false if mutation is intergenic or off target
# =============================================================================

muts_uncalled['Lr_flag'] = False
muts_uncalled.loc[(muts_uncalled['Log_ratio']<=lr_max)|(muts_uncalled['GENE'].isin(cna_exclusion)),'Lr_flag']=True

# =============================================================================
# Flag mutation on allosomes. True if on autosome. 
# =============================================================================

muts_uncalled['Allosome_flag'] = True
muts_uncalled.loc[muts_uncalled['CHROM'].isin(['chrX','chrY']), 'Allosome_flag'] = False

# =============================================================================
# Flag patients with missing gDNA. False if no gDNA
# =============================================================================

muts_uncalled['gDNA_flag'] = True
muts_uncalled.loc[muts_uncalled['Patient ID'].isin(missing_gDNA.get(cohort)), 'gDNA_flag'] = False

# =============================================================================
# Calculate tumor content
# =============================================================================

muts_uncalled['mut_TC'] = 2 / (1+1/muts_uncalled['Allele_frequency'])
muts_uncalled.loc[muts_uncalled['Allosome_flag'] == False, 'mut_TC'] = muts_uncalled['Allele_frequency']

# =============================================================================
# Merge exlcusion list onto uncalled
# =============================================================================

muts_uncalled = muts_uncalled.merge(exclude, on = ['Sample ID','GENE','POSITION','CHROM','REF','ALT','EFFECT'], how = 'left')

# =============================================================================
# Remerge and sort the mut_tc table
# =============================================================================

mut_tc = mut_tc[~mut_tc['Sample ID'].isin(muts_uncalled['Sample ID'])]
muts_uncalled = muts_uncalled[mut_tc.columns.tolist()]

mut_tc = mut_tc.append(muts_uncalled, ignore_index = True)

mut_tc['Manual_curation'] = mut_tc['Manual_curation'].fillna(1)

mut_tc = mut_tc.sort_values(['Sample ID','Manual_curation','mut_TC','Allosome_flag','Depth_flag','Lr_flag', 'gDNA_flag'], ascending = [True,False,False,False,False,False,False])

mut_tc = mut_tc[['Cohort', 'Patient ID', 'Sample ID','mut_TC', 'CHROM', 'POSITION', 'REF', 'ALT','GENE', 'EFFECT', 'Allele_frequency', 'Read_depth', 'NOTES', 'Log_ratio', 'Allosome_flag','Depth_flag','Lr_flag','gDNA_flag','Manual_curation']]

mut_tc['Manual_curation'] = mut_tc['Manual_curation'].replace(1, np.nan)

# =============================================================================
# Label first mutation
# =============================================================================

mut_tc['TC_call'] = False

unique_samples = mut_tc['Sample ID'].drop_duplicates().tolist()
sample_id_dict = dict(zip(unique_samples, [0]*len(unique_samples)))

for index, row in mut_tc.iterrows():
    if sample_id_dict.get(row['Sample ID']) == 0:
        mut_tc.at[index, 'TC_call'] = True
        sample_id_dict[row['Sample ID']] = 1

# =============================================================================
# Save to excel
# =============================================================================

mut_tc.to_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_allMuts.xlsx' % cohort, index = None)

# =============================================================================
# Keep only top mutation per patient. All flags must be true
# =============================================================================

mut_tc = mut_tc[mut_tc['gDNA_flag'] == True]
mut_tc = mut_tc[mut_tc['Manual_curationl'] != 0]
mut_tc = mut_tc.drop_duplicates('Sample ID')

mut_tc = mut_tc[['Cohort', 'Patient ID', 'Sample ID', 'mut_TC', 'CHROM', 'POSITION', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','Log_ratio']]
mut_tc.columns = ['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position','Gene','Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio']

mut_tc = mut_tc.sort_values('Sample ID')

mut_tc.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction.tsv' % cohort, sep = '\t', index = None)