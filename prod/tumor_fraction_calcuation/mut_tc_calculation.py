# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:38:39 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import sys

if len(sys.argv) > 1:
    cohort = sys.argv[1]
    min_reads = int(sys.argv[2])
    lr_max = float(sys.argv[3])
    min_truncal = 0.75
else:
    cohort = 'M1RP'
    min_reads = 30
    lr_max = 0.3
    min_truncal = 0.75

# =============================================================================
# Import data
# =============================================================================

called = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations.tsv', sep = '\t')
muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/%s_mutations_inclDependent.tsv' % cohort, sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/%s_FFPE_cna.tsv' % cohort, sep = '\t')
cn_tmp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/%s_cfDNA_cna.tsv' % cohort, sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['mut_TC'].str.split('%').str[0].astype(np.float64) / 100
muts['Cohort'] = cohort

tc = tc[tc['Cohort'] == cohort]

cn = cn.append(cn_tmp, ignore_index = True)
del cn_tmp

muts.loc[muts['GENE'] == 'BIVM-ERCC5;ERCC5', 'GENE'] = 'ERCC5'

# =============================================================================
# Merge log_ratio onto mutations
# =============================================================================

muts = muts.merge(cn[['Sample ID','GENE','Log_ratio']], on = ['Sample ID','GENE'], how = 'left')

# =============================================================================
# Flag mutation with < min reads. Depth flag = False if < min_reads
# =============================================================================

muts['Depth_flag'] = True
muts.loc[muts['Read_depth'] < min_reads, 'Depth_flag'] = False

# =============================================================================
# Flag mutations on amplified genes (log_ratio > lr_max). Exclude mutations on 
# genes in cna_exclusion
# Lr_flag is false if mutation is on an included amplified gene. Lr_flag also 
# gets false if mutation is intergenic or off target
# =============================================================================

muts['Lr_flag'] = False
muts.loc[muts['Log_ratio']<=lr_max,'Lr_flag']=True

# =============================================================================
# Flag mutation on allosomes. True if on autosome. 
# =============================================================================

muts['Allosome_flag'] = True
muts.loc[muts['CHROM'].isin(['chrX','chrY']), 'Allosome_flag'] = False

# =============================================================================
# Create % column and flag mutations in less that 75% of patient matched samples
# =============================================================================

sample_counts = called[['Sample ID','Patient ID']].drop_duplicates().groupby('Patient ID').count().reset_index()
mutation_counts = called.groupby(['Patient ID','CHROM','POSITION','GENE','EFFECT']).count()[['Cohort']].reset_index()

mutation_counts = mutation_counts.merge(sample_counts, on = 'Patient ID', how = 'left')

mutation_counts = mutation_counts.rename(columns = {'Cohort':'Mutation count','Sample ID':'Sample count'})

mutation_counts['Percent'] = mutation_counts['Mutation count'] / mutation_counts['Sample count']

mutation_counts = mutation_counts.drop(['Mutation count','Sample count'], axis = 1)

muts = muts.merge(mutation_counts, on = ['Patient ID','CHROM','POSITION','GENE','EFFECT'], how = 'left' )

muts['Non-truncal flag'] = True
muts.loc[muts['Percent'] < min_truncal, 'Non-truncal flag'] = False

# =============================================================================
# Calculate muts for each mutation
# =============================================================================

muts['mut_TC'] = 2/(1+1/muts['Allele_frequency'])
muts.loc[muts['Allosome_flag'] == False, 'mut_TC'] = muts['Allele_frequency']

# =============================================================================
# Sort mutations based on sample id, flags, tc. Merge in samples with no mutations
# =============================================================================

muts = muts.sort_values(['Sample ID','Lr_flag','Allosome_flag','Depth_flag','Non-truncal flag','mut_TC'], ascending = [True,False,False,False,False,False])

muts = muts[['Cohort', 'Patient ID', 'Sample ID','mut_TC', 'CHROM', 'POSITION', 'REF', 'ALT','GENE', 'EFFECT', 'Allele_frequency', 'Read_depth', 'NOTES', 'Log_ratio', 'Allosome_flag','Depth_flag','Lr_flag','Non-truncal flag']]

# =============================================================================
# Label first mutation
# =============================================================================

muts['TC_call'] = False

unique_samples = muts['Sample ID'].drop_duplicates().tolist()
sample_id_dict = dict(zip(unique_samples, [0]*len(unique_samples)))

for index, row in muts.iterrows():
    if sample_id_dict.get(row['Sample ID']) == 0:
        muts.at[index, 'TC_call'] = True
        sample_id_dict[row['Sample ID']] = 1

# =============================================================================
# Save to excel
# =============================================================================

muts.to_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_mut_tc_allMuts.xlsx' % cohort, index = None)

# =============================================================================
# Keep only top mutation per patient. All flags must be true
# =============================================================================

muts = muts[muts['Lr_flag'] != False]
muts = muts[muts['Depth_flag'] != False]
muts = muts.drop_duplicates('Sample ID')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])
samples = samples[samples['Cohort'] == cohort]
muts = muts.merge(samples[['Cohort','Patient ID','Sample ID']], on = ['Cohort','Patient ID','Sample ID'], how = 'right')

muts['mut_TC'] = muts['mut_TC'].fillna(0)

muts = muts[['Cohort', 'Patient ID', 'Sample ID', 'mut_TC', 'CHROM', 'POSITION', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','Log_ratio','Non-truncal flag']]
muts.columns = ['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position','Gene','Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio','Non-truncal flag']

muts = muts.sort_values('Sample ID')

muts.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_mut_tc.tsv' % cohort, sep = '\t', index = None)