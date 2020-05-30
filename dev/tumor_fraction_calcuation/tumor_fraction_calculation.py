# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:29:23 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

missing_gDNA = {'M1B':['ID22','ID27','ID1','ID31', 'ID19','ID21'],
                'M1RP':[]}
cohort = 'M1RP'
min_reads = 75
lr_max = 0.3
cna_exclusion = ['FOXA1']

# =============================================================================
# Import data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/%s_mutations.tsv' % cohort, sep = '\t')
tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/%s_cna_melted.tsv' % cohort, sep = '\t')

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction['mut_TC'].str.split('%').str[0].astype(np.float64) / 100
muts['Cohort'] = cohort

tumor_fraction = tumor_fraction[tumor_fraction['Cohort'] == cohort]

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
muts.loc[(muts['Log_ratio']<=lr_max)|(muts['GENE'].isin(cna_exclusion)),'Lr_flag']=True

# =============================================================================
# Flag mutation on allosomes. True if on autosome. 
# =============================================================================

muts['Allosome_flag'] = True
muts.loc[muts['CHROM'].isin(['chrX','chrY']), 'Allosome_flag'] = False

# =============================================================================
# Flag patients with missing gDNA. False if no gDNA
# =============================================================================

muts['gDNA_flag'] = True
muts.loc[muts['Patient ID'].isin(missing_gDNA.get(cohort)), 'gDNA_flag'] = False

# =============================================================================
# Calculate mut_TC for each mutation
# =============================================================================

muts['mut_TC'] = 2/(1+1/muts['Allele_frequency'])
muts.loc[muts['Allosome_flag'] == False, 'mut_TC'] = muts['Allele_frequency']

# =============================================================================
# Sort mutations based on sample id, flags, tc. Merge in samples with no mutations
# =============================================================================

muts = muts.merge(tumor_fraction[['Cohort','Patient ID','Sample ID']], on = ['Cohort','Patient ID','Sample ID'], how = 'right')

muts['mut_TC'] = muts['mut_TC'].fillna(0)

muts = muts.sort_values(['Sample ID','Depth_flag','Lr_flag','Allosome_flag','gDNA_flag', 'mut_TC'], ascending = [True,False,False,False,False,False])

# =============================================================================
# Reorder columns and save as full mutation file
# =============================================================================

muts = muts[['Cohort', 'Patient ID', 'Sample ID','mut_TC', 'CHROM', 'POSITION', 'REF', 'ALT','GENE', 'EFFECT', 'Allele_frequency', 'Read_depth', 'NOTES', 'Log_ratio', 'Depth_flag', 'Lr_flag', 'Allosome_flag', 'gDNA_flag']]

muts.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_allMuts.tsv' % cohort, sep = '\t', index = None)

# =============================================================================
# Keep only top mutation per patient. All flags must be true
# =============================================================================

muts = muts[(muts['Depth_flag'] == True)&(muts['Lr_flag'] == True)&(muts['Allosome_flag'] == True)&(muts['gDNA_flag'] == True)]
muts = muts.drop_duplicates('Sample ID')

muts = muts.merge(tumor_fraction[['Cohort','Patient ID','Sample ID']], on = ['Cohort','Patient ID','Sample ID'], how = 'right')

muts = muts[['Cohort', 'Patient ID', 'Sample ID', 'mut_TC', 'CHROM', 'POSITION', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','Log_ratio']]
muts.columns = ['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position','Gene',	'Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio']

muts.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction.tsv' % cohort, sep = '\t', index = None)