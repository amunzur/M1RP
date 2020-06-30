# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:29:23 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import sys

missing_gDNA = {'M1B':['ID22','ID27','ID1','ID31', 'ID19','ID21'],
                'M1RP':[]}

if len(sys.argv) > 1:
    cohort = sys.argv[1]
    min_reads = int(sys.argv[2])
    lr_max = float(sys.argv[3])
else:
    cohort = 'M1RP'
    min_reads = 75
    lr_max = 0.3
    
cna_exclusion = []

# =============================================================================
# Import data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/%s_mutations.tsv' % cohort, sep = '\t')
tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_FFPE_cna.tsv' % cohort, sep = '\t')
cn_tmp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_cfDNA_cna.tsv' % cohort, sep = '\t')
exclude = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_allMuts.xlsx' % cohort)

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction['mut_TC'].str.split('%').str[0].astype(np.float64) / 100
muts['Cohort'] = cohort

tumor_fraction = tumor_fraction[tumor_fraction['Cohort'] == cohort]

cn = cn.append(cn_tmp, ignore_index = True)
del cn_tmp

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
# Merge old exclusion column onto dataframe
# =============================================================================

exclude = exclude[['Sample ID','GENE','POSITION','CHROM','REF','ALT','EFFECT', 'Manual_curation']]

exclude['Manual_curation'] = exclude['Manual_curation'].fillna(1)

exclude.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_exclusion.tsv' % cohort, sep = '\t', index = None)

muts = muts.merge(exclude, on = ['Sample ID','GENE','POSITION','CHROM','REF','ALT','EFFECT'], how = 'left')

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

muts = muts.sort_values(['Sample ID','Manual_curation','mut_TC','Allosome_flag', 'Depth_flag', 'Lr_flag', 'gDNA_flag'], ascending = [True,False, False,False,False,False, False])

muts['Manual_curation'] = muts['Manual_curation'].replace(1, np.nan)

# =============================================================================
# Reorder columns and save as full mutation file
# =============================================================================

muts = muts[['Cohort', 'Patient ID', 'Sample ID','mut_TC', 'CHROM', 'POSITION', 'REF', 'ALT','GENE', 'EFFECT', 'Allele_frequency', 'Read_depth', 'NOTES', 'Log_ratio', 'Allosome_flag','Depth_flag','Lr_flag','gDNA_flag','Manual_curation']]

muts.to_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_allMuts.xlsx' % cohort, index = None)

# =============================================================================
# Keep only top mutation per patient. All flags must be true
# =============================================================================

muts = muts[muts['gDNA_flag'] == True]
muts = muts.drop_duplicates('Sample ID')

muts = muts.merge(tumor_fraction[['Cohort','Patient ID','Sample ID']], on = ['Cohort','Patient ID','Sample ID'], how = 'right')

muts = muts[['Cohort', 'Patient ID', 'Sample ID', 'mut_TC', 'CHROM', 'POSITION', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','Log_ratio']]
muts.columns = ['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position','Gene',	'Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio']

muts = muts.sort_values('Sample ID')

muts.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction.tsv' % cohort, sep = '\t', index = None)