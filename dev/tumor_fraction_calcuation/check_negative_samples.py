# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:35:56 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction['mut_TC'] = tumor_fraction['mut_TC'].str.split('%').str[0].astype(np.float64) / 100

tumor_fraction = tumor_fraction[tumor_fraction['mut_TC'] == 0]

muts = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/betastasis_M1B_noncoding.xlsx', index_col=None)
muts = pd.melt(muts, id_vars=['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES'])
muts.rename(columns={'value': 'Allele_frequency'}, inplace=True)
muts.rename(columns={'variable': 'Sample ID'}, inplace=True)
muts['Patient ID'] = muts['Sample ID'].str.split('_').str[1]
muts['Read_depth'] = muts['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[1]
muts = muts[muts['Read_depth'].str.contains("*", regex=False)]
muts['Read_depth'] = muts['Read_depth'].replace('\(','', regex=True)
muts['Read_depth'] = muts['Read_depth'].replace('\)','', regex=True)
muts['Read_depth'] = muts['Read_depth'].replace('\*','', regex=True)
muts['Read_depth'] = muts['Read_depth'].replace('\*','', regex=True)
muts['Read_depth'] = muts['Read_depth'].replace("\[[^]]*\]",'', regex=True)
muts['Allele_frequency'] = muts['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[0]
muts = muts[['Patient ID','Sample ID', 'CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','NOTES']]
muts[['Read_depth','Allele_frequency']] = muts[['Read_depth','Allele_frequency']].apply(pd.to_numeric)

muts = muts[muts['Sample ID'].isin(tumor_fraction['Sample ID'])]