# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:42:21 2021

@author: amurtha
"""

import pandas as pd
import string

cohort = 'M1RP'
snp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/ghent_M1RP_hetz_snps.vcf', sep = '\t')

snp = pd.melt(snp, id_vars = ['CHROM','POSITION','REF','ALT','NOTES'])

snp = snp.rename(columns = {'variable':'Sample ID'})
snp['Read depth'] = snp['value'].str.split(':').str[1].astype(float)
snp['Mutant reads'] = snp['value'].str.split(':').str[0].astype(float)
snp['VAF'] = snp['Mutant reads'] / snp['Read depth']
snp['Patient ID'] = snp['Sample ID'].str.split('_').str[1]
snp['Sample type'] = snp['Sample ID'].str.split('_').str[2].str.strip(string.digits)
snp.loc[snp['Sample ID'].str.contains('WBC'), 'Sample type'] = 'gDNA'

gl = snp.copy()
gl = gl[gl['Sample type'].isin(['NL','gDNA','WBC'])]
gl = gl[gl['value'].str.contains('\*')]
gl['Called'] = True

snp = snp[~snp['Sample type'].isin(['gDNA','WBC'])]

gl = gl[['CHROM','POSITION','REF','ALT','Patient ID','Called']]
snp = snp.merge(gl, on = ['CHROM','POSITION','REF','ALT','Patient ID'])

def label_genes(snp, gene):
    snp['GENE'] = None
    for index, row in gene.iterrows():
        start = row['START']-1000
        end = row['END']+1000
        chrom = row['CHROMOSOME']
        snp.loc[(snp['CHROM']==chrom)&(snp['POSITION']<=end)&(snp['POSITION']>=start), 'GENE'] = row['GENE']
    return snp;

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_FFPE_cna.tsv', sep = '\t')

cn = cn[['GENE','CHROMOSOME','START','END']].drop_duplicates()

snp = label_genes(snp,cn)

snp.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/M1RP_snps_melted.tsv', sep = '\t')