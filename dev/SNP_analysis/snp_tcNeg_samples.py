# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:03:16 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np


# =============================================================================
# Constants
# =============================================================================

min_reads = 100
alpha = 0.05
min_snps = 5
window_size = 25 # In each direction
cohort = 'm1rp'

# =============================================================================
# Helpers
# =============================================================================

def label_genes(snp, gene):
    snp['GENE'] = None
    for index, row in gene.iterrows():
        start = row['START']-1000
        end = row['END']+1000
        chrom = row['CHROMOSOME']
        snp.loc[(snp['CHROM']==chrom)&(snp['POSITION']<=end)&(snp['POSITION']>=start), 'GENE'] = row['GENE']
    return snp;

def create_depth_std(gl, window_size):
    stds = []
    means= []
    for i in np.arange(min_reads, 3000, 1):
        tmp = gl[(gl['Read depth'] > i- window_size)&(gl['Read depth']<= i + window_size)]
        stds.append(tmp['VAF'].std())
        means.append(tmp['VAF'].mean())
    fig,ax = plt.subplots()
    ax.plot(np.arange(min_reads, 3000, 1), stds)
    return pd.DataFrame({'Read depth':np.arange(min_reads,3000,1),'std':stds, 'mean':means})


# =============================================================================
# Import SNPs
# =============================================================================

snp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/ghent_%s_hetz_snps.vcf' % cohort, sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Ghent M1B Copy Number Plots/M1B_gene_cna.tsv', sep = '\t')

# =============================================================================
# Melt SNPS and create new columns. Keep only tissue samples
# =============================================================================

snp = pd.melt(snp, id_vars = ['CHROM','POSITION','REF','ALT','NOTES'])

snp = snp.rename(columns = {'variable':'Sample ID'})
snp['Read depth'] = snp['value'].str.split(':').str[1].astype(float)
snp['Mutant reads'] = snp['value'].str.split(':').str[0].astype(float)
snp['VAF'] = snp['Mutant reads'] / snp['Read depth']
snp['Patient ID'] = snp['Sample ID'].str.split('_').str[1]
snp['Sample type'] = snp['Sample ID'].str.split('_').str[2].str.replace('\d+', '')
# snp = snp[snp['Sample type'] != 'cfDNA']

# =============================================================================
# Select only gDNA snps called
# =============================================================================

gl = snp.copy()
gl = gl[gl['Sample type'].isin(['NL','gDNA'])]
gl = gl[gl['value'].str.contains('\*')]
gl['Called'] = True

snp = snp[~snp['Sample type'].isin(['gDNA'])]

# =============================================================================
# Calcuate moving STD
# =============================================================================

snp = snp.merge(create_depth_std(gl, window_size), on = 'Read depth', how = 'left')

# =============================================================================
# Merge germline back onto snp dataframe
# =============================================================================

gl_grouped = gl.copy()

gl = gl[['CHROM','POSITION','REF','ALT','Patient ID','Called']]
snp = snp.merge(gl, on = ['CHROM','POSITION','REF','ALT','Patient ID'])

# =============================================================================
# Get snp counts per patient
# =============================================================================

gl_grouped = gl_grouped.drop_duplicates(['Patient ID','CHROM','POSITION','REF','ALT'])
gl_grouped = gl_grouped.groupby(['Patient ID']).count()[['CHROM']]
gl_grouped.columns = ['Patient SNP count']
gl_grouped = gl_grouped.reset_index()

snp = snp.merge(gl_grouped, on = 'Patient ID')

# =============================================================================
# Keep snps with min reads > x
# =============================================================================

snp = snp[snp['Read depth'] >= min_reads]

# =============================================================================
# Merge Gene onto position
# =============================================================================

cn = cn[['GENE','CHROMOSOME','START','END']]
snp = label_genes(snp, cn)

# =============================================================================
# Keep only genes with >= min snps
# =============================================================================

sample_gene_counts = snp.groupby(['Sample ID','GENE']).count()[['CHROM']]
sample_gene_counts.columns = ['SNP count']
sample_gene_counts = sample_gene_counts.reset_index()
snp = snp.groupby(['Sample ID','GENE']).median()
snp = snp.merge(sample_gene_counts, on = ['Sample ID','GENE'])
snp = snp[snp['SNP count'] >= min_snps]

# =============================================================================
# Get number of genes evaluable per sample
# =============================================================================

gene_counts = snp.copy()
gene_counts = gene_counts.drop_duplicates(['GENE','Sample ID'])
gene_counts = gene_counts.groupby(["Sample ID"]).count()[['GENE']]
gene_counts.columns = ['Gene count']
gene_counts = gene_counts.reset_index()
snp = snp.merge(gene_counts, on = ['Sample ID'], how = 'left')

# =============================================================================
# Create p-values for each snp
# =============================================================================

snp['VAF'] = (snp['VAF'] - 0.5).abs() + 0.5
snp['p-val'] = snp.apply(lambda x: stats.norm.sf(x['VAF'], loc = x['mean'], scale = x['std']), axis = 1)

snp['Significant'] = False
snp.loc[snp['p-val'] < alpha/72, 'Significant'] = True

sig_counts_per_sample = snp.groupby(['Sample ID','Significant']).count()

sig_counts_per_sample = sig_counts_per_sample.reset_index()

pivot = pd.pivot_table(sig_counts_per_sample, index='Sample ID', columns = 'Significant', values = 'GENE', fill_value = 0)
pivot.columns = ['Copy neutral','LOH']

negative_samples = pivot[pivot['LOH'] == 0]

tumor_fraction = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tumor_fraction.columns = tumor_fraction.iloc[0]
tumor_fraction = tumor_fraction.drop(tumor_fraction.index[0])
tumor_fraction['mut_TC'] = tumor_fraction['mut_TC'].str.split('%').str[0].astype(float)

negative_samples = negative_samples.merge(tumor_fraction[['Sample ID','mut_TC']], on = 'Sample ID', how = 'left')

negative_samples.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_tc_negative_snp.tsv' % cohort, sep = '\t')

snp.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_all_snp.tsv' % cohort, sep = '\t', index = None)