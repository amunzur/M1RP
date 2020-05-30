# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:02:16 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import math


# =============================================================================
# Constants
# =============================================================================

min_reads = 100
alpha = 0.001
min_snps = 3
min_loh = 4
window_size = 25 # In each direction
cna_lr_threshold = 0.3
top_fraction = 0.75 # Take the top 75% of genes with allelic imbalance
sig_gene_fraction = 0.5 # 1/2 of the snps in a gene must be significant to count it
cohort = 'm1b'

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
    # fig,ax = plt.subplots()
    # ax.plot(np.arange(min_reads, 3000, 1), stds)
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
if cohort != 'abienza':
    snp['Patient ID'] = snp['Sample ID'].str.split('_').str[1]
    snp['Sample type'] = snp['Sample ID'].str.split('_').str[2].str.replace('\d+', '')
else:
    snp['Patient ID'] = snp['Sample ID'].str.split('-').str[0]+'_'+snp['Sample ID'].str.split('-').str[1]
    snp['Sample type'] = 'cfDNA'
    snp.loc[snp['Sample ID'].str.contains('WBC'), 'Sample type'] = 'gDNA'
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

snp['combined p-val'] = snp.groupby(['Sample ID','GENE'])['p-val'].transform(lambda x: [stats.combine_pvalues(x, method = 'fisher')[1]]*len(x))

snp['snp_significant'] = None
snp.loc[snp['p-val'] < alpha/(snp['SNP count']), 'snp_significant'] = True

snp['Significant'] = False
snp.loc[snp['combined p-val'] < alpha/(snp['SNP count']*snp['Gene count']), 'Significant'] = True

# =============================================================================
# Calculate number of significant SNPs per gene. Keep only genes that have >= min_snp
# =============================================================================

sig_counts_per_gene = snp.copy()
sig_counts_per_gene = snp[['Sample ID','GENE','snp_significant']].groupby(['Sample ID','GENE']).count().reset_index()
sig_counts_per_gene.columns = ['Sample ID','GENE','sig_count_gene']

snp = snp.merge(sig_counts_per_gene, on = ['Sample ID','GENE'])
snp.loc[snp['sig_count_gene'] < sig_gene_fraction * snp['SNP count'].clip(lower = 2), 'Significant'] = False

# =============================================================================
# Calculate number of regions of LOH
# =============================================================================

snp = snp.groupby(['Sample ID','GENE']).median().reset_index()

sig_counts_per_sample = snp.groupby(['Sample ID','Significant']).count()
sig_counts_per_sample = sig_counts_per_sample.reset_index()

sig_counts_per_sample = sig_counts_per_sample[sig_counts_per_sample['Significant']==True]
sig_counts_per_sample = sig_counts_per_sample[['Sample ID', 'GENE']]
sig_counts_per_sample.columns = ['Sample ID', 'Num LOH']

snp = snp.merge(sig_counts_per_sample, on = 'Sample ID', how = 'left')


# =============================================================================
# Merge chromosome back into dataframe
# =============================================================================

cn.columns = ['GENE','CHROM','START','END']
snp = snp.merge(cn[['GENE','CHROM']], on  = 'GENE')

# =============================================================================
# Keep LOH and eliminate sex chromosomes
# =============================================================================

tc = snp.copy()
tc = tc[tc['Significant'] == True]
tc = tc[tc['Num LOH'] >= min_loh]
tc = tc[~tc['CHROM'].isin(['chrX','chrY'])]

# =============================================================================
# Merge copy number data onto cohort. Keep only gene < cna_lr_threshold
# =============================================================================

# cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/%s_cna_melted.tsv' % cohort, sep = '\t')
# cn = cn[['Sample ID','GENE','Log_ratio']]

# tc = tc.merge(cn, on = ['Sample ID','GENE'], how = 'left')
# tc = tc[tc['Log_ratio'] < cna_lr_threshold]

# =============================================================================
# Keep upper half of LOH genes
# =============================================================================

tc = tc.sort_values('combined p-val', ascending = True)
count_dict = dict(zip(tc['Sample ID'].drop_duplicates(),[0]*len(tc['Sample ID'].drop_duplicates())))
                 
tc['num'] = None
for index, row in tc.iterrows():
    tc.at[index,'num'] = count_dict.get(row['Sample ID'])
    count_dict[row['Sample ID']] = count_dict.get(row['Sample ID']) + 1
    
tc['num'] = tc['num'] / tc['Num LOH']
tc = tc[tc['num'] < top_fraction]

# =============================================================================
# Calculate TC based on SNPs
# =============================================================================

tc = tc.groupby('Sample ID').mean()[['VAF', 'Gene count']]

tc['snp_TC'] = 2 - (1 / tc['VAF'])
tc = tc.reset_index()

tc = snp[['Sample ID']].merge(tc, on = 'Sample ID', how = 'left')
tc['snp_TC'] = tc['snp_TC'].fillna(0)
tc['Gene count'] = tc['Gene count'].fillna(0)
tc = tc.drop_duplicates('Sample ID')[['Sample ID', 'VAF','snp_TC','Gene count']]
tc.columns = ['Sample ID','median LOH VAF','snp_TC','LOH gene count']

# =============================================================================
# Merge mut_TC calls onto snp_TC
# =============================================================================
if cohort in ['m1b','m1rp']:
    mut_TC = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
    
    mut_TC.columns = mut_TC.iloc[0]
    mut_TC = mut_TC.drop(mut_TC.index[0])
    mut_TC['mut_TC'] = mut_TC['mut_TC'].str.split('%').str[0].astype(float) / 100
    mut_TC = mut_TC[['Sample ID','mut_TC']]
elif cohort == 'lum':
    mut_TC = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/lm_tc.tsv', sep = '\t')
    mut_TC['mut_TC'] = mut_TC['mut_TC'] / 100
    mut_TC = mut_TC[['Sample ID','mut_TC']]
else:
    mut_TC = pd.read_csv('C:/Users/amurtha/Dropbox/2017 - Abi-enza second manuscript/Data freeze/ctdna_fractions.tsv', sep = '\t', header = None, names = ['Sample ID','mut_TC'])
    mut_TC['mut_TC'] = mut_TC['mut_TC'] / 100

tc = tc.merge(mut_TC, on = 'Sample ID')

pos_tc = tc.copy()
pos_tc = pos_tc[(pos_tc['snp_TC'] != 0)&(pos_tc['mut_TC'] != 0)]

pos_tc['residual'] = pos_tc['snp_TC'] - pos_tc['mut_TC']

tc['Sample type'] = tc['Sample ID'].str.split('_').str[2]
tc['order'] = 1
if cohort != 'abienza':
    tc.loc[tc['Sample type'] == 'cfDNA', 'order'] = 2

tc = tc.sort_values(['order','snp_TC','mut_TC','Sample ID'], ascending = True)
tc = tc.merge(snp[['Sample ID','Gene count']].drop_duplicates(), on = ['Sample ID'], how = 'left')
tc['Copy neutral gene count'] = tc['Gene count'] - tc['LOH gene count']
tc = tc[['Sample ID','Sample type','median LOH VAF','LOH gene count','Copy neutral gene count','snp_TC','mut_TC']]

tc.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_tc_snp.tsv' % cohort, sep = '\t', index = None)

snp = snp.sort_values(['Sample ID','Significant'], ascending = [True, False])
snp.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_all_snp.tsv' % cohort, sep = '\t', index = None)
# =============================================================================
# Plot mut_TC vs snp_TC
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(tc['mut_TC'],tc['snp_TC'], s = 10, c = 'k', alpha = 0.5)
ax.plot([0,1],[0,1])
ax.set_ylabel('snp_TC')
ax.set_ylim(0,1)
ax.set_xlabel('mut_TC')
ax.set_xlim(0,1)

ax.set_title(cohort)

linregress = stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC'])
slope = str(round(linregress[0],2))
r = str(round(linregress[2],2))
std_xy = round(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum()/len(pos_tc['residual'])),2)
std_xy = str(std_xy)

ax.plot([0,1],[linregress[1],linregress[0]+linregress[1]], linestyle = 'dashed')

ax.text(x = 0.01, y = 0.95, s = 'm: %s, r: %s\nstd from y=x: %s\nn=%i' % (slope, r, std_xy, len(tc)))


print(stats.linregress(pos_tc['mut_TC'], pos_tc['snp_TC']))
print(math.sqrt(pos_tc['residual'].apply(lambda x: x**2).sum() / len(pos_tc['residual'])))

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/SNP_analysis/%s_mutsnpTC_scatter.pdf' % cohort)
