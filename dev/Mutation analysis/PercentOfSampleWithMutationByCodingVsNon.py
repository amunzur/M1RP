# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:52:38 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import string
import seaborn as sns

# =============================================================================
# Constants
# =============================================================================

s_type_dict = {'PB':'Primary',
               'RP':'Primary',
               'MLN':'Metastatic',
               'MB':'Metastatic',
               'MT':'Metastatic',
               'cfDNA':'cfDNA',}

mut_cat_dict = {'Missense':'Missense',
                "5'-UTR":'Silent',
                'Intronic':'Silent',
                'Intergenic':'Silent',
                'Frameshift':'Truncation',
                'Stopgain':'Truncation',
                'Synonymous':'Silent',
                "3'-UTR":'Silent',
                'Non-frameshift':'Non-framshift',
                'Exonic':'Silent',
                'Splice':'Truncation'}

mut_color_dict = {'Missense':'green',
                  'Truncation':'yellow',
                  'Non-frameshift':'black',
                  'Silent':'grey'}

# =============================================================================
# Helpers
# =============================================================================

def fullMutationTable(muts, tc):
    df = pd.DataFrame(index = pd.MultiIndex.from_product([[],[]], names = ['Sample ID','Mutation ID']))
    for pt in muts['Patient ID'].unique().tolist():
        samples = tc[tc['Patient ID'] == pt]['Sample ID'].unique().tolist()
        pt_muts = muts[muts['Patient ID'] == pt].copy()
        pt_muts = pt_muts[['GENE','EFFECT','POSITION']].drop_duplicates()
        pt_muts = [tuple(x) for x in pt_muts.to_numpy()]
        tmp = pd.DataFrame(index = pd.MultiIndex.from_product([samples, pt_muts], names = ['Sample ID','Mutation ID']))
        df = df.append(tmp)
    
    df = df.reset_index()
    df['GENE'] = df['Mutation ID'].apply(lambda x: x[0])
    df['EFFECT'] = df['Mutation ID'].apply(lambda x: x[1])
    df['POSITION'] = df['Mutation ID'].apply(lambda x: x[2])
    df = df.drop('Mutation ID', axis = 1)
    return df;

def get_VAF(bet, x):
    if not np.isnan(x['Allele_frequency']):
        return x;
    else:
        gene = x['GENE'] if x['GENE'] != 'intergenic' else np.nan
        if float(bet.at[(gene,x['EFFECT'],x['POSITION']), x['Sample ID']].split(')')[0].split('(')[1]) > 30:
            x['Allele_frequency'] = float(bet.at[(gene,x['EFFECT'],x['POSITION']), x['Sample ID']].split('%')[0]) / 100
        else:
            x['Allele_frequency'] = np.nan
    return x;
    
# =============================================================================
# Main: 
#   Import mutations and tumor purity data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
bet = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/M1RP_betastasis_all.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Limit to TC+ samples
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0.15]
tc = tc[tc['Patient ID'] != "ID8"]
muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]

# =============================================================================
# Create full mutation dataframe. This means that for every TC positive sample
# there is a row for each mutation found in patient pt's samples 
# =============================================================================

# For syntax reasons, intergeneic mutaitons with nan gene must be changed to a value
muts['GENE'] = muts['GENE'].fillna('intergenic')

full_muts = fullMutationTable(muts, tc)

full_muts = full_muts.merge(muts[['Sample ID','GENE','EFFECT','POSITION', 'Allele_frequency']], on = ['Sample ID','GENE','EFFECT','POSITION'], how = 'left')
full_muts = full_muts.merge(tc[['Sample ID', 'Final tNGS_TC']], on = 'Sample ID', how = 'left')

# =============================================================================
# Fill in nan allele frequency values
# =============================================================================

bet = bet.set_index(['GENE','EFFECT','POSITION'])

full_muts = full_muts.apply(lambda x: get_VAF(bet, x), axis = 1)
full_muts = full_muts[~full_muts['Allele_frequency'].isnull()]

# =============================================================================
# Label sample types, patient ID. Add adjusted vaf (vaf/tc)
# =============================================================================

full_muts['Sample type'] = full_muts['Sample ID'].apply(lambda x: s_type_dict.get(x.split('_')[2].rstrip(string.digits)))
full_muts['Patient ID'] = full_muts['Sample ID'].str.split('_').str[1]

full_muts = full_muts[full_muts['Allele_frequency'] > 0]

# =============================================================================
# Group by patient, sample type, and mutation
# =============================================================================

grouped_muts = full_muts.groupby(['Patient ID','GENE','EFFECT','POSITION']).count().reset_index()
grouped_muts = grouped_muts[['Patient ID','GENE','EFFECT','POSITION','Sample ID']].rename(columns = {'Sample ID':'mut_count'})

grouped_muts['Total_samples'] = grouped_muts['Patient ID'].apply(lambda x: len(tc[(tc['Patient ID'] == x)&(tc['Cohort'] == 'M1RP')]))

grouped_muts['Percent_mut'] = grouped_muts['mut_count'] / grouped_muts['Total_samples']

grouped_muts = grouped_muts[grouped_muts['Total_samples'] > 5]

# =============================================================================
# Color muts: Coding, noncoding
# =============================================================================

grouped_muts['color'] = grouped_muts['EFFECT'].apply(lambda x: 'black' if mut_cat_dict.get(x.split(' ')[0]) == 'Silent' else 'green')

# =============================================================================
# Plot dotplot
# =============================================================================

fig,[ax1,ax2,ax3] = plt.subplots(nrows = 3, gridspec_kw = {'height_ratios':[2,1,0.8]}, figsize = (4,3), sharex=True)

ax3.scatter(grouped_muts['Percent_mut'], np.random.random(len(grouped_muts)), s = 6, alpha = 0.4, c = grouped_muts['color'])

coding = grouped_muts[grouped_muts['color'] == 'green'].copy()
noncoding = grouped_muts[grouped_muts['color'] != 'green'].copy()

#Plot KDE
# sns.kdeplot(grouped_muts['Percent_mut'], ax = ax2, cut = 0.001, bw = 0.1, shade = True)
sns.kdeplot(coding['Percent_mut'], ax = ax2, cut = 0.001, bw = 0.1, shade = False, color = 'green')
sns.kdeplot(noncoding['Percent_mut'], ax = ax2, cut = 0.001, bw = 0.1, shade = False, color = 'black')

# Plot cdf

# PDF
grouped_muts['pdf'] = grouped_muts['Percent_mut'] / sum(grouped_muts['Percent_mut'])
coding['pdf'] = coding['Percent_mut'] / sum(coding['Percent_mut'])
noncoding['pdf'] = noncoding['Percent_mut'] / sum(noncoding['Percent_mut'])

# CDF
grouped_muts = grouped_muts.sort_values('Percent_mut')
grouped_muts['cdf'] = grouped_muts['pdf'].cumsum()

coding = coding.sort_values('Percent_mut')
coding['cdf'] = coding['pdf'].cumsum()

noncoding = noncoding.sort_values('Percent_mut')
noncoding['cdf'] = noncoding['pdf'].cumsum()

ax1.plot(grouped_muts['Percent_mut'],grouped_muts['cdf'])
ax1.plot(coding['Percent_mut'],coding['cdf'], color = 'green')
ax1.plot(noncoding['Percent_mut'],noncoding['cdf'], color = 'black')
ax1.plot([0,2],[0.5,0.5],lw = 0.8, c = 'k', ls = '--')

ax1.set_xlim(0,1.03)
ax2.set_xlim(0,1.03)
ax3.set_xlim(0,1.03)

ax1.set_ylim(0,1)

ax2.tick_params(left = False, labelleft = False)
ax3.yaxis.set_visible(False)

ax2.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)

ax1.set_ylabel("CDF")
ax2.set_ylabel("KDE")

ax3.set_xlabel("Percent of patient-matched samples with a shared mutation")

ax2.get_legend().set_visible(False)


s = stats.mannwhitneyu(coding['Percent_mut'], noncoding['Percent_mut'], alternative = 'greater')
fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/FrequencyOfMutationAcrossSamplesByCvsNC.pdf")