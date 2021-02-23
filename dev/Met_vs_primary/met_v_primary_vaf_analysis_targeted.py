# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:15:34 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

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

# =============================================================================
# Merge TC and VAF back onto mutations
# =============================================================================

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
full_muts['aVAF'] = full_muts['Allele_frequency'] / full_muts['Final tNGS_TC']

# =============================================================================
# Group by patient, sample type, and mutation
# =============================================================================

grouped_muts = full_muts.groupby(['Patient ID','Sample type','GENE','EFFECT','POSITION']).mean().reset_index()

# =============================================================================
# Pivot the table to include primary, metastatic, and cfDNA as titles
# =============================================================================

pivot_muts = grouped_muts.pivot(index = ['Patient ID','GENE','POSITION','EFFECT'], columns = 'Sample type', values = 'aVAF')

# =============================================================================
# Adjust color of genes of interest
# =============================================================================

pivot_muts = pivot_muts.reset_index()
pivot_muts['mut_category'] = pivot_muts['EFFECT'].apply(lambda x: mut_cat_dict.get(x.split(' ')[0]))
pivot_muts['Color'] = pivot_muts['mut_category'].apply(lambda x: mut_color_dict.get(x))

# =============================================================================
# Create Truncal score = mP*mM
# =============================================================================

pivot_muts['Truncal_score'] = pivot_muts['Primary'] * pivot_muts['Metastatic']

# =============================================================================
# Plot figure
# =============================================================================

fig,ax = plt.subplots(figsize = (2.5,2.5))

p_muts = pivot_muts.copy().dropna()

ax.scatter(p_muts['Primary'],p_muts['Metastatic'], c = p_muts['Color'], s = 8, alpha = 0.6)

ax.set_xlabel("Mean corrected VAF in primary samples", fontsize = 8)
ax.set_ylabel("Mean corrected VAF in metastatic samples", fontsize = 8)

ax.tick_params(labelsize = 6)

ax.set_xlim(-0.03, 1)
ax.set_ylim(-0.03, 1)

import scipy.stats as stats

s = stats.linregress(p_muts['Primary'], p_muts['Metastatic'])
fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/meanCorrectedVaf_PrimaryVsMet.pdf")

# =============================================================================
# Plot truncal score by mutation type
# =============================================================================

fig,ax = plt.subplots()

boxplots = [p_muts[p_muts['mut_category'] == x]['Truncal_score'].tolist() for x in p_muts['mut_category'].unique().tolist()]
labels = [x for x in p_muts['mut_category'].unique().tolist()]

ax.boxplot(boxplots, labels = labels)

ax.set_xticklabels(['%s\nn=%i' % (s,len(n)) for (s,n) in zip(labels,boxplots)])
ax.set_ylabel('Truncal score')

fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/truncal_score_byMutEffect.pdf")

'''
# =============================================================================
# Plot truncal score by gene for CODING mutation only
# =============================================================================

coding = p_muts[p_muts['mut_category'] != 'Silent'].copy()
boxplots = [coding[coding['GENE'] == x]['Truncal_score'].tolist() for x in coding['GENE'].unique().tolist() if len(coding[coding['GENE'] == x]['Truncal_score'].tolist()) > 2]
labels = [x for x in coding['GENE'].unique().tolist() if len(coding[coding['GENE'] == x]['Truncal_score'].tolist()) > 2]

fig,ax = plt.subplots(figsize = (4,1))

ax.boxplot(boxplots, labels = labels)
ax.tick_params(labelsize = 8)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
'''