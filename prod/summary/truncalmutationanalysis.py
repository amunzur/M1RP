# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:27:34 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# =============================================================================
# Import mutation data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')

wes_bet = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/betastasis/june2021_wes_mutations_betastasis.tsv', sep = '\t')
wes_samples = wes_bet.columns[7:].tolist()

del wes_bet

# =============================================================================
# Create dictionary of dataframes with cluster data
# =============================================================================

secondary_tumors = ['M1RP_ID30_RP3','M1RP_ID30_UCC','M1RP_ID19_cfDNA_2017Jan13']
wes_samples = [s for s in wes_samples if s not in secondary_tumors]

path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/pyclone_wes'
pts = [f.name for f in os.scandir(path) if f.is_dir() and not 'old' in f.name and not f in ['ID9','ID20']]

# =============================================================================
# Get per patient cluster and loci files
# =============================================================================

pts_to_remove = []

clusters = {}
loci = {}
for pt in pts:
    if pt == 'ID4':
        pt = pt
    if not os.path.exists('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/pyclone_wes/%s/pyclone_output/tables/cluster.tsv' % pt):
        pts_to_remove.append(pt)
        continue;
    c = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/pyclone_wes/%s/pyclone_output/tables/cluster.tsv' % pt, sep = '\t')
    c = c[(c['size'] > 1)&(c['mean'] > 0.05)]
    pivot = c.pivot(index = 'cluster_id', columns = 'sample_id', values = 'mean')
    pivot = pivot[[col for col in pivot.columns.tolist() if col in wes_samples]]
    pivot = pivot.dropna(how='any')
    clusters[pt] = pivot
    
    l = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/pyclone_wes/%s/pyclone_output/tables/loci_filtered.tsv' % pt, sep = '\t')
    l = l[l['cluster_id'].isin(pivot.index)]
    loci[pt] = l

del c, pivot

pts = [pt for pt in pts if pt not in pts_to_remove]


# =============================================================================
# Merge cluster information onto mutations
# =============================================================================

all_loci = pd.concat(list(loci.values()), ignore_index = True)

muts['EFFECT'] = muts['EFFECT'].str.split(' ').str[0]
muts['mutation_id'] = muts['CHROM']+'_'+muts['GENE']+'_'+muts['POSITION'].astype(str)+'_'+muts['REF']+'_'+muts['ALT']+'_'+muts['EFFECT']
muts = muts.rename(columns = {'Sample ID':'sample_id'})
muts = muts[muts['sample_id'].isin(wes_samples)]
muts = muts[muts['Patient ID'].isin(pts)]

muts = muts.merge(all_loci[['mutation_id','sample_id','cluster_id']], on = ['mutation_id','sample_id'], how = 'left')

# =============================================================================
# Create columns for truncal cluster
# =============================================================================

muts['truncal'] = 1
muts.loc[muts['cluster_id'].isna(), 'truncal'] = 0

muts = muts[muts['Patient ID'] != 'ID8']


truncal_counts = muts.copy()
truncal_counts = truncal_counts.drop_duplicates('mutation_id')
truncal_counts = truncal_counts.groupby(['GENE','truncal']).count().reset_index()