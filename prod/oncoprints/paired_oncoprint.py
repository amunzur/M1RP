# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:57:10 2021

@author: amurtha
"""
# =============================================================================
# Import Packages and set parameters
# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

#What genes do you want to show?
gene_list = ['TP53', 'PTEN', 'RB1', 'SPOP', 'FOXA1', 'BRCA2', 'ATM', 'CDK12', 'NKX3-1',
             'CLU', 'NCOA2', 'MYC', 'PIK3CA']

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
mut = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')

# =============================================================================
# Keep coding mutations
# =============================================================================

mut = mut[mut['EFFECT'].str.split(' ').str[0].isin(['Missense','Frameshift','Stopgain','Splice','Non-frameshift'])]

# =============================================================================
# Keep TC > 20.
# =============================================================================

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

tc = tc[tc['Final tNGS_TC']>0.15]

# =============================================================================
# Keep patients with >= 2 metastatic samples. 
# =============================================================================

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

samples = samples[samples['Sample ID'].isin(tc['Sample ID'].tolist())]

samples['Sample Category'] = samples['Sample Category'].replace(dict.fromkeys(['MLN','cfDNA','MB','MT'], 'Met'))
samples['Sample Category'] = samples['Sample Category'].replace(dict.fromkeys(['RP','PB'], 'Primary'))

met_count = samples[samples['Sample Category'] == 'Met']
met_count = met_count[['Patient ID']]
met_count['met_count'] = 1
met_count = met_count.groupby('Patient ID').count().reset_index()

samples = samples.merge(met_count, on = 'Patient ID', how = 'left')
samples['met_count'] = samples['met_count'].fillna(0)

pri_count = samples[samples['Sample Category'] == 'Primary']
pri_count = pri_count[['Patient ID']]
pri_count['primary_count'] = 1
pri_count = pri_count.groupby('Patient ID').count().reset_index()

samples = samples.merge(pri_count, on = 'Patient ID', how = 'left')
samples['primary_count'] = samples['primary_count'].fillna(0)

samples = samples[(samples['met_count'] >= 1)&(samples['primary_count'] >= 1)]
patients = samples['Patient ID'].unique().tolist()

# =============================================================================
# Keep patients in samples
# =============================================================================

cn = cn[cn['Sample ID'].isin(samples['Sample ID'].tolist())]
mut = mut[mut['Sample ID'].isin(samples['Sample ID'].tolist())]

# =============================================================================
# Keep copy number data present in > 2 samples
# =============================================================================

cn = cn[cn['GENE'].isin(gene_list)]
mut = mut[mut['GENE'].isin(gene_list)]

cn['Copy_num'] = cn['Copy_num'].replace(dict.fromkeys([-1, -2], -1))
cn['Copy_num'] = cn['Copy_num'].replace(dict.fromkeys([1, 2], 1))

cn_count = cn[['Patient ID','GENE','Copy_num']].copy()
cn_count = cn_count[cn_count['Copy_num'] != 0]
cn_count['Count'] = 1
cn_count = cn_count.groupby(['Patient ID','GENE','Copy_num']).count().reset_index()
cn_count = cn_count[cn_count['Count'] >= 1]

cn = cn.merge(cn_count, on = ['Patient ID','GENE','Copy_num'])

# =============================================================================
# Combine germline and somatic mutations
# =============================================================================

gl = gl[gl['GENE'].isin(gene_list)]
gl = gl[(gl['NOTES'].fillna('').str.contains('ClinVar:Pathogenic'))|(gl['EFFECT'] == 'Stopgain')]

mut['Somatic'] = True

for index, row in gl.iterrows():
    pt_samples = tc[tc['Patient ID'] == row['Patient ID']]['Sample ID']
    pt_tc = tc[tc['Patient ID'] == row['Patient ID']]['Final tNGS_TC']
    mut = mut.append(pd.DataFrame({'Cohort':'M1RP', 
                   'Patient ID':row['Patient ID'], 
                   'Sample ID':pt_samples,
                   'CHROM':row['CHROM'], 
                   'POSITION':row['POSITION'], 
                   'REF':row['REF'], 
                   'ALT':row['ALT'],
                   'GENE':row['GENE'], 
                   'EFFECT':row['EFFECT'], 
                   'Allele_frequency':row['Allele_frequency'], 
                   'Read_depth':row['Read_depth'], 
                   'NOTES':row['NOTES'],
                   'Independent':True, 
                   'Final tNGS_TC':pt_tc, 
                   'Clonal':True,
                   'Somatic':False}), ignore_index = True)

# =============================================================================
# Assign Y values
# =============================================================================

ys = dict(zip(gene_list[::-1], np.arange(0,len(gene_list),1)))
cn['y'] = cn['GENE'].map(ys)
mut['y'] = mut['GENE'].map(ys)

# =============================================================================
# Assign color to muts and cn
# =============================================================================

cn['color'] = cn['Copy_num'].map({-1:'#9CC5E9',0:'#efefef',1:'#F59496'})
mut['color'] = mut['EFFECT'].str.split(' ').str[0].map({'Missense':'#79B443',
                                                        'Frameshift':'#FFC907',
                                                        'Stopgain':'#FFC907',
                                                        'Splice':'#BD4398',
                                                        'Non-frameshift':'#bbbbbb'})


# =============================================================================
# order
# =============================================================================

order = mut.copy()
order = order.drop_duplicates(['Patient ID','GENE','POSITION','EFFECT'])
order['val'] = 2
order['val'] = order.apply(lambda x: x['val']**ys.get(x['GENE']), axis = 1)

order = order.groupby(['Patient ID']).sum().reset_index()
order = order.sort_values(['val', 'Patient ID'], ascending = [False, True])

patients = order['Patient ID'].tolist() + [p for p in patients if p not in order['Patient ID'].tolist()]

# =============================================================================
# Add X coordinates
# =============================================================================

factor = 2.5
xs = dict(zip(patients, np.arange(0,len(patients)*factor,factor)))

cn['x'] = cn['Patient ID'].map(xs)
mut['x'] = mut['Patient ID'].map(xs)

# =============================================================================
# Find double mutations
# =============================================================================

mut_count = mut[['Patient ID','GENE','EFFECT']].drop_duplicates().copy()
mut_count = mut_count.groupby(['Patient ID','GENE']).count()
mut_count.columns = ['count']

offset = 0.15

for index, row in mut_count.iterrows():
    pt = index[0]
    gene = index[1]
    if row['count'] > 1:
        effects = mut[(mut['GENE'] == gene)&(mut['Patient ID'] == pt)]['EFFECT'].drop_duplicates().to_list()
        if row['count'] == 2:
            offset = 0.15
            for i, effect in enumerate(effects):
                if i == 0:
                    mut.loc[(mut['GENE'] == gene) & (mut['Patient ID'] == pt) & (mut['EFFECT'] == effect), 'y'] = mut['y']+offset
                elif i == 1:
                    mut.loc[(mut['GENE'] == gene) & (mut['Patient ID'] == pt) & (mut['EFFECT'] == effect), 'y'] = mut['y']-offset
        elif row['count'] == 3:
            offset= 0.25
            for i, effect in enumerate(effects):
                if i == 0:
                    mut.loc[(mut['GENE'] == gene) & (mut['Patient ID'] == pt) & (mut['EFFECT'] == effect), 'y'] = mut['y']+offset
                elif i == 1:
                    mut.loc[(mut['GENE'] == gene) & (mut['Patient ID'] == pt) & (mut['EFFECT'] == effect), 'y'] = mut['y']
                else:
                    mut.loc[(mut['GENE'] == gene) & (mut['Patient ID'] == pt) & (mut['EFFECT'] == effect), 'y'] = mut['y'] - offset


# =============================================================================
# Split by primary v metastatic
# =============================================================================

cn = cn.merge(samples[['Sample ID','Sample Category']], on = 'Sample ID')
mut = mut.merge(samples[['Sample ID','Sample Category']], on = 'Sample ID')

cn_primary = cn[cn['Sample Category'] == 'Primary'].copy()
cn_met = cn[cn['Sample Category'] == 'Met'].copy()

mut_primary = mut[mut['Sample Category'] == 'Primary'].copy()
mut_met = mut[mut['Sample Category'] == 'Met'].copy()

del cn, tc, samples, met_count, pri_count

# =============================================================================
# 
# =============================================================================

cn_primary = cn_primary[['Patient ID','GENE','Copy_num','y','color','x']]
cn_primary = cn_primary.drop_duplicates()

cn_met = cn_met[['Patient ID','GENE','Copy_num','y','color','x']]
cn_met = cn_met.drop_duplicates()

cn_met['x'] = cn_met['x'] + 1
mut_met['x'] = mut_met['x'] + 1

# =============================================================================
# Create dataframe of genes/pts with both gain and loss called
# =============================================================================

cn_primary_count = cn_primary[['GENE','Patient ID','color']].groupby(['GENE','Patient ID']).count()
cn_primary_count.columns = ['Count']

cn_primary = cn_primary.merge(cn_primary_count, left_on = ['GENE','Patient ID'], right_index = True, how = 'left')
cn_primary.loc[cn_primary['Count'] == 2, 'color'] = 'white'


# =============================================================================
# Create plot
# =============================================================================

fig,ax = plt.subplots(figsize = (5.625,2.75))

# =============================================================================
# plot cna data
# =============================================================================

for x in list(xs.values()):
    for gene, y in list(ys.items()):
        if gene != 'FOXA1':
            ax.bar(x,0.8,bottom= y, color = '#efefef', zorder = 0)
            ax.bar(x+1,0.8,bottom= y, color = '#efefef', zorder = 0)
        else:
            ax.bar(x,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')
            ax.bar(x+1,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')

ax.bar(cn_primary['x'],0.8,bottom = cn_primary['y'], color = cn_primary['color'], zorder = 10)
ax.bar(cn_met['x'],0.8,bottom = cn_met['y'], color = cn_met['color'], zorder = 10)

cn_both = cn_primary[cn_primary['Count'] == 2]
for index, row in cn_both.iterrows():
    x1 = row['x'] - 0.4
    x2 = row['x'] + 0.4
    y1 = row['y']
    y2 = row['y'] + 0.8
    p1 = Polygon(np.array([[x1,y1],[x1,y2], [x2,y1]]), color = '#9CC5E9', lw = 0, zorder = 100)
    p2 = Polygon(np.array([[x2,y2],[x1,y2], [x2,y1]]), color = '#F59496', lw = 0, zorder = 100)
    ax.add_patch(p1)
    ax.add_patch(p2)

# =============================================================================
# Plot mutation data
# =============================================================================

somatic_p = mut_primary[mut_primary['Somatic'] == True]
germline_p = mut_primary[mut_primary['Somatic'] == False]

somatic_m = mut_met[mut_met['Somatic'] == True]
germline_m = mut_met[mut_met['Somatic'] == False]

ax.scatter(somatic_p['x'], somatic_p['y']+0.4, c = somatic_p['color'], lw = 0, zorder = 1000, marker = 's', s = 8)
ax.scatter(somatic_m['x'], somatic_m['y']+0.4, c = somatic_m['color'], lw = 0, zorder = 1000, marker = 's', s = 8)

ax.scatter(germline_p['x'], germline_p['y']+0.4, c = germline_p['color'], lw = 0, zorder = 1000, marker = '*', s = 22)
ax.scatter(germline_m['x'], germline_m['y']+0.4, c = germline_m['color'], lw = 0, zorder = 1000, marker = '*', s = 22)

# =============================================================================
# Label x and y axis
# =============================================================================

ax.set_xticks(np.array(list(xs.values()))+0.5)
ax.set_yticks(np.array(list(ys.values()))+0.4)

ax.set_xticklabels(patients, rotation = 90, fontsize = 6)
ax.set_yticklabels(gene_list[::-1], fontsize = 6)

ax.set_ylim(-0.02, len(gene_list))
ax.set_xlim(-0.52, len(patients)*factor-0.5)

ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(left = False, bottom = False, pad = 0)

# ax.invert_yaxis()

fig.tight_layout()
plt.savefig('G:/Andy Murtha/Ghent/M1RP/prod/summary/paired_oncoprint.pdf')
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Oncoprints/paired_oncoprint.pdf')
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Oncoprints/paired_oncoprint.png')