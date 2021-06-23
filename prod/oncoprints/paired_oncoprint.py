# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:57:10 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import string

gene_list = ['TP53', 'PTEN', 'RB1', 'SPOP', 'CHD1','FOXA1', 'BRCA2', 'ATM', 'CDK12', 'NKX3-1',
             'CLU', 'NCOA2', 'MYC']

effect_dict = {'Missense':'Missense', 'Stopgain':'Truncating', 'Non-frameshift':'Non-frameshift', 'Frameshift':'Truncating', 'Splice':'Truncating'}

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

mut_colors = {'Missense':'#79B443',
              'Truncating':'#FFC907',
              'Non-frameshift':'#BD4398'}

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# Import data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

tc = tc[~tc['Sample ID'].str.contains('cfDNA')]

# =============================================================================
# Keep mutation and CN from TC positive samples
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0]
tc_high = tc[tc['Final tNGS_TC'] > 0.2].copy()

muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
cn = cn[cn['Sample ID'].isin(tc_high['Sample ID'].tolist())]

# =============================================================================
# Keep mutations and CNA in genes
# =============================================================================

muts = keepCodingMutations(muts)

muts = muts[muts['GENE'].isin(gene_list)]
cn = cn[cn['GENE'].isin(gene_list)]
gl = gl[gl['GENE'].isin(gene_list)]

# =============================================================================
# Get sample counts per patient for CNA
# =============================================================================

high_tc_s_counts = tc_high.groupby('Patient ID').count()[['Sample ID']]
high_tc_s_counts.columns = ['pt_sample_count']

cn = cn.merge(high_tc_s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Add CNA counts. Keep CNA where count > 3. If not enough samples are present,
# remove it. Including copy neutral
# =============================================================================

cn_counts = cn.groupby(['Patient ID','GENE','Copy_num']).count()[['Log_ratio']]
cn_counts.columns = ['cn_count']
cn_counts = cn_counts.reset_index()
cn_counts = cn_counts[(cn_counts['cn_count'] >= 3)|(cn_counts['Copy_num'].isin([-2,2]))]

cn = cn.merge(cn_counts, on = ['Patient ID','GENE','Copy_num'])

# =============================================================================
# Separate into cn_m and cn_p
# =============================================================================

cn['Sample category'] = cn['Sample ID'].str.split('_').str[2].str.strip(string.digits)
cn['Sample category'] = cn['Sample category'].map({'RP':'Primary', 'MLN':'Metastatic', 'PB':'Primary', 'MT':'Metastatic', 'MB':'Metastatic', 'UCC':'UCC'})

cn = cn[cn['Sample category'] != 'UCC']

cn_m = cn[cn['Sample category'] == 'Metastatic'].copy()
cn_p = cn[cn['Sample category'] == 'Primary'].copy()

# =============================================================================
# Get denominators for both
# =============================================================================

cn = cn[(cn['Patient ID'].isin(cn_m['Patient ID'].tolist()))&(cn['Patient ID'].isin(cn_p['Patient ID'].tolist()))]

pts = cn['Patient ID'].unique().tolist()

cn_m = cn_m[cn_m['Patient ID'].isin(pts)]
cn_p = cn_p[cn_p['Patient ID'].isin(pts)]

# =============================================================================
# Get metastatic cn counts
# =============================================================================

cn_m = cn_m[['Patient ID','GENE','Copy_num']]
cn_m = cn_m[cn_m['Copy_num'] != 0]
cn_m = cn_m.drop_duplicates()

cn_m_dd = cn_m[cn_m['Copy_num'] == -2].copy()
for index, row in cn_m_dd.iterrows():
    cn_m = cn_m[(cn_m['Patient ID'] != row['Patient ID'])|(cn_m['GENE'] != row['GENE'])|(cn_m['Copy_num'] != -1)]
    
    
cn_m_amp = cn_m[cn_m['Copy_num'] == 2].copy()
for index, row in cn_m_amp.iterrows():
    cn_m = cn_m[(cn_m['Patient ID'] != row['Patient ID'])|(cn_m['GENE'] != row['GENE'])|(cn_m['Copy_num'] != 1)]


# =============================================================================
# Get primary cn counts
# =============================================================================

cn_p = cn_p[['Patient ID','GENE','Copy_num']]
cn_p = cn_p[cn_p['Copy_num'] != 0]
cn_p = cn_p.drop_duplicates()

cn_p_dd = cn_p[cn_p['Copy_num'] == -2].copy()
for index, row in cn_p_dd.iterrows():
    cn_p = cn_p[(cn_p['Patient ID'] != row['Patient ID'])|(cn_p['GENE'] != row['GENE'])|(cn_p['Copy_num'] != -1)]
    
cn_p_amp = cn_p[cn_p['Copy_num'] == 2].copy()
for index, row in cn_p_amp.iterrows():
    cn_p = cn_p[(cn_p['Patient ID'] != row['Patient ID'])|(cn_p['GENE'] != row['GENE'])|(cn_p['Copy_num'] != 1)]

# =============================================================================
# Create dataframe of genes/pts with both gain and loss called
# =============================================================================

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

cn_p = cn_p.drop_duplicates()
cn_m = cn_m.drop_duplicates()

cn_p['color'] = cn_p['Copy_num'].map(cn_color)
cn_m['color'] = cn_m['Copy_num'].map(cn_color)

cn_primary_count = cn_p[['GENE','Patient ID','color']].groupby(['GENE','Patient ID']).count()
cn_primary_count.columns = ['Count']

cn_p = cn_p.merge(cn_primary_count, left_on = ['GENE','Patient ID'], right_index = True, how = 'left')
# cn_p.loc[cn_p['Count'] == 2, 'color'] = 'white'

cn_met_count = cn_m[['GENE','Patient ID','color']].groupby(['GENE','Patient ID']).count()
cn_met_count.columns = ['Count']

cn_m = cn_m.merge(cn_met_count, left_on = ['GENE','Patient ID'], right_index = True, how = 'left')
# cn_m.loc[cn_m['Count'] == 2, 'color'] = 'white'

# =============================================================================
# Start plotting
# =============================================================================

factor = 2.5
x_dict = dict(zip(pts, np.arange(0,len(pts)*factor,factor)))

cn_m['x'] = cn_m['Patient ID'].map(x_dict) + 1
cn_p['x'] = cn_p['Patient ID'].map(x_dict)

y_dict = dict(zip(gene_list, np.arange(len(gene_list))[::-1]))

cn_m['y'] = cn_m['GENE'].map(y_dict)
cn_p['y'] = cn_p['GENE'].map(y_dict)

# =============================================================================
# Plot copy number
# =============================================================================

fig,ax = plt.subplots(figsize = (5,2.5))

for x in list(x_dict.values()):
    for gene, y in list(y_dict.items()):
        if gene != 'FOXA1':
            ax.bar(x,0.8,bottom=y, color = '#efefef', zorder = 0)
            ax.bar(x+1,0.8,bottom= y, color = '#efefef', zorder = 0)
        else:
            ax.bar(x,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')
            ax.bar(x+1,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')

ax.bar(cn_p['x'],0.8,bottom = cn_p['y'], color = cn_p['color'], zorder = 10)
ax.bar(cn_m['x'],0.8,bottom = cn_m['y'], color = cn_m['color'], zorder = 10)

cn_both = cn_p[cn_p['Count'] == 2]
for index, row in cn_both.iterrows():
    c = cn_both[(cn_both['Patient ID'] == row['Patient ID'])&(cn_both['GENE'] == row['GENE'])]['color'].tolist()
    x1 = row['x'] - 0.4
    x2 = row['x'] + 0.4
    y1 = row['y']
    y2 = row['y'] + 0.8
    p1 = Polygon(np.array([[x1,y1],[x1,y2], [x2,y1]]), color = c[0], lw = 0, zorder = 100)
    p2 = Polygon(np.array([[x2,y2],[x1,y2], [x2,y1]]), color = c[1], lw = 0, zorder = 100)
    ax.add_patch(p1)
    ax.add_patch(p2)


# =============================================================================
# Do mutations
# =============================================================================

# =============================================================================
# Get sample counts per patient for mutations
# =============================================================================

s_counts = tc.groupby('Patient ID').count()[['Sample ID']]
s_counts.columns = ['pt_sample_count']

muts = muts.merge(s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Combine germline and somatic mutations
# =============================================================================

gl = gl[gl['GENE'].isin(gene_list)]
gl = gl[(gl['NOTES'].fillna('').str.contains('ClinVar:Pathogenic'))|(gl['EFFECT'] == 'Stopgain')]

muts['Somatic'] = True

for index, row in gl.iterrows():
    pt_samples = tc[tc['Patient ID'] == row['Patient ID']]['Sample ID']
    pt_tc = tc[tc['Patient ID'] == row['Patient ID']]['Final tNGS_TC']
    muts = muts.append(pd.DataFrame({'Cohort':'M1RP', 
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
                   'pt_sample_count':len(pt_samples),
                   'Somatic':False}), ignore_index = True)

# =============================================================================
# Add mutation counts. Keep mutations where count > 3 or = num_samples
# =============================================================================

m_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
m_counts.columns = ['Mutation count']
m_counts = m_counts.reset_index()

muts = muts.merge(m_counts, on = ['Patient ID','GENE','POSITION','EFFECT'], how = 'left')

muts = muts[(muts['Mutation count'] >= 3)|((muts['pt_sample_count'] < 3)&(muts['pt_sample_count'] == muts['Mutation count']))]

# =============================================================================
# Split into metastatic and primary
# =============================================================================

tc['Sample category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
tc['Sample category'] = tc['Sample category'].map({'RP':'Primary', 'MLN':'Metastatic', 'PB':'Primary', 'MT':'Metastatic', 'MB':'Metastatic', 'UCC':'UCC'})

tc = tc[tc['Sample category'] != 'UCC']

tc_m = tc[tc['Sample category'] == 'Metastatic'].copy()
tc_p = tc[tc['Sample category'] == 'Primary'].copy()

# =============================================================================
# Get denominators for both
# =============================================================================

mut_pts = tc.groupby(['Patient ID','Sample category']).count()[['Sample ID']].reset_index()
mut_pts = mut_pts[mut_pts['Sample category'] == 'Metastatic']['Patient ID'].unique().tolist()

muts = muts[muts['Patient ID'].isin(mut_pts)]

tmp = pts + [pt for pt in mut_pts if pt not in pts]

factor = 2.5
x_dict = dict(zip(tmp, np.arange(0,len(tmp)*factor,factor)))

muts['x'] = muts['Patient ID'].map(x_dict)
muts['x'] = muts['Patient ID'].map(x_dict)

muts['y'] = muts['GENE'].map(y_dict)
muts['y'] = muts['GENE'].map(y_dict)

for pt in [pt for pt in mut_pts if pt not in pts]:
    for gene, y in list(y_dict.items()):
        x = x_dict.get(pt)
        ax.bar(x,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')
        ax.bar(x+1,0.8,bottom= y, color = '#efefef', zorder = 0, hatch = '//////', edgecolor = 'w')


# =============================================================================
# Get double/triple mutations
# =============================================================================

mut_count = muts.copy()
mut_count = mut_count.drop_duplicates(['Patient ID','GENE','EFFECT'])
mut_count = mut_count.groupby(['Patient ID','GENE']).count()[['Cohort']]
mut_count.columns = ['count']

mut_count = mut_count[mut_count['count'] >= 2]

offset = 0.15

for index, row in mut_count.iterrows():
    pt = index[0]
    gene = index[1]
    if row['count'] > 1:
        effects = muts[(muts['GENE'] == gene)&(muts['Patient ID'] == pt)]['EFFECT'].drop_duplicates().to_list()
        if row['count'] == 2:
            offset = 0.15
            for i, effect in enumerate(effects):
                if i == 0:
                    muts.loc[(muts['GENE'] == gene) & (muts['Patient ID'] == pt) & (muts['EFFECT'] == effect), 'y'] = muts['y']+offset
                elif i == 1:
                    muts.loc[(muts['GENE'] == gene) & (muts['Patient ID'] == pt) & (muts['EFFECT'] == effect), 'y'] = muts['y']-offset
        elif row['count'] == 3:
            offset= 0.25
            for i, effect in enumerate(effects):
                if i == 0:
                    muts.loc[(muts['GENE'] == gene) & (muts['Patient ID'] == pt) & (muts['EFFECT'] == effect), 'y'] = muts['y'] + offset
                elif i == 1:
                    muts.loc[(muts['GENE'] == gene) & (muts['Patient ID'] == pt) & (muts['EFFECT'] == effect), 'y'] = muts['y']
                else:
                    muts.loc[(muts['GENE'] == gene) & (muts['Patient ID'] == pt) & (muts['EFFECT'] == effect), 'y'] = muts['y'] - offset

# =============================================================================
# Get counts for metastatic and primary
# =============================================================================

mut_p = muts[muts['Sample ID'].isin(tc_p['Sample ID'])].copy()
mut_p['EFFECT'] = mut_p['EFFECT'].str.split(' ').str[0].map(effect_dict)
mut_p = mut_p[['Patient ID','GENE','EFFECT','CHROM','Somatic','x','y']]


mut_m = muts[muts['Sample ID'].isin(tc_m['Sample ID'])].copy()
mut_m['EFFECT'] = mut_m['EFFECT'].str.split(' ').str[0].map(effect_dict)
mut_m = mut_m[['Patient ID','GENE','EFFECT','CHROM','Somatic','x','y']]

# =============================================================================
# Set up color and x coordinates
# =============================================================================

mut_m['color'] = mut_m['EFFECT'].map(mut_colors)
mut_p['color'] = mut_p['EFFECT'].map(mut_colors)

# =============================================================================
# Start plotting
# =============================================================================

mut_m['x'] = mut_m['x'] + 1

# =============================================================================
# Split by somatic and germline. Then plot scatters
# =============================================================================

somatic_p = mut_p[mut_p['Somatic'] == True]
germline_p = mut_p[mut_p['Somatic'] == False]

somatic_m = mut_m[mut_m['Somatic'] == True]
germline_m = mut_m[mut_m['Somatic'] == False]

ax.scatter(somatic_p['x'], somatic_p['y']+0.4, c = somatic_p['color'], lw = 0, zorder = 1000, marker = 's', s = 8)
ax.scatter(somatic_m['x'], somatic_m['y']+0.4, c = somatic_m['color'], lw = 0, zorder = 1000, marker = 's', s = 8)

ax.scatter(germline_p['x'], germline_p['y']+0.4, c = germline_p['color'], lw = 0, zorder = 1000, marker = '*', s = 22)
ax.scatter(germline_m['x'], germline_m['y']+0.4, c = germline_m['color'], lw = 0, zorder = 1000, marker = '*', s = 22)

# =============================================================================
# Label x and y axis
# =============================================================================

ax.set_xticks(np.array(list(x_dict.values()))+0.5)
ax.set_yticks(np.array(list(y_dict.values()))+0.4)

ax.set_xticklabels(list(x_dict.keys()), rotation = 90, fontsize = 6)
ax.set_yticklabels(gene_list, fontsize = 6)

ax.set_ylim(-0.02, len(gene_list))
ax.set_xlim(-0.52, len(x_dict)*factor-0.5)

ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(left = False, bottom = False, pad = 0)

fig.tight_layout()
fig.subplots_adjust(right = 0.98, top = 0.98)

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/paired_oncoprint.pdf')
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/paired_oncoprint.png')
