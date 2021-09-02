# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:43:20 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# Constants
# =============================================================================

genes = ['TP53','PTEN','RB1','FOXA1','SPOP','AR', 'BRCA2','CDK12','ATM']
sample_type_color = {0:'green',
                     1:'yellow',
                     2:'purple',
                     3:'red',
                     np.nan:'white'}

tc_color = { 1:'#ffcccc',
             3:'#ff6666',
             5:'#ff0000'}

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

muts_color = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907'}

def get_tc_color(x):
    if x == 0:
        return '#ffcccc';
    elif x < 0.2:
        return '#ff6666'
    else:
        return '#ff0000'


n = 2

s_order = ['ID3','ID43','ID10','ID14','ID15','ID19','ID21','ID23','ID33','ID38','ID40', 'ID1', 'ID4', 'ID5','ID6', 'ID7', 'ID9', 'ID16', 'ID17', 'ID18', 'ID20', 'ID22', 'ID24', 'ID25', 'ID27', 'ID28', 'ID31', 'ID32', 'ID34', 'ID35', 'ID36', 'ID37', 'ID39', 'ID42', 'ID2', 'ID8', 'ID11', 'ID12', 'ID13', 'ID26', 'ID29', 'ID30', 'ID41']

# =============================================================================
# Import copy number and mutation data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')
clin = pd.read_excel("C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx")

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

tc = tc[~tc['Sample ID'].str.contains('cfDNA')]

# =============================================================================
# Clinicl
# =============================================================================

clin['study_id'] = clin['study_id'].str.split('-').str[1]
clin = clin[clin['cohort'] == 'M1RP']
clin = clin[['study_id','latitude']]
clin['color'] = clin['latitude'].map({'Low-risk':'#E6E7E8','High-risk':'k'})

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

muts = muts[muts['GENE'].isin(genes)]
cn = cn[cn['GENE'].isin(genes)]
gl = gl[gl['GENE'].isin(genes)]

# =============================================================================
# Add germline mutations
# =============================================================================

gl = gl[(gl['NOTES'].str.contains('ClinVar:Pathogenic')) | (gl['EFFECT'].isin(['Frameshift','Stopgain','Splice']))]
gl = gl[~gl['NOTES'].str.contains('Benign')]

muts['Somatic'] = True

for index, row in gl.iterrows():
    muts = muts.append(pd.DataFrame({'Cohort':'M1RP', 
                   'Patient ID':row['Patient ID'], 
                   'Sample ID':tc['Sample ID'],
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
                   'Final tNGS_TC':tc['Final tNGS_TC'], 
                   'Clonal':True,
                   'Somatic':False}), ignore_index = True)

# =============================================================================
# Get sample counts per patient for mutations
# =============================================================================

s_counts = tc.groupby('Patient ID').count()[['Sample ID']]
s_counts.columns = ['pt_sample_count']

muts = muts.merge(s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Get sample counts per patient for CNA
# =============================================================================

high_tc_s_counts = tc_high.groupby('Patient ID').count()[['Sample ID']]
high_tc_s_counts.columns = ['pt_sample_count']

cn = cn.merge(high_tc_s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Add mutation counts. Keep mutations where count > 3 or = num_samples
# =============================================================================

m_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
m_counts.columns = ['Mutation count']
m_counts = m_counts.reset_index()

muts = muts.merge(m_counts, on = ['Patient ID','GENE','POSITION','EFFECT'], how = 'left')

muts = muts[(muts['Mutation count'] >= 3)|((muts['pt_sample_count'] < 3)&(muts['pt_sample_count'] == muts['Mutation count']))]

# =============================================================================
# Add CNA counts. Keep CNA where count > 3. If not enough samples are present,
# remove it. Including copy neutral
# =============================================================================

cn_counts = cn.groupby(['Patient ID','GENE','Copy_num']).count()[['Log_ratio']]
cn_counts.columns = ['cn_count']
cn_counts = cn_counts.reset_index()

cn = cn.merge(cn_counts, on = ['Patient ID','GENE','Copy_num'])

## If using 3 as limit for AMP/DeepDel, then use this code so no shallow dels are lost
'''
for index, row in cn.iterrows():
    pt = row['Patient ID']
    g = row['GENE']
    if row['Copy_num'] == -2 and row['cn_count'] < 3:
        cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == -1), 'cn_count'] = cn['cn_count']+row['cn_count']
    elif row['Copy_num'] == 2 and row['cn_count'] < 3:
        cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == 1), 'cn_count'] = cn['cn_count'] + row['cn_count']
cn = cn[cn['cn_count'] >= 3]
'''
# Keep amps and deep dels and gains/shallow losses if count >= 3
cn = cn[(cn['Copy_num'].isin([-2,2]))|(cn['cn_count'] >= 3)]


# =============================================================================
# Assign x and y coordinates to mutations and CN events 
# =============================================================================

x_dict = dict(zip(s_order,np.arange(len(s_order))))
y_dict = dict(zip(genes[::-1], np.arange(len(genes))))

cn['y'] = cn['GENE'].map(y_dict)
cn['x'] = cn['Patient ID'].map(x_dict)

muts['y'] = muts['GENE'].map(y_dict)
muts['x'] = muts['Patient ID'].map(x_dict)

clin['x'] = clin['study_id'].map(x_dict)
clin['y'] = len(y_dict)

# =============================================================================
# Drop dupbliate mutations
# =============================================================================

muts = muts.drop_duplicates(['Patient ID','GENE','EFFECT'])

# =============================================================================
# Find multiple mutations
# =============================================================================

m_count_pt = muts.groupby(['Patient ID','GENE']).count()[['Sample ID']]
m_count_pt.columns = ['gene_mut_count']
m_count_pt = m_count_pt.reset_index()

muts = muts.merge(m_count_pt, on = ['Patient ID','GENE'], how = 'left')

offset = 0.15

for index, row in m_count_pt.iterrows():
    if row['gene_mut_count'] > 1:
        gene = row['GENE']
        pt = row['Patient ID']
        effects = muts[(muts['GENE'] == gene)&(muts['Patient ID'] == pt)]['EFFECT'].drop_duplicates().to_list()
        for i, effect in enumerate(effects):
            if i == 0:
                # muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']+offset
                muts.loc[(muts['GENE'] == gene) & (muts['EFFECT'] == effect) & (muts['Patient ID'] == pt), 'y'] = muts['y']+offset
            else:
                # muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']-offset
                muts.loc[(muts['GENE'] == gene) & (muts['EFFECT'] == effect) & (muts['Patient ID'] == pt), 'y'] = muts['y']-offset

# =============================================================================
# Add color
# =============================================================================

cn['color'] = cn['Copy_num'].map(cn_color)
muts['color'] = muts['EFFECT'].str.split(' ').str[0].map(muts_color)

# =============================================================================
# Plot
# =============================================================================

fig, ax = plt.subplots(figsize = (5,2.25))

for i in np.arange(len(x_dict)):
    for j in np.arange(len(y_dict)):
        ax.bar(i,0.8,bottom = j, color = cn_color.get(0), hatch = '\\\\\\\\\\', zorder = 0, edgecolor = 'w')

# =============================================================================
# CN 
# =============================================================================

cn_0 = cn[cn['Copy_num'] == 0].copy()
ax.bar(cn_0['x'], 0.8, bottom = cn_0['y'], color = cn_0['color'], zorder = 10)

cn_1 = cn[cn['Copy_num'].isin([-1,1])].copy()

ax.bar(cn_1['x'], 0.8, bottom = cn_1['y'], color = cn_1['color'], zorder = 50)

cn_2 = cn[cn['Copy_num'].isin([-2,2])]

ax.bar(cn_2['x'], 0.8, bottom = cn_2['y'], color = cn_2['color'], zorder = 100)

# =============================================================================
# Clinical risk criteria
# =============================================================================

ax.bar(clin['x'], 0.8, bottom = clin['y'], color = clin['color'])

# =============================================================================
# Muts
# =============================================================================

muts_s = muts[muts['Somatic'] == True]
muts_gl = muts[muts['Somatic'] == False]

ax.scatter(muts_s['x'], muts_s['y']+0.4, c = muts_s['color'], marker = 's', lw = 0, zorder = 1000, s = 17)
ax.scatter(muts_gl['x'], muts_gl['y']+0.4, c = muts_gl['color'], marker = '*', lw = 0, zorder = 1000)

# =============================================================================
# Aethetics
# =============================================================================

ax.set_xticks(np.arange(len(x_dict)))
ax.set_yticks(np.arange(0.4,len(y_dict)+1,1))

ax.set_xticklabels(s_order, rotation = 90, fontsize = 6)
ax.set_yticklabels(genes[::-1]+['High risk'], fontsize = 6)

ax.set_xlim(-0.5, len(x_dict)-0.4)
ax.set_ylim(-0.1, len(y_dict)+1)

ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(left = False, bottom = False, pad = 0)

plt.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Oncoprints/aggregate_oncoprint.pdf')