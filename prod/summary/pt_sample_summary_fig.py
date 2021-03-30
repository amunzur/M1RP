# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 10:05:37 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

# =============================================================================
# Constants
# =============================================================================

age = [60,70]
psa = [35, 100, 350]
n_scoring = {'N0':8,
             'N1':7,
             'N1c':7,
             'NX':-1}

m_scoring = {'M1a':8,
             'M1a/b':8,
             'M1ab':8,
             'M1b':8,
             'M1b/c':7,
             'M1c':7}

pt_summary_colors = ['white', # NE
                     '#b3cde3','#8c96c6','#88419d', # AGE 1,2,3
                     '#fbb4b9','#f768a1','#c51b8a','#7a0177', # PSA 1-4
                     'k','lightgrey',] # YES, NO
                     #'#fd8d3c','#e6550d','#a63603'] # M1a,M1b,M1c


# =============================================================================
# Helpers
# =============================================================================

def get_TC_color(tc):
    if tc == 0:
        return 'lightgrey'
    elif tc < 0.2:
        return 'grey'
    else: 
        return '#ff6666'
    
def get_TC_cat(tc):
    if tc == 0:
        return 0
    elif tc < 0.2:
        return 1
    else: 
        return 2

# =============================================================================
# Import data from web
# =============================================================================

pt_data = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=33162023')

# =============================================================================
# Set up data for patient matrix. include age, psa at dx, bone mets, lung mets,
# other mets, and high volume disease.
# Using stampede data to split age and PSA
# =============================================================================

# YES = 7, NO = 8, UNAVAILABLE/UNEVALUABLE = -1

pt_data = pt_data[pt_data['Cohort'] == 'M1RP']

# Age uses 0-2
pt_data['Age at diagnosis (years)'] = pt_data['Age at diagnosis (years)'].apply(lambda x: 0 if x < age[0] else (1 if x >= age[0] and x < age[1] else 2))

# PSA uses 3-6
pt_data['Serum PSA at diagnosis'] = pt_data['Serum PSA at diagnosis'].apply(lambda x: 3 if x < psa[0] else (4 if x < psa[1] else (5 if x < psa[2] else 6)))

# LN mets uses TNM. If either p or c is N1, patient is N1. 
pt_data['cN'] = pt_data['cTNM'].str.split(' ').str[1].apply(lambda x: n_scoring.get(x))
pt_data['pN'] = pt_data['pTNM'].str.split(' ').str[1].apply(lambda x: n_scoring.get(x))
pt_data['N'] = np.nan
for index, row in pt_data.iterrows():
    if row['cN'] == 7 or row['pN'] == 7:
        pt_data.at[index, 'N'] = 7
    elif row['cN'] == 8 or row['pN'] == 8:
        pt_data.at[index,'N'] = 8
    else:
        pt_data.at[index,'N'] = -1

# Bone mets. Yes or no
pt_data['No. bone metastases'] = pt_data['No. bone metastases'].apply(lambda x: 8 if x == 0 else 7)

# M1 scoring
pt_data['cM'] = pt_data['cTNM'].str.split(' ').str[2].apply(lambda x: m_scoring.get(x))

# low volume (8) vs high volume (7)
pt_data['CHAARTED classification'] = pt_data['CHAARTED classification'].apply(lambda x: 7 if x == 'High volume' else 8)

# Keep desired columns
pt_data = pt_data.set_index('Patient ID')
pt_data = pt_data[['Age at diagnosis (years)','Serum PSA at diagnosis','CHAARTED classification','N','No. bone metastases','cM',]]
pt_data.columns = ['Age at Dx.','PSA at Dx.','High volume disease','Pelvic lymph nodes','Bone mets.','Other organ involvement',]

# =============================================================================
# Get sample counts
# =============================================================================

additional_samples = pd.read_excel('C:/Users/amurtha/Downloads/GSC samples submission - February 2021 - M1RP only.xlsx')
additional_samples = additional_samples[~additional_samples['Sample ID'].str.contains('cfDNA')]

s_data = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
s_data.columns = s_data.iloc[0]
s_data = s_data.drop(s_data.index[0])
s_data = s_data[s_data['Cohort'] == 'M1RP']

s_counts = s_data.groupby(['Patient ID','Sample Category']).count().reset_index()
s_counts = s_counts[['Patient ID','Sample Category','Sample ID']]
s_counts.columns = ['Patient ID','Sample Category','s_type_count']
s_counts = s_counts.pivot_table(values = 's_type_count', columns = 'Sample Category', index = 'Patient ID')
s_counts = s_counts[['PB', 'RP', 'MLN', 'cfDNA']]
s_counts = s_counts.reset_index().fillna(0)

pt_data = pt_data.merge(s_counts, on = 'Patient ID', how = 'left') 
pt_data['HasMetSample'] = pt_data['MLN'].apply(lambda x: 7 if x > 0 else 8)

# =============================================================================
# Add complete sequencign data
# =============================================================================

pt_data['WES complete'] = pt_data['Patient ID'].apply(lambda x: 8 if int(x.split('ID')[1]) in [9,20] or int(x.split('ID')[1]) >= 24 else 7)

# =============================================================================
# Sort data
# =============================================================================

pt_data = pt_data.sort_values(['WES complete','HasMetSample','High volume disease','RP', 'PB', 'MLN', 'cfDNA','Pelvic lymph nodes','Bone mets.','Other organ involvement','Age at Dx.','PSA at Dx.'], ascending = [True,True,True, False, False, False, False,True, True, False, False, False])
tmp = pt_data.copy()
pt_data = pt_data.drop(['PB', 'RP', 'MLN', 'cfDNA','HasMetSample'], axis = 1)
pt_data = pt_data.set_index('Patient ID')
pt_data = pt_data.transpose()

pt_order = pt_data.columns.tolist()

# =============================================================================
# Create matricies for sample information
# =============================================================================

fig,[ax1,ax2,ax3,ax4,ax5] = plt.subplots(nrows = 5, sharex = True, figsize = (7.5, 7), gridspec_kw={'height_ratios':[8,8,13,9,2]})

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Cohort'] == 'M1RP']
s_data = s_data.set_index('Sample ID')
tc = tc.merge(s_data[['GGG']], right_index = True, left_on = 'Sample ID', how = 'left')
tc['TC_cat'] = tc['Final tNGS_TC'].apply(lambda x: get_TC_cat(x))
tc = tc.drop_duplicates()

tc = tc.merge(s_data[['Whole-exome sequencing']], right_index = True, left_on = 'Sample ID', how = 'left')
tc['Whole-exome sequencing'] = tc['Whole-exome sequencing'].replace({'Completed':1}).fillna(0)

max_cfdna = 0

for x, pt in enumerate(pt_order):
    pb_count = 0
    rp_count = 0
    mln_count = 0
    bl_cfdna = 0
    samples = tc[tc['Patient ID'] == pt].copy().sort_values(['Whole-exome sequencing','TC_cat','GGG'], ascending = False)
    sample_tc = zip(samples['Sample ID'].tolist(), samples['Final tNGS_TC'].tolist())
    for s,s_tc in sample_tc:
        s_type = s.split('_')[2]
        c = get_TC_color(s_tc)
        ggg = s_data.at[s, 'GGG']
        wes = s_data.at[s, 'Whole-exome sequencing'] == 'Completed'
        hatch = None
        if 'PB' in s_type:
            if len(ggg) > 1:
                ggg = '5'
                wes = False
            ax2.text(x,pb_count+0.45, ggg, va = 'center', ha = 'center', fontsize = 7,zorder = 10)
            if wes:
                lw = 1
            else:
                lw = 0
            ax2.bar(x, 0.8, bottom = pb_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            pb_count += 1
        elif 'RP' in s_type:
            ggg = '' if pd.isnull(ggg) else ggg
            ax3.text(x,rp_count+0.45, ggg, va = 'center', ha = 'center', fontsize = 7, zorder = 10)
            if wes:
                lw = 1
            else:
                lw = 0
            ax3.bar(x, 0.8, bottom = rp_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            rp_count += 1
        elif 'MLN' in s_type:
            if wes:
                lw = 1
            else:
                lw = 0
            ax4.bar(x, 0.8, bottom = mln_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10, hatch = hatch)
            mln_count += 1
        elif 'cfDNA' in s_type:
            if s_data.at[s, 'Status at collection'] in ['Radical prostatectomy','Diagnostic biopsy']:
                if wes:
                    lw = 1
                else:
                    lw = 0
                ax5.bar(x, 0.8, bottom = bl_cfdna, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10)
                bl_cfdna += 1
            else:
                # print(s, s_data.at[s, 'Status at collection'])
                s=s
        else: 
            print(s)
    
ax2.set_ylabel('Prostate biopsy', fontsize = 8)
ax3.set_ylabel('Radical prostatectomy', fontsize = 8)
ax4.set_ylabel('Metastatic\nlymph node', fontsize = 8)
ax5.set_ylabel('Pre-OP\ncfDNA', fontsize = 8)

ax1.tick_params(left = False, bottom = False, labelbottom = False, labeltop = True)
ax2.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
ax3.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
ax4.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
ax5.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)

for spine in ['left','bottom']:
    ax1.spines[spine].set_visible(False)
    ax2.spines[spine].set_visible(False)
    ax3.spines[spine].set_visible(False)
    ax4.spines[spine].set_visible(False)
    ax5.spines[spine].set_visible(False)

ax2.invert_yaxis()
ax3.invert_yaxis()
ax4.invert_yaxis()
ax5.invert_yaxis()

# =============================================================================
# Plot patient summary matrix
# =============================================================================

for y, (index, row) in enumerate(pt_data.iterrows()):
    for x, col in enumerate(pt_data.columns.tolist()):
        ax1.bar(x,0.8, bottom = y, color = pt_summary_colors[int(pt_data.at[index,col])+1])
        
ax1.set_yticks(np.arange(0.4,7.4,1))
ax1.set_yticklabels(pt_data.index.tolist(), fontsize = 8)

ax1.invert_yaxis()

ax1.set_xlim(-0.55, len(pt_data.columns)-0.45)

ax1.set_xlabel('Patient')
ax1.xaxis.set_label_position('top')
ax1.set_xticks(np.arange(0,len(pt_order),1))
ax1.set_xticklabels(pt_order, fontsize = 6, rotation = 90)

fig.tight_layout()
fig.subplots_adjust(hspace = 0, top = 0.93, bottom = 0.08)

pt_summary_colors = ['white', # NE
                     '#b3cde3','#8c96c6','#88419d', # AGE 1,2,3
                     '#fbb4b9','#f768a1','#c51b8a','#7a0177', # PSA 1-4
                     'k','lightgrey', # YES, NO
                     '#fd8d3c','#e6550d','#a63603'] # M1a,M1b,M1c

handles = [
    Patch(color = 'k'),
    Patch(color = 'lightgrey'),
    Patch(color = 'white'),Patch(color = 'white'),
    Patch(color = '#b3cde3'), Patch(color = '#8c96c6'), Patch(color = '#88419d'),Patch(color = 'white'),
    Patch(color = '#fbb4b9'), Patch(color = '#f768a1'), Patch(color = '#c51b8a'), Patch(color = '#7a0177'),
    Patch(facecolor='#ff6666', edgecolor = 'k', linewidth=1.0),
    Patch(facecolor='lightgrey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='grey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='#ff6666', edgecolor = 'k', linewidth=0.0),
    ]

labels = [
    'Yes',
    'No','','',
    '<60','60-70','>70','',
    '<35','35-100','100-350','>350','Selected for WES',
    'No tumor detected','0-20% tumor','>20% tumor']

fig.legend(handles, labels, fontsize = 6, handlelength = 0.8, ncol = 4, loc = 'lower right')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/prod/summary/pt_sample_summary_fig.pdf')