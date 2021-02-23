# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:06:14 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import string

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

gleason_dict = {'5+5':5,
                '5+4':5,
                '4+5':5,
                '4+4':4, 
                '4+3':3,
                '4+3+mini5':3,
                '3+4':2, 
                '3+3':1, 
    }



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

pt_data = pd.read_csv('https://docs.google.com/spreadsheets/d/1QHxb-zXtmIihVSvFXne6loL6hXT3aKbLIynk3N2uKzA/export?format=csv&gid=1432498324')

pt_data = pt_data[~pt_data['Patient'].isin(['Median','Q1','Q3','IQR'])]
pt_data = pt_data[~pt_data['Patient'].isnull()]

# =============================================================================
# Set up data for patient matrix. include age, psa at dx, bone mets, lung mets,
# other mets, and high volume disease.
# Using stampede data to split age and PSA
# =============================================================================

# YES = 7, NO = 8, UNAVAILABLE/UNEVALUABLE = -1

# Age uses 0-2
pt_data['Age at PCa diagnosis (years)'] = pt_data['Age at PCa diagnosis (years)'].apply(lambda x: 0 if x < age[0] else (1 if x >= age[0] and x < age[1] else 2))

# PSA uses 3-6
pt_data['PSA at diagnosis'] = pt_data['PSA at diagnosis'].apply(lambda x: 3 if x < psa[0] else (4 if x < psa[1] else (5 if x < psa[2] else 6)))

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
pt_data['Bone metastases (yes/no)'] = pt_data['Bone metastases (yes/no)'].apply(lambda x: 8 if x == 'no' else 7)

# M1 scoring
pt_data['cM'] = pt_data['cTNM'].str.split(' ').str[2].apply(lambda x: m_scoring.get(x))

# low volume (8) vs high volume (7)

# Keep desired columns
pt_data = pt_data[['Patient','Age at PCa diagnosis (years)','PSA at diagnosis','N','Bone metastases (yes/no)',]]
pt_data.columns = ['Patient_ID','Age at Dx.','PSA at Dx.','Pelvic lymph nodes','Bone mets.']

# =============================================================================
# Get sample counts
# =============================================================================

s_data = pd.read_csv('https://docs.google.com/spreadsheets/d/1QHxb-zXtmIihVSvFXne6loL6hXT3aKbLIynk3N2uKzA/export?format=csv&gid=1917404413')

s_data['Sample Category'] = s_data['Oncoprint_ID'].str.split('_').str[-1].str.strip(string.digits)

s_counts = s_data.groupby(['Patient_ID','Sample Category']).count().reset_index()
s_counts = s_counts[['Patient_ID','Sample Category','Oncoprint_ID']]
s_counts.columns = ['Patient_ID','Sample Category','s_type_count']
s_counts = s_counts.pivot_table(values = 's_type_count', columns = 'Sample Category', index = 'Patient_ID')
s_counts = s_counts.reset_index().fillna(0)

pt_data = pt_data.merge(s_counts, on = 'Patient_ID', how = 'left')

# =============================================================================
# Add complete sequencign data
# =============================================================================

pt_data['Sequencing complete'] =  pt_data['Patient_ID'].apply(lambda x: 8 if x in ['ID8','ID5'] else 7)
pt_data['gDNA complete'] = pt_data['Patient_ID'].apply(lambda x: 8 if x in ['ID8'] else 7)

# =============================================================================
# Sort data
# =============================================================================

pt_data = pt_data.sort_values(['gDNA complete','P', 'ML', 'MLN', 'MT','Pelvic lymph nodes','Bone mets.','Age at Dx.','PSA at Dx.'], ascending = [True, False, False, False, False,True, True, False, False])
tmp = pt_data.copy()
pt_data = pt_data.drop(['P', 'ML', 'MLN', 'MT'], axis = 1)
pt_data = pt_data.set_index('Patient_ID')
pt_data = pt_data.transpose()

pt_order = pt_data.columns.tolist()

# =============================================================================
# Set up GGG score in s_data
# =============================================================================

s_data['GGG'] = s_data['Gleason score'].replace(gleason_dict)

# =============================================================================
# Create matricies for sample information
# =============================================================================

# fig,[ax1,ax2,ax3,ax4,ax5] = plt.subplots(nrows = 5, sharex = True, figsize = (7.5, 7), gridspec_kw={'height_ratios':[6,7,8,3,1]})
fig,[ax1,ax2,ax3,ax4] = plt.subplots(nrows = 4, sharex = True, figsize = (7.5, 7), gridspec_kw={'height_ratios':[6,7,8,4]})

s_data['WES completed?'] = s_data['WES completed?'].replace({'Completed':1}).fillna(0)
s_data['TC_cat'] = s_data['Tumor fraction (tNGS)'].apply(lambda x: get_TC_cat(x))

max_cfdna = 0

for x, pt in enumerate(pt_order):
    p_count = 0
    ml_count = 0
    m_count = 0
    samples = s_data[s_data['Patient_ID'] == pt].copy().sort_values(['WES completed?','TC_cat','GGG'], ascending = False)
    sample_tc = zip(samples['Oncoprint_ID'].tolist(), samples['Tumor fraction (tNGS)'].tolist())
    for s,s_tc in sample_tc:
        s_type = s.split('_')[-1].strip(string.digits)
        c = get_TC_color(s_tc)
        ggg = s_data.set_index('Oncoprint_ID').at[s, 'GGG']
        wes = s_data.set_index('Oncoprint_ID').at[s, 'WES completed?'] == 1
        if s in ['ID5_P4']:
            hatch = '///////'
            c = 'white'
        else:
            hatch = None
        if 'P' == s_type:
            ax2.text(x,p_count+0.45, ggg, va = 'center', ha = 'center', fontsize = 7,zorder = 10)
            if wes:
                lw = 1
            else:
                lw = 0
            ax2.bar(x, 0.8, bottom = p_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            p_count += 1
        elif 'ML' == s_type:
            if wes:
                lw = 1
            else:
                lw = 0
            ax3.bar(x, 0.8, bottom = ml_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            ml_count += 1
        elif 'MLN' in s_type or 'MT' == s_type:
            if wes:
                lw = 1
            else:
                lw = 0
            ax4.bar(x, 0.8, bottom = m_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10, hatch = hatch)
            m_count += 1
        elif 'MT' == s_type:
            if wes:
                lw = 1
            else:
                lw = 0
            ax4.bar(x, 0.8, bottom = m_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10, hatch = hatch)
            m_count += 1
        else: 
            print("Missing %s" % s)
    
ax2.set_ylabel('Prostate sample', fontsize = 8)
ax3.set_ylabel('Metastatic lung', fontsize = 8)
ax4.set_ylabel('Metastatic\nsamples', fontsize = 8)
# ax5.set_ylabel('Metastatic\ntissue', fontsize = 8)

ax1.tick_params(left = False, bottom = False, labelbottom = False, labeltop = True)
ax2.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
ax3.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
ax4.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
# ax5.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)

for spine in ['left','bottom']:
    ax1.spines[spine].set_visible(False)
    ax2.spines[spine].set_visible(False)
    ax3.spines[spine].set_visible(False)
    ax4.spines[spine].set_visible(False)
    # ax5.spines[spine].set_visible(False)

ax2.invert_yaxis()
ax3.invert_yaxis()
ax4.invert_yaxis()
# ax5.invert_yaxis()

# =============================================================================
# Plot patient summary matrix
# =============================================================================

for y, (index, row) in enumerate(pt_data.iterrows()):
    for x, col in enumerate(pt_data.columns.tolist()):
        ax1.bar(x,0.8, bottom = y, color = pt_summary_colors[int(pt_data.at[index,col])+1])
        
ax1.set_yticks(np.arange(0.4,len(pt_data),1))
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
    Patch(color = '#b3cde3'), Patch(color = '#8c96c6'), Patch(color = '#88419d'),
    Patch(color = '#fbb4b9'), Patch(color = '#f768a1'), Patch(color = '#c51b8a'), Patch(color = '#7a0177'),
    Patch(facecolor='#ff6666', edgecolor = 'k', linewidth=1.0),
    Patch(facecolor='lightgrey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='grey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='#ff6666', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='white', edgecolor = 'k', linewidth=0, hatch = '///////')
    ]

labels = [
    'Yes',
    'No',
    '<60','60-70','>70',
    '<35','35-100','100-350','>350','Selected for WES',
    'No tumor detected','0-20% tumor','>20% tumor','Undergoing targeted sequencing']

fig.legend(handles, labels, fontsize = 6, handlelength = 0.8, ncol = 3, loc = 'lower right')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/pt_sample_summary_fig_lum.pdf')