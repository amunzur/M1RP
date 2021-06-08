# -*- coding: utf-8 -*-
"""
Created on Wed May 12 11:05:07 2021

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
psa = [10, 20, 350]
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
                     'k','lightgrey','grey',
                     'red','orange'] # YES, NO
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

pt_data = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx')

# =============================================================================
# Set up data for patient matrix. include age, psa at dx, bone mets, lung mets,
# other mets, and high volume disease.
# Using stampede data to split age and PSA
# =============================================================================

# YES = 7, NO = 8, UNAVAILABLE/UNEVALUABLE = -1

pt_data = pt_data[pt_data['cohort'] == 'M1RP']

# Age uses 0-2
pt_data['age_dx'] = pt_data['age_dx'].apply(lambda x: 0 if x < age[0] else (1 if x >= age[0] and x < age[1] else 2))

# PSA uses 3-6
pt_data['ipsa'] = pt_data['ipsa'].apply(lambda x: 3 if x < psa[0] else (4 if x < psa[1] else (5 if x < psa[2] else 6)))

# LN mets uses TNM. If either p or c is N1, patient is N1. 
pt_data['N'] = 7
pt_data.loc[pt_data['pelvic_ln_pos'] == 'No', 'N'] = 8
pt_data.loc[pt_data['plnd_performed'] == 0, 'N'] = -1

pt_data['crpc_status'] = pt_data['crpc_status'].apply(lambda x: 8 if x == 0 else 7)

# Bone mets. Yes or no
pt_data['bone_mets_num'] = pt_data['bone_mets_num'].apply(lambda x: 8 if x == 0 else (7 if x >= 3 else 9))

# low volume (8) vs high volume (7)
pt_data['latitude'] = pt_data['latitude'].apply(lambda x: 7 if x == 'High-risk' else 8)

# other ordan involvement
pt_data['visceral'] = pt_data['visceral'].apply(lambda x: 7 if x == 1 else 8)

pt_data['neoadj_tx'] = pt_data['neoadj_tx'].apply(lambda x: 7 if x == 1 else 8)

pt_data['Treatment intensification'] = pt_data['mHSPC_chemo_administered'].apply(lambda x: 10 if x == 1 else 8)
pt_data.loc[(pt_data['mCRPC_L1_arpi_administered'] == 1)&((pt_data['crpc_status'] == 0)|((pt_data['mCRPC_L1_arpi_start'] <= pt_data['bcr_date'])&(pt_data['mCRPC_L1_arpi_start'] <= pt_data['crpc_date']))), 'Treatment intensification'] = 11

# Keep desired columns
pt_data = pt_data.set_index('patient_id')
pt_data = pt_data[['age_dx','ipsa','latitude','N','bone_mets_num','visceral','crpc_status', 'neoadj_tx','Treatment intensification']]
pt_data.columns = ['Age at Dx.','PSA at Dx.','High risk disease','Pelvic lymph nodes','Bone mets.','Other organ involvement','CRPC','Neo-adj Tx.','Treatment intensification']

# =============================================================================
# Get sample counts
# =============================================================================

s_data = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
s_data.columns = s_data.iloc[0]
s_data = s_data.drop(s_data.index[0])
s_data = s_data[s_data['Cohort'] == 'M1RP']

s_counts = s_data.groupby(['Patient ID','Sample Category']).count().reset_index()
s_counts = s_counts[['Patient ID','Sample Category','Sample ID']]
s_counts.columns = ['Patient ID','Sample Category','s_type_count']
s_counts = s_counts.pivot_table(values = 's_type_count', columns = 'Sample Category', index = 'Patient ID')
s_counts['M'] = s_counts['MB'].fillna(0)+s_counts['MLN'].fillna(0)+s_counts['MT'].fillna(0)
s_counts = s_counts[['PB', 'RP', 'M', 'cfDNA']]
s_counts = s_counts.reset_index().fillna(0)

pt_data = pt_data.merge(s_counts, left_index = True, right_on = 'Patient ID', how = 'left') 
pt_data['HasMetSample'] = pt_data['M'].apply(lambda x: 7 if x > 0 else 8)

# =============================================================================
# Sort data
# =============================================================================

pt_data = pt_data.sort_values(['HasMetSample','CRPC','High risk disease','RP', 'PB', 'M', 'cfDNA','Pelvic lymph nodes','Bone mets.','Other organ involvement','Age at Dx.','PSA at Dx.'], ascending = [True,True, True,False, False, False, False,True, True, False, False, False])
tmp = pt_data.copy()
pt_data = pt_data.drop(['PB', 'RP', 'M', 'cfDNA','HasMetSample'], axis = 1)
pt_data = pt_data.set_index('Patient ID')
pt_data = pt_data.transpose()

pt_order = pt_data.columns.tolist()

# =============================================================================
# Create matricies for sample information
# =============================================================================

fig,[ax1,ax2,ax3,ax4,ax5] = plt.subplots(nrows = 5, sharex = True, figsize = (5, 3), gridspec_kw={'height_ratios':[9/3,8/4,13/4,9/4,2/4]})

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Cohort'] == 'M1RP']

# =============================================================================
# Merge ctDNA data onto tc data
# =============================================================================

ct = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/cfDNA timepoints.xlsx')[['Sample ID','Status']]
s_data = s_data.merge(ct, on = 'Sample ID', how = 'left')

# =============================================================================
# 
# =============================================================================

s_data = s_data.set_index('Sample ID')
tc = tc.merge(s_data[['GGG']], right_index = True, left_on = 'Sample ID', how = 'left')
tc['TC_cat'] = tc['Final tNGS_TC'].apply(lambda x: get_TC_cat(x))
tc = tc.drop_duplicates()

tc = tc.merge(s_data[['Whole-exome sequencing']], right_index = True, left_on = 'Sample ID', how = 'left')
tc['Whole-exome sequencing']=tc['Whole-exome sequencing'].replace({'Completed':1,'Selected':1}).fillna(0)

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
        wes = s_data.at[s, 'Whole-exome sequencing'] in ['Completed','Selected']
        hatch = None
        if 'PB' in s_type:
            if len(ggg) > 1:
                ggg = '5'
                wes = False
            # ax2.text(x,pb_count+0.45, ggg, va = 'center', ha = 'center', fontsize = 7,zorder = 10)
            if wes:
                lw = 0
                hatch = '\\\\\\\\\\\\\\\\\\'
            else:
                lw = 0
            ax2.bar(x, 0.8, bottom = pb_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            pb_count += 1
        elif 'RP' in s_type:
            ggg = '' if pd.isnull(ggg) else ggg
            # ax3.text(x,rp_count+0.45, ggg, va = 'center', ha = 'center', fontsize = 7, zorder = 10)
            if wes:
                lw = 0
                hatch = '\\\\\\\\\\\\\\\\\\'
            else:
                lw = 0
            ax3.bar(x, 0.8, bottom = rp_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 1, hatch = hatch)
            rp_count += 1
        elif 'MLN' in s_type or 'MT' in s_type or 'MB' in s_type:
            if wes:
                lw = 0
                hatch = '\\\\\\\\\\\\\\\\\\'
            else:
                lw = 0
            ax4.bar(x, 0.8, bottom = mln_count, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10, hatch = hatch)
            mln_count += 1
        elif 'cfDNA' in s_type:
            if s_data.at[s, 'Status'] in ['Tx naive']:
                if wes:
                    lw = 0
                    hatch = '\\\\\\\\\\\\\\\\\\'
                else:
                    lw = 0
                ax5.bar(x, 0.8, bottom = bl_cfdna, facecolor = c, edgecolor = 'k', lw = lw, zorder = 10, hatch = hatch)
                bl_cfdna += 1
            else:
                # print(s, s_data.at[s, 'Status at collection'])
                s=s
        else: 
            print(s)
    
ax2.set_ylabel('Prostate biopsy', fontsize = 6, rotation = 0, ha = 'right', va = 'center')
ax3.set_ylabel('Radical prostatectomy', fontsize = 6, rotation = 0, ha = 'right', va = 'center')
ax4.set_ylabel('Metastatic tissue', fontsize = 6, rotation = 0, ha = 'right', va = 'center')
ax5.set_ylabel('Treatment naive cfDNA', fontsize = 6, rotation = 0, ha = 'right', va = 'center')

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
        
ax1.set_yticks(np.arange(0.4,len(pt_data)+0.4,1))
ax1.set_yticklabels(pt_data.index.tolist(), fontsize = 6)
ax1.tick_params(axis = 'both', pad = 0)

ax1.invert_yaxis()

ax1.set_xlim(-0.55, len(pt_data.columns)-0.45)


ax1.xaxis.set_label_position('top')
ax1.set_xticks(np.arange(0,len(pt_order),1))
ax1.set_xticklabels(pt_order, fontsize = 6, rotation = 90)

fig.tight_layout()
fig.subplots_adjust(hspace = 0, top = 0.9, bottom = 0.08, right = 0.85)

pt_summary_colors = ['white', # NE
                     '#b3cde3','#8c96c6','#88419d', # AGE 1,2,3
                     '#fbb4b9','#f768a1','#c51b8a','#7a0177', # PSA 1-4
                     'k','lightgrey', # YES, NO
                     '#fd8d3c','#e6550d','#a63603'] # M1a,M1b,M1c

handles = [
    Patch(facecolor = 'white', lw = 0),
    Patch(facecolor = 'k', lw = 0),
    Patch(facecolor = 'lightgrey', lw = 0),
    Patch(facecolor = 'white', lw = 0),
    Patch(facecolor = 'lightgrey', lw = 0),
    Patch(facecolor = 'grey', lw = 0),
    Patch(facecolor = 'k', lw = 0),
    Patch(facecolor = 'white', lw = 0),
    Patch(facecolor = '#b3cde3', lw = 0), Patch(facecolor = '#8c96c6', lw = 0), Patch(facecolor = '#88419d', lw = 0),Patch(facecolor = 'white', lw = 0),
    Patch(facecolor = '#fbb4b9', lw = 0), Patch(facecolor = '#f768a1', lw = 0), Patch(facecolor = '#c51b8a', lw = 0), Patch(facecolor = '#7a0177', lw = 0),
    Patch(facecolor='white', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='red', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='orange', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor = 'white'),
    Patch(facecolor='#ff6666', edgecolor = 'k', hatch='\\\\\\\\\\\\\\\\\\', lw = 0),
    Patch(facecolor='lightgrey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='grey', edgecolor = 'k', linewidth=0.0),
    Patch(facecolor='#ff6666', edgecolor = 'k', linewidth=0.0),
    ]

labels = [
    "$\\bf{Clinical\ features}$",'Yes','No',
    "$\\bf{Bone\ mets}$",'0','1-2',"$\\geq$3",
    '$\\bf{Age}$','<60','60-70','>70',
    '$\\bf{PSA}$','<35','35-100','100-350','>350',
    "$\\bf{Rx\ intensification}$",'ARPI','Chemo',
    '$\\bf{Sample\ Summary}$',
    'Selected for WES',
    'No tumor detected','0-20% tumor','>20% tumor',]

fig.legend(handles, labels, fontsize = 6, handlelength = 0.8, ncol = 1, loc = 'upper right', bbox_to_anchor = (1.015,1), frameon = False)

fig.tight_layout()

fig.subplots_adjust(left = 0.2, right = 0.8, top = 0.92, bottom = 0.01, hspace = 0.01)

fig.savefig('G:/Andy Murtha/Ghent/M1RP/prod/summary/pt_sample_summary_fig_compressed.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/summary/pt_sample_summary_fig_compressed.pdf')