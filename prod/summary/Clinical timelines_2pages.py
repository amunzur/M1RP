# -*- coding: utf-8 -*-
"""
Created on Wed May 12 11:11:39 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import datetime
from matplotlib.lines import Line2D


# =============================================================================
# Constants
# =============================================================================

chemo_1l_dict = {1:'Docetaxel',
              2:'Carboplatinum/Etoposide',
              3:'Carboplatinum/Cabazitaxel'}

chemo_2l_dict = {1: 'Docetaxel', 2: 'Cabazitaxel', 3: 'Cisplatinuml/Etoposide', 4:'Xofigo+xgeva', 5: 'carboplatinum'}

hormonal_dict = {1:'Abi',
                 2:'Enza',
                 3:'Abi+Apa',
                 4: 'Abi+Ipa',
                 5:'Abi->Enza', 
                 6: 'Apa'}

labels = ['CRPC dx.', 'Last FU','Death', 'cfDNA sample', 'Prostate biopsy']
handles = [Line2D([0],[0],lw = 0, marker = "$\u00A9$", color = 'k', markeredgewidth = 0),
           Line2D([0],[0],lw = 0, marker = '>', color = 'k', markeredgewidth = 0),
           Line2D([0],[0],lw = 0, marker = 'x', color = 'k', markeredgewidth = 0),
           Line2D([0],[0],lw = 0, marker = 'o', markerfacecolor = 'w', markeredgecolor = 'r'),
           Line2D([0],[0],lw = 0, marker = r"$\mathcal{Bx.}$", color = 'k', markeredgewidth = 0, )]

# =============================================================================
# Import clinical data
# =============================================================================

clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Clinical data tables/Clinical data_16.02.2021.xlsx')
clin = clin[clin['Cohort'] == 'M1RP']
clin_cpy  = clin.copy()


# =============================================================================
# Import samples
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Keep cfDNA samples
# =============================================================================

tc = tc[tc['Sample ID'].str.contains('cfDNA')]
tc['Date collected'] = tc['Sample ID'].str.split('_').str[3]

tc['Date collected'] = pd.to_datetime(tc['Date collected'], format = '%Y%b%d')


# =============================================================================
# Create plot    
# =============================================================================
    
fig,axs = plt.subplots(figsize = (8.5,11), ncols = 2, nrows = 11, sharex = False, sharey = False)

# =============================================================================
# loop over patients
# =============================================================================

for index, row in clin.iterrows():
    if index == 22:
        fig.legend(handles, labels, loc = 'lower right', fontsize = 6, ncol = 2)
        plt.tight_layout()
        fig.subplots_adjust(left = 0.02, right = 0.98, top = 0.98, bottom = 0.08)
        plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Clinical timelines/pt1_22.pdf')
        fig,axs = plt.subplots(figsize = (8.5,11), ncols = 2, nrows = 11, sharex = False, sharey = False)
    if index >= 22:
        index = index - 22
                
    
    ax = axs[index//2][index%2]
    
    # =============================================================================
    # Plot date PB    
    # =============================================================================

    date_pb = (row['Date PB'] - row['Date RP']).days / 30
    ax.scatter([date_pb],[0.2], c = 'k', s = 80, marker = r"$\mathcal{Bx}$", lw = 0)
    # ax.text(date_pb, 0.2, 'PB', ha = 'center', va = 'center', fontsize = 6)
    
    # =============================================================================
    # Plot date CRPC    
    # =============================================================================
        
    if row['CRPC'] == 1:
        date_crpc = (row['Date CRPC'] - row['Date RP']).days / 30
        ax.scatter([date_crpc],[0.20], marker = "$\u00A9$", color = 'k', lw = 0)
        # ax.text(date_crpc, 0.2015, 'CRPC Dx', ha = 'center', va = 'center', fontsize = 6)
        
    # =============================================================================
    # Plot last followup
    # =============================================================================
    
    if row['Survival Status     (1: dead, 0:alive)'] == 1:
        death_date = (row['Date of last FU'] - row['Date RP']).days / 30
        ax.scatter([death_date],[0.2], marker = 'x', c = 'k', lw = 0)
        # ax.text(death_date, 0.2015, 'Death', ha = 'center', va = 'center', fontsize = 6)
    else:
        last_fu = (row['Date of last FU'] - row['Date RP']).days / 30
        ax.scatter([last_fu],[0.2], marker = '>', c = 'k',lw = 0)
        # ax.text(last_fu, 0.2015, 'Last FU', ha = 'center', va = 'center', fontsize = 6)
        
    y = 0.2055
    
    ct_samples = tc.copy()
    ct_samples = ct_samples[ct_samples['Patient ID'] == row['Patient ID']]
    
    for ct_index, ct_row in ct_samples.iterrows():
        s_date = (ct_row['Date collected']-row['Date RP']).days / 30
        ax.scatter([s_date],[y], facecolor = 'w', edgecolor = 'r', s = 60)
        ax.text(s_date,y,'%i' % int(ct_row['Final tNGS_TC']*100), ha = 'center', va='center', fontsize = 5)
    if len(ct_samples) > 0: 
        y=y+0.005
    
    # =============================================================================
    # Plot ADT
    # =============================================================================
    if row['1° line antiandrogen therapy (1:yes, 2:no)'] == 1:
        adt_start = (row['Start of 1° line antiandrogen therapy'] - row['Date RP']).days / 30
        if str(row['Stop of 1° line antiandrogen therapy']).strip() == 'continued':
            adt_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([adt_start,adt_end],[y,y], marker = None, lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.scatter([adt_start],[y], marker = '|', lw = 0, c = 'k')
            ax.scatter([adt_end],[y], marker = '>', lw = 0, c = 'k')
            ax.text((adt_start+adt_end)/2, y+0.0015, 'ADT', ha = 'center', fontsize = 6)
        else:    
            adt_end = (row['Stop of 1° line antiandrogen therapy'] - row['Date RP']).days / 30
            ax.plot([adt_start,adt_end],[y,y], marker = '|', lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.text((adt_start+adt_end)/2, y+0.0015, 'ADT', ha = 'center', fontsize = 6)
        y = y + 0.005
    
    # =============================================================================
    # Plot 1st line cytotoxic therapy
    # =============================================================================
    
    chemo_1l = row['1° line cytotoxic therapy (0:no, 1: Docetaxel, 2:Carboplatinum/Etoposide 3:Carboplatinum/Cabazitaxel)']
    if not pd.isnull(chemo_1l) and chemo_1l not in ['-',0]:
        chemo_1l_start = (row['Start of 1° line cytotoxic therapy'] - row['Date RP']).days / 30
        if isinstance(row['Stop of 1° line cytotoxic therapy'],datetime.datetime):
            chemo_1l_end = (row['Stop of 1° line cytotoxic therapy'] - row['Date RP']).days / 30
            ax.plot([chemo_1l_start,chemo_1l_end],[y,y], marker = '|', lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.text((chemo_1l_start+chemo_1l_end)/2, y+0.0015, chemo_1l_dict.get(chemo_1l), ha = 'center', fontsize = 6)
        else:
            chemo_1l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_1l_start,chemo_1l_end],[y,y], marker = None, lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.scatter([chemo_1l_start],[y], marker = '|', lw = 0, c = 'k')
            ax.scatter([chemo_1l_end],[y], marker = '>', lw = 0, c = 'k')
            ax.text((chemo_1l_start+chemo_1l_end)/2, y+0.0015, chemo_1l_dict.get(chemo_1l), ha = 'center', fontsize = 6)
        y = y + 0.005
    
    # =============================================================================
    # Plot 2nd line hormonal therapy    
    # =============================================================================
    
    if row['2° line hormonal therapy (0:no, 1: abiraterone, 2: enzalutamide, 3: abiraterone + apalutamide, 4: abiraterone + ipatasertib 5:abiraterone followed by enzalutamide), 6: apalutamide'] not in ['-',0]:
        hor_2l_start = (row['Start 2° line hormonal therapy'] - row['Date RP']).days / 30
        hor_2l = row['2° line hormonal therapy (0:no, 1: abiraterone, 2: enzalutamide, 3: abiraterone + apalutamide, 4: abiraterone + ipatasertib 5:abiraterone followed by enzalutamide), 6: apalutamide']
        if isinstance(row['Stop 2° line hormonal therapy'],datetime.datetime):
            hor_2l_end = (row['Stop 2° line hormonal therapy'] - row['Date RP']).days / 30
            ax.plot([hor_2l_start,hor_2l_end],[y,y], marker = '|', lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.text((hor_2l_start+hor_2l_end)/2, y+0.0015, hormonal_dict.get(hor_2l), ha = 'center', fontsize = 6)
        else:
            hor_2l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([hor_2l_start,hor_2l_end],[y,y], marker = None, lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.scatter([hor_2l_start],[y], marker = '|', lw = 0, c = 'k')
            ax.scatter([hor_2l_end],[y], marker = '>', lw = 0, c = 'k')
            ax.text((hor_2l_start+hor_2l_end)/2, y+0.0015, hormonal_dict.get(hor_2l), ha = 'center', fontsize = 6)
        y = y + 0.005
            
    # =============================================================================
    # Plot second line Cytotoxic therapy
    # =============================================================================
    chemo_2l = row['2° line cytotoxic therapy  (0:no, 1: Docetaxel, 2: Cabazitaxel 3: Cisplatinuml/Etoposide 4:Xofigo+xgeva, 5: carboplatinum)']       
    if isinstance(chemo_2l, int) and chemo_2l > 0:
        chemo_2l_start = (row['Start 2° line cytotoxic therapy'] - row['Date RP']).days / 30
        
        if isinstance(row['Stop 2° line cytotoxic therapy'],datetime.datetime):
            chemo_2l_end = (row['Stop 2° line cytotoxic therapy'] - row['Date RP']).days / 30
            ax.plot([chemo_2l_start,chemo_2l_end],[y,y], marker = '|', lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.text((chemo_2l_start+chemo_2l_end)/2, y+0.0015, chemo_2l_dict.get(chemo_2l), ha = 'center', fontsize = 6)
        else:
            chemo_2l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_2l_start,chemo_2l_end],[y,y], marker = None, lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.scatter([chemo_2l_start],[y], marker = '|', lw = 0, c = 'k')
            ax.scatter([chemo_2l_end],[y], marker = '>', lw = 0, c = 'k')
            ax.text((chemo_2l_start+chemo_2l_end)/2, y+0.0015, chemo_2l_dict.get(chemo_2l), ha = 'center', fontsize = 6)
        y = y + 0.005
            
            
    # =============================================================================
    # Plot third line Cytotoxic therapy
    # =============================================================================
            
    if row['3° line cytotoxic treatment'] > 0:
        chemo_3l_start = (row['Start 3° line treatment '] - row['Date RP']).days / 30
        chemo_3l = row['Treatment3']
        if isinstance(row['Stop 3° line cytotoxic treatment'],datetime.datetime):
            chemo_3l_end = (row['Stop 3° line cytotoxic treatment'] - row['Date RP']).days / 30
            ax.plot([chemo_3l_start,chemo_3l_end],[y,y], marker = '|', lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.text((chemo_3l_start+chemo_3l_end)/2,y+0.0015, chemo_3l, ha = 'center', fontsize = 6)
        else:
            chemo_3l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_3l_start,chemo_3l_end],[y,y], marker = None, lw = 0.5, c = 'k', markeredgewidth = 0)
            ax.scatter([chemo_3l_start],[y], marker = '|', lw = 0, c = 'k')
            ax.scatter([chemo_3l_end],[y], marker = '>', lw = 0, c = 'k')
            ax.text((chemo_3l_start+chemo_3l_end)/2, y+0.0015, chemo_3l, ha = 'center', fontsize = 6)
        y = y + 0.005
            
    # =============================================================================
    # Plot fourth line Cytotoxic therapy
    # =============================================================================
            
    if row['4° line cytotoxic treatment'] > 0:
        chemo_4l_start = (row['Start 4° line cytotoxic traetment'] - row['Date RP']).days / 30
        chemo_4l = row['Treatment2']
        if isinstance(row['Stop 4° line cytotoxic treatment'],datetime.datetime):
            chemo_4l_end = (row['Stop 4° line cytotoxic treatment'] - row['Date RP']).days / 30
            ax.plot([chemo_4l_start,chemo_4l_end],[y,y], marker = '|', lw = 0.5, color = 'k', markeredgewidth = 0)
            ax.text((chemo_4l_start+chemo_4l_end)/2, y+0.0015, chemo_4l, ha = 'center', fontsize = 6)
        else:
            chemo_4l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_4l_start,chemo_4l_end],[y,y], marker = None, lw = 0.5, markeredgewidth = 0)
            ax.scatter([chemo_4l_start],[y], marker = '|', c = 'k', lw = 0)
            ax.scatter([chemo_4l_end],[y], marker = '>', c = 'k', lw = 0)
            ax.text((chemo_4l_start+chemo_4l_end)/2, y+0.0015, chemo_4l, ha = 'center', fontsize = 6)
        y = y + 0.005
    '''
    # =============================================================================
    # Plot BCR
    # =============================================================================
    
    if row['BCR failure (0:no; 1:yes)'] == 1:
        date_bcr = (row['Date of BCR FU'] - row['Date RP']).days / 30
        ax.scatter([date_bcr],[y])
        ax.text(date_bcr, y+0.0015, 'BCR\nPSA=%s' % row['PSA at BCR'], ha = 'center', va = 'center', fontsize = 6)
        y = y + 0.005
            
    # =============================================================================
    # Plot radiological progression event    
    # =============================================================================
    
    if row['Radiologic progression event (0: no radiologic progression; 1: radiologic new lesions)'] == 1:
        date_rpe = (row['Date radiologic progression event'] - row['Date RP']).days / 30
        ax.scatter([date_rpe],[y], c = 'k', )
        ax.text(date_rpe, y+0.0015, 'RPE\nPSA=%s' % row['PSA at radiologic progression '], ha = 'center', va = 'center', fontsize = 6)
        y = y + 0.005
    '''
        
    ax.plot([0,0],[0,y], ls = 'dashed',lw=0.5, alpha = 0.5, zorder = 0)
    ax.set_ylim(0.197, y)
    if index//2 == 10: 
        ax.set_xlabel('Months since RP', fontsize = 6)

    
    # ax.set_yticklabels()
    
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis = 'y', left = False, labelleft = False)
    ax.set_title(row['Patient ID'], fontsize = 6)
    ax.tick_params(axis = 'x', labelsize = 6)
    
# =============================================================================
# 
# =============================================================================

axs[10][1].set_visible(False)

fig.legend(handles, labels, loc = 'lower right', fontsize = 6, ncol = 2)

plt.tight_layout()
fig.subplots_adjust(left = 0.02, right = 0.98, top = 0.98, bottom = 0.05)
    
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Clinical timelines/pt23_43.pdf')