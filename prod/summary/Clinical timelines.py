# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:54:12 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import datetime


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
# loop over patients
# =============================================================================

for index, row in clin.iterrows():
    # if row['Patient ID'] != 'ID5':
    #     continue;
    # =============================================================================
    # Create plot    
    # =============================================================================
    
    fig,ax = plt.subplots()
    
    # =============================================================================
    # Plot date PB    
    # =============================================================================

    date_pb = (row['Date PB'] - row['Date RP']).days / 30
    ax.scatter([date_pb],[0.2])
    ax.text(date_pb, 0.2015, 'PB', ha = 'center')
    
    # =============================================================================
    # Plot date CRPC    
    # =============================================================================
        
    if row['CRPC'] == 1:
        date_crpc = (row['Date CRPC'] - row['Date RP']).days / 30
        ax.scatter([date_crpc],[0.20])
        ax.text(date_crpc, 0.2015, 'CRPC Dx', ha = 'center', va = 'center')
        
    # =============================================================================
    # Plot last followup
    # =============================================================================
    
    if row['Survival Status     (1: dead, 0:alive)'] == 1:
        death_date = (row['Date of last FU'] - row['Date RP']).days / 30
        ax.scatter([death_date],[0.2], marker = 'x', c = 'k')
        ax.text(death_date, 0.2015, 'Death', ha = 'center', va = 'center')
    else:
        last_fu = (row['Date of last FU'] - row['Date RP']).days / 30
        ax.scatter([last_fu],[0.2], marker = '>', c = 'k')
        ax.text(last_fu, 0.2015, 'Last FU', ha = 'center', va = 'center')
        
    y = 0.205
    
    ct_samples = tc.copy()
    ct_samples = ct_samples[ct_samples['Patient ID'] == row['Patient ID']]
    
    for ct_index, ct_row in ct_samples.iterrows():
        s_date = (ct_row['Date collected']-row['Date RP']).days / 30
        ax.scatter([s_date],[y], c = 'red')
        ax.text(s_date,y+0.0015,'%i%%' % int(ct_row['Final tNGS_TC']*100), ha = 'center', va='center')
    if len(ct_samples) > 0: 
        y=y+0.005
    
    # =============================================================================
    # Plot ADT
    # =============================================================================
    if row['1° line antiandrogen therapy (1:yes, 2:no)'] == 1:
        adt_start = (row['Start of 1° line antiandrogen therapy'] - row['Date RP']).days / 30
        if str(row['Stop of 1° line antiandrogen therapy']).strip() == 'continued':
            adt_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([adt_start,adt_end],[y,y], marker = None)
            ax.scatter([adt_start],[y], marker = '|')
            ax.scatter([adt_end],[y], marker = '>')
            ax.text((adt_start+adt_end)/2, y+0.0015, 'ADT', ha = 'center')
        else:    
            adt_end = (row['Stop of 1° line antiandrogen therapy'] - row['Date RP']).days / 30
            ax.plot([adt_start,adt_end],[y,y], marker = '|')
            ax.text((adt_start+adt_end)/2, y+0.0015, 'ADT', ha = 'center')
        y = y + 0.005
    
    # =============================================================================
    # Plot 1st line cytotoxic therapy
    # =============================================================================
    
    chemo_1l = row['1° line cytotoxic therapy (0:no, 1: Docetaxel, 2:Carboplatinum/Etoposide 3:Carboplatinum/Cabazitaxel)']
    if not pd.isnull(chemo_1l) and chemo_1l not in ['-',0]:
        chemo_1l_start = (row['Start of 1° line cytotoxic therapy'] - row['Date RP']).days / 30
        if isinstance(row['Stop of 1° line cytotoxic therapy'],datetime.datetime):
            chemo_1l_end = (row['Stop of 1° line cytotoxic therapy'] - row['Date RP']).days / 30
            ax.plot([chemo_1l_start,chemo_1l_end],[y,y], marker = '|')
            ax.text((chemo_1l_start+chemo_1l_end)/2, y+0.0015, chemo_1l_dict.get(chemo_1l), ha = 'center')
        else:
            chemo_1l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_1l_start,chemo_1l_end],[y,y], marker = None)
            ax.scatter([chemo_1l_start],[y], marker = '|')
            ax.scatter([chemo_1l_end],[y], marker = '>')
            ax.text((chemo_1l_start+chemo_1l_end)/2, y+0.0015, chemo_1l_dict.get(chemo_1l), ha = 'center')
        y = y + 0.005
    
    # =============================================================================
    # Plot 2nd line hormonal therapy    
    # =============================================================================
    
    if row['2° line hormonal therapy (0:no, 1: abiraterone, 2: enzalutamide, 3: abiraterone + apalutamide, 4: abiraterone + ipatasertib 5:abiraterone followed by enzalutamide), 6: apalutamide'] not in ['-',0]:
        hor_2l_start = (row['Start 2° line hormonal therapy'] - row['Date RP']).days / 30
        hor_2l = row['2° line hormonal therapy (0:no, 1: abiraterone, 2: enzalutamide, 3: abiraterone + apalutamide, 4: abiraterone + ipatasertib 5:abiraterone followed by enzalutamide), 6: apalutamide']
        if isinstance(row['Stop 2° line hormonal therapy'],datetime.datetime):
            hor_2l_end = (row['Stop 2° line hormonal therapy'] - row['Date RP']).days / 30
            ax.plot([hor_2l_start,hor_2l_end],[y,y], marker = '|')
            ax.text((hor_2l_start+hor_2l_end)/2, y+0.0015, hormonal_dict.get(hor_2l), ha = 'center')
        else:
            hor_2l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([hor_2l_start,hor_2l_end],[y,y], marker = None)
            ax.scatter([hor_2l_start],[y], marker = '|')
            ax.scatter([hor_2l_end],[y], marker = '>')
            ax.text((hor_2l_start+hor_2l_end)/2, y+0.0015, hormonal_dict.get(hor_2l), ha = 'center')
        y = y + 0.005
            
    # =============================================================================
    # Plot second line Cytotoxic therapy
    # =============================================================================
    chemo_2l = row['2° line cytotoxic therapy  (0:no, 1: Docetaxel, 2: Cabazitaxel 3: Cisplatinuml/Etoposide 4:Xofigo+xgeva, 5: carboplatinum)']       
    if isinstance(chemo_2l, int) and chemo_2l > 0:
        chemo_2l_start = (row['Start 2° line cytotoxic therapy'] - row['Date RP']).days / 30
        
        if isinstance(row['Stop 2° line cytotoxic therapy'],datetime.datetime):
            chemo_2l_end = (row['Stop 2° line cytotoxic therapy'] - row['Date RP']).days / 30
            ax.plot([chemo_2l_start,chemo_2l_end],[y,y], marker = '|')
            ax.text((chemo_2l_start+chemo_2l_end)/2, y+0.0015, chemo_2l_dict.get(chemo_2l), ha = 'center')
        else:
            chemo_2l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_2l_start,chemo_2l_end],[y,y], marker = None)
            ax.scatter([chemo_2l_start],[y], marker = '|')
            ax.scatter([chemo_2l_end],[y], marker = '>')
            ax.text((chemo_2l_start+chemo_2l_end)/2, y+0.0015, chemo_2l_dict.get(chemo_2l), ha = 'center')
        y = y + 0.005
            
            
    # =============================================================================
    # Plot third line Cytotoxic therapy
    # =============================================================================
            
    if row['3° line cytotoxic treatment'] > 0:
        chemo_3l_start = (row['Start 3° line treatment '] - row['Date RP']).days / 30
        chemo_3l = row['Treatment3']
        if isinstance(row['Stop 3° line cytotoxic treatment'],datetime.datetime):
            chemo_3l_end = (row['Stop 3° line cytotoxic treatment'] - row['Date RP']).days / 30
            ax.plot([chemo_3l_start,chemo_3l_end],[y,y], marker = '|')
            ax.text((chemo_3l_start+chemo_3l_end)/2,y+0.0015, chemo_3l, ha = 'center')
        else:
            chemo_3l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_3l_start,chemo_3l_end],[y,y], marker = None)
            ax.scatter([chemo_3l_start],[y], marker = '|')
            ax.scatter([chemo_3l_end],[y], marker = '>')
            ax.text((chemo_3l_start+chemo_3l_end)/2, y+0.0015, chemo_3l, ha = 'center')
        y = y + 0.005
            
    # =============================================================================
    # Plot fourth line Cytotoxic therapy
    # =============================================================================
            
    if row['4° line cytotoxic treatment'] > 0:
        chemo_4l_start = (row['Start 4° line cytotoxic traetment'] - row['Date RP']).days / 30
        chemo_4l = row['Treatment2']
        if isinstance(row['Stop 4° line cytotoxic treatment'],datetime.datetime):
            chemo_4l_end = (row['Stop 4° line cytotoxic treatment'] - row['Date RP']).days / 30
            ax.plot([chemo_4l_start,chemo_4l_end],[y,y], marker = '|')
            ax.text((chemo_4l_start+chemo_4l_end)/2, y+0.0015, chemo_4l, ha = 'center')
        else:
            chemo_4l_end = (row['Date of last FU'] - row['Date RP']).days / 30
            ax.plot([chemo_4l_start,chemo_4l_end],[y,y], marker = None)
            ax.scatter([chemo_4l_start],[y], marker = '|')
            ax.scatter([chemo_4l_end],[y], marker = '>')
            ax.text((chemo_4l_start+chemo_4l_end)/2, y+0.0015, chemo_4l, ha = 'center')
        y = y + 0.005
    
    # =============================================================================
    # Plot BCR
    # =============================================================================
    
    if row['BCR failure (0:no; 1:yes)'] == 1:
        date_bcr = (row['Date of BCR FU'] - row['Date RP']).days / 30
        ax.scatter([date_bcr],[y])
        ax.text(date_bcr, y+0.0015, 'BCR\nPSA=%s' % row['PSA at BCR'], ha = 'center', va = 'center')
        y = y + 0.005
            
    # =============================================================================
    # Plot radiological progression event    
    # =============================================================================
    
    if row['Radiologic progression event (0: no radiologic progression; 1: radiologic new lesions)'] == 1:
        date_rpe = (row['Date radiologic progression event'] - row['Date RP']).days / 30
        ax.scatter([date_rpe],[y])
        ax.text(date_rpe, y+0.0015, 'RPE\nPSA=%s' % row['PSA at radiologic progression '], ha = 'center', va = 'center')
        y = y + 0.005
        
    ax.plot([0,0],[0,y], ls = 'dashed',lw=0.5, alpha = 0.5)
    ax.set_ylim(0.199, y)
    ax.set_xlabel('Months since RP')
    
    # ax.set_yticklabels()
    
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis = 'y', left = False, labelleft = False)
    ax.set_title(row['Patient ID'])
    
    plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Clinical timelines/%s.pdf' % row['Patient ID'])