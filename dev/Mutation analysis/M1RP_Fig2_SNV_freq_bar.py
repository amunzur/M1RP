# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 10:31:19 2020

@author: Sarah
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy import stats
from matplotlib import gridspec
import seaborn as sns


include_wes = False
sort_by_median = False


targeted = pd.read_csv('Y:/users/sng/m1rp_mutations.tsv', sep='\t')
ids = pd.read_csv('Y:/users/sng/m1rp_patient_sample_IDs.tsv', sep='\t')

# Exclude ID8 hypermutant
ids = ids[ids['Patient ID'] != 'ID8']
targeted = targeted[targeted['Patient ID'] != 'ID8']

#Filter for only primary tissue samples
targeted = targeted.loc[(targeted['Sample ID'].str.contains('_PB'))|(targeted['Sample ID'].str.contains('_RP'))]
ids = ids.loc[(ids['Sample ID'].str.contains('_PB'))|(ids['Sample ID'].str.contains('_RP'))]

# Number of samples per patient
sample_count = ids.groupby(['Patient ID']).size().reset_index(name='sample_num')

#Targeted SNV_counts, how many patients in cohort have a SNV in gene
tcounts = targeted.groupby(['Patient ID', 'GENE']).size().reset_index(name='p_count')

tfreq = tcounts.groupby(['GENE']).size().reset_index(name='SNV_count')
tfreq = tfreq.sort_values('SNV_count', ascending = False)

tfreq = tfreq.reset_index(drop=True)
tfreq = tfreq[tfreq['SNV_count'] >= 3]

temp = pd.DataFrame(tfreq.GENE)

tcounts = temp.merge(tcounts, how='left', on='GENE')
tcounts = tcounts.merge(sample_count, how='left', on='Patient ID')
tcounts.p_count = tcounts.p_count/tcounts.sample_num

if sort_by_median:
    med = tcounts.groupby('GENE').median()
    med = med.sort_values('p_count', ascending = False)
    med = med.reset_index()
    
    tcounts = med.merge(tcounts, how='left', on='GENE')
    tfreq = med.merge(tfreq, how='left', on='GENE')



if include_wes:
    wes = pd.read_csv('Y:/eritch/projects/ghent_m1rp/wxs/all_files_nov1st2020/all_maf.tsv', sep='\t')
    wes_snv = pd.read_csv('Y:/eritch/projects/ghent_m1rp/wxs/all_files_nov1st2020/all_snv.tsv', sep='\t')
    wes = wes[~wes['Tumor_Sample_Barcode'].str.contains('ID8')]
    wes = wes.loc[(wes['Tumor_Sample_Barcode'].str.contains('_PB'))|(wes['Tumor_Sample_Barcode'].str.contains('_RP'))]
    
    #WES SNV counts, how many patients in cohort have a SNV in gene
    wcounts = wes.copy()
    wcounts['Patient ID'] = wes['Tumor_Sample_Barcode'].str.split('_')
    wcounts['Patient ID'] = wcounts['Patient ID'].str[0] + '_' + wcounts['Patient ID'].str[1]
    wcounts = wcounts.groupby(['Patient ID', 'Hugo_Symbol']).size().reset_index(name='p_count')
    
    wfreq = wcounts.groupby(['Hugo_Symbol']).size().reset_index(name='SNV_count')
    wfreq = wfreq.sort_values('SNV_count', ascending = False)
    wfreq = wfreq.reset_index(drop=True)
    wfreq.columns = ['GENE','SNV_count']
    wfreq = temp.merge(wfreq, how='left', on='GENE')
    wfreq = med.merge(wfreq, how='left', on='GENE')

##############################################################################################
 
# creating the bar plot 
fig = plt.figure(figsize=(15,10))
gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [1,1])
plt.subplots_adjust(hspace=0.5)

ax = fig.add_subplot(gs[0,0])

ax.bar(tfreq['GENE'], tfreq['SNV_count'], color ='#424242', width = 0.8) 

if include_wes:
       ax.bar(wfreq['GENE'], wfreq['SNV_count'], color ='#b3b3b3', width = 0.4)  
       
ax.set_ylabel("# Patients with SNV") 
           
ax.set_ylim(0,15)
plt.xticks(rotation=90) 

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


# creating the bar plot 
ax2 = fig.add_subplot(gs[1,0],sharex=ax)

sns.stripplot(x="GENE",y="p_count",data=tcounts, ax=ax2, color ='#424242', zorder=0, dodge=True, size = 3)
sns.boxplot(x="GENE",y="p_count",data=tcounts, whis=np.inf,ax=ax2,color='grey',boxprops=dict(alpha=.1),zorder=10,dodge=True, width=0.7, linewidth = 1)
       
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.invert_yaxis() 
ax2.set_ylabel("% Samples with SNV per patient") 

 
fig.savefig('C:\\Users\\Sarah\\Desktop\\SNV_frequency_sorted_frequency.png', bbox_inches='tight', dpi = 600)

















