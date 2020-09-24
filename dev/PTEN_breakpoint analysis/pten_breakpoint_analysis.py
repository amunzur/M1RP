# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 12:54:40 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec

# =============================================================================
# Get PTEN location
# =============================================================================

chrom = 'chr10'
pten_refSeqID = 'NM_000314'

# =============================================================================
# Import log-ratio of whole gene from targeted and WES
# =============================================================================

t_vs_wes = pd.read_csv('https://docs.google.com/spreadsheets/d/1zi7UiPteDA3VU4jtcp1dxET4UGCMKvtLqab-7G13_oU/export?format=csv&gid=1176347443')
t_vs_wes['Patient ID'] = t_vs_wes['Patient ID'].str.split('_').str[1]

t_vs_wes = t_vs_wes[t_vs_wes['Gene'] == 'PTEN']
t_vs_wes['lrd_diff'] = (t_vs_wes['Log_ratio_Targeted'] - t_vs_wes['Log_ratio_WES']).abs()
t_vs_wes['cn_diff'] = (t_vs_wes['CNV_Targeted'] - t_vs_wes['CN_WES']).abs()

pts = t_vs_wes[(t_vs_wes['lrd_diff'] > 0.5)]['Patient ID'].unique().tolist()

t_vs_wes = t_vs_wes[t_vs_wes['Patient ID'].isin(pts)]
t_vs_wes = t_vs_wes[t_vs_wes['Patient ID'] != 'ID8']

max_samples = 0
sample_dict = {}
for pt in pts:
    nsamples = len(t_vs_wes[t_vs_wes['Patient ID'] == pt]['Sample ID'].tolist())
    if nsamples > max_samples:
        max_samples = nsamples
    sample_dict[pt] = t_vs_wes[t_vs_wes['Patient ID'] == pt]['Sample ID'].tolist()

# =============================================================================
# Import refGene and get exon start/end points
# =============================================================================

refgene = pd.read_csv('G:/Andy Murtha/panelDesign/refGene.txt', delimiter = '\t', low_memory = False, names = ['bin','name','chr', 'strand', 'txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], usecols = ['bin','name','chr','strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], dtype = {'bin': np.int32,'name': np.str,'chr': np.str,'strand': np.str, 'txStart': np.int32, 'txEnd': np.int32, 'cdsStart': np.int32, 'cdsEnd': np.int32, 'exonCount': np.int32, 'exonStarts': np.str,'exonEnds': np.str,'score': np.int32, 'name2': np.str, 'cdsStartStat': np.str, 'cdsEndStat': np.str,'exonFrames': np.str})

refgene = refgene[refgene['name'] == pten_refSeqID].reset_index()

exons = pd.DataFrame({'starts':refgene['exonStarts'][0].split(',')[:-1], 'ends':refgene['exonEnds'][0].split(',')[:-1]})

exons = exons.astype(float)

start = refgene['txStart'][0]
end = refgene['txEnd'][0]

# =============================================================================
# Create plot of PTEN log-ratios across gene
# =============================================================================

fig = plt.figure(figsize = (2.5*len(sample_dict)+1,1.5*max_samples))
gs = gridspec.GridSpec(max_samples,len(sample_dict))

for x,pt in enumerate(sample_dict.keys()):
    for y,sample in enumerate(sample_dict.get(pt)):
        ax = fig.add_subplot(gs[y,x])
        cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)/%s_logratio.igv' % sample, sep = '\t')
        cn = cn.rename(columns = {sample:'log_ratio'})
        
        cn = cn[cn['CHROM'] == chrom]
        cn = cn[cn['START'] >= start]
        cn = cn[cn['END'] <= end]
        
        ax.scatter(cn['START'], cn['log_ratio'], color = 'k', s = 6)
        
        ax.set_xlim(start,end)
        ax.set_ylim(math.floor(cn['log_ratio'].min()),0)
        
        if y == 0: 
            ax.set_title(pt+'\n'+sample)
        else:
            ax.set_title(sample)
        ax.xaxis.set_visible(False)
        
        ax.barh(0, exons['ends']-exons['starts'], left = exons['starts'], height = math.floor(cn['log_ratio'].min()) / 5)
        
        pt_comparison = t_vs_wes[t_vs_wes['Sample ID'] == sample]
        pt_comparison = pt_comparison[pt_comparison['Gene'] == 'PTEN'].reset_index()
        print(len(pt_comparison))
        
        targeted = min(0,pt_comparison['Log_ratio_Targeted'][0])
        wes = min(0,pt_comparison['Log_ratio_WES'][0])
        ax.plot([start,end],[targeted,targeted], color = 'red')
        ax.plot([start,end],[wes,wes], color = 'red', linestyle = 'dashed')
        
        
        
fig.tight_layout()
gs.update(wspace = 0.2, right = 0.98)
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/PTEN_breakpoint analysis/pten_breakpoint_analysis.pdf')