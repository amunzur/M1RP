# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:54:57 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import string

def get_expected_lr(cn, tc):
    return math.log2(cn*tc+2-2*tc)-1

# =============================================================================
# Obtain PTEN coords
# =============================================================================
chrom = 'chr10'
pten_refSeqID = 'NM_000314'

refgene = pd.read_csv('G:/Andy Murtha/panelDesign/refGene.txt', delimiter = '\t', low_memory = False, names = ['bin','name','chr', 'strand', 'txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], usecols = ['bin','name','chr','strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], dtype = {'bin': np.int32,'name': str,'chr': str,'strand': str, 'txStart': np.int32, 'txEnd': np.int32, 'cdsStart': np.int32, 'cdsEnd': np.int32, 'exonCount': np.int32, 'exonStarts': str,'exonEnds': str,'score': np.int32, 'name2': str, 'cdsStartStat': str, 'cdsEndStat': str,'exonFrames': str})

refgene = refgene[refgene['name'] == pten_refSeqID].reset_index()

exons = pd.DataFrame({'starts':refgene['exonStarts'][0].split(',')[:-1], 'ends':refgene['exonEnds'][0].split(',')[:-1]})

exons = exons.astype(float)

start = refgene['txStart'][0] - 1000
end = refgene['txEnd'][0] + 1000

# =============================================================================
# Constants
# =============================================================================

sample_cat = {'RP':'Primary',
              'PB':'Primary',
              'MLN': 'Met',
              'cfDNA': 'cfDNA'}

probe_path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Probe-level logratios (targeted panel)'
out_dir = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Probe level CNA analysis/PTEN loss analysis'

# =============================================================================
# import cna and tc
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] >= 0.1]

tc['Sample Category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits).apply(lambda x: sample_cat.get(x))
tc_copy = tc.copy().set_index('Sample ID')

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Import SNPS and limit to PTEN
# =============================================================================

snp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/M1RP_snps_melted.tsv', sep = '\t')
snp = snp[(snp['CHROM'] == 'chr10')&(snp['POSITION'] >= start) & (snp['POSITION'] <= end)]
snp['VAF'] = snp['VAF'].astype(float)
snp['VAF'] = (snp['VAF']-0.5).abs()+0.5
snp = snp[snp['Read depth'] >= 20]

# =============================================================================
# Loop over patients and plot
# =============================================================================

def plot_probe_lr(ax, path, s, start, end):
    probes = pd.read_csv(os.path.join(path, '%s_logratio.igv' % s), sep = '\t')
    probes = probes.rename(columns = {s:'log_ratio'})
    
    probes = probes[probes['CHROM'] == chrom]
    probes = probes[probes['START'] >= start]
    probes = probes[probes['END'] <= end]
    
    ax.scatter(probes['START'], probes['log_ratio'], color = 'k', s = 6, zorder = 1000)
    
    ax.set_xlim(start,end)
    ax.set_ylim(math.floor(probes['log_ratio'].min()),0.15)
    
    # if y == 0: 
    #     ax.set_title(pt+'\n'+sample)
    # else:
    #     ax.set_title(sample)
    # ax.xaxis.set_visible(False)
    
    ax.barh(0, exons['ends']-exons['starts'], left = exons['starts'], height = 0.15, zorder = 100,align = 'edge')
    
    lr_0 = get_expected_lr(0, tc_copy.at[s,'Final tNGS_TC'])
    ax.plot([start,end],[lr_0,lr_0], marker = None, lw = 0.8, color = '#3F60AC', zorder = 5100)
    
    lr_1 = get_expected_lr(1, tc_copy.at[s,'Final tNGS_TC'])
    ax.plot([start,end],[lr_1,lr_1], marker = None, lw = 0.8, color = '#9CC5E9', zorder = 500)
    
    ax.tick_params(bottom = False, labelbottom = False, labelsize = 6)
    ax.set_xlabel(s+', TC:%.2f' % (tc_copy.at[s,'Final tNGS_TC']), fontsize = 6)
    
def plot_SNP_vaf(ax, snps, s):
    snps['VAF'] = snps['VAF'].astype(float)
    ax.scatter(snps['POSITION'],snps['VAF'], lw = 0, c = 'k', s = 10, clip_on = False)
    
    ax.set_ylim(0.5,1)
    ax.tick_params(bottom = False, labelbottom = False, labelsize = 6)


# for pt in ['ID%s' % (i+1) for i in np.arange(43)]:
for pt in ['ID38']:
    pt_samples = tc[(tc['Patient ID'] == pt)&(tc['Final tNGS_TC'] >= 0.2)]
    if pt == 'ID28':
        pt = pt
    pt_snp = snp[snp['Patient ID'] == pt].copy()
    if len(pt_samples) == 0: continue;
    elif len(pt_samples) == 1:
        if len(pt_snp) == 0:
            fig,ax = plt.subplots(ncols = len(pt_samples))
            plot_probe_lr(ax, probe_path, pt_samples['Sample ID'].tolist()[0], start, end)
        else:
            fig,axs = plt.subplots(ncols = 1, nrows = 2, sharex = True)
            s = pt_samples['Sample ID'].tolist()[0]
            plot_probe_lr(axs[0], probe_path, s, start, end)
            plot_SNP_vaf(axs[1], pt_snp[pt_snp['Sample ID'] == s].copy(), s)
    else:
        if len(pt_snp) == 0:
            fig,axs = plt.subplots(ncols = len(pt_samples), figsize = (1.75*len(pt_samples),2), sharex = True)
            for ax, s in zip(axs, pt_samples['Sample ID'].tolist()):
                plot_probe_lr(ax, probe_path, s, start, end)
        else:
            fig,axs = plt.subplots(ncols = len(pt_samples), figsize = (1.75*len(pt_samples),2), nrows = 2, sharex = True)
            for i, s in enumerate(pt_samples['Sample ID'].tolist()):
                plot_probe_lr(axs[0][i], probe_path, s, start, end)
                plot_SNP_vaf(axs[1][i], pt_snp[pt_snp['Sample ID'] == s].copy(), s)
        
        
    fig.tight_layout()
    fig.suptitle('%s' % pt)
    fig.savefig(os.path.join(out_dir, '%s.PTEN.pdf' % pt))

    
    