# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 11:42:42 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pt = 'ID1'

chr_size = pd.read_excel('G:/Andy Murtha/panelDesign/chr_size.xlsx')
chr_size = chr_size[~chr_size['CHR'].isin(['chrY','chrX'])]
chr_size = chr_size.set_index('CHR')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
    
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc.set_index('Sample ID')

fig, ax_matrix = plt.subplots(ncols = 22, nrows = 43, figsize = (8.5,11), sharex = 'col', gridspec_kw={'width_ratios':chr_size['Total length (bp)']})

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC','M1RP_ID30_RP3','M1RP_ID15_PB5','M1RP_ID18_PB8']

truncal_losses = []
truncal_gains = []

branch_losses = []
branch_gains = []

for i in np.arange(43):
    pt = 'ID%i' % (i+1)
    segs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 1.25)/M1RP_%s_-segments.txt' % pt, sep = '\t')

    samples = segs.columns.tolist()[3:]
    segs = segs.drop([s for s in samples if tc.at[s, 'Final tNGS_TC'] < 0.2 or s in secondary_tumors], axis = 1)
    samples = [s for s in samples if tc.at[s, 'Final tNGS_TC'] >= 0.2 and s not in secondary_tumors]
    
    segs = segs[~segs['CHR'].isin(['chrM','chrY','chrX'])]
    segs['P'] = 2
    segs.loc[segs['CHR'].isin(['chrX']), 'P'] = 1
    
    # =============================================================================
    # 
    # =============================================================================
    
    for s in samples:
        segs[s+'_lr'] = segs[s].str.split(':').str[0].astype(float)
        segs[s+'_snp'] = segs[s].str.split(':').str[1].astype(float)
        segs[s+'_tc'] = tc.at[s, 'Final tNGS_TC']
        segs[s+'_cn'] = (segs['P']*(segs[s+'_tc']+2**segs[s+'_lr']-1))/segs[s+'_tc']
        segs[s+'_cnRounded'] = segs[s+'_cn'].round().clip(lower = 1, upper = 3)
    
    segs = segs.drop(samples, axis = 1)
    
    # =============================================================================
    # 
    # =============================================================================
    
    loss_counts = segs.copy()
    loss_counts = loss_counts.set_index(['CHR','START','END'])
    loss_counts = loss_counts[[s+'_cnRounded' for s in samples]]
    gain_counts = segs.copy()
    gain_counts = gain_counts.set_index(['CHR','START','END'])
    gain_counts = gain_counts[[s+'_cnRounded' for s in samples]]
    
    for s in samples:
        s = s+'_cnRounded'
        loss_counts.loc[loss_counts[s] != 1, s] = 0
        gain_counts.loc[gain_counts[s] != 3, s] = 0
        gain_counts.loc[gain_counts[s] == 3, s] = 1
        
    loss_counts['sum'] = loss_counts.sum(axis = 1)
    gain_counts['sum'] = gain_counts.sum(axis = 1)
    
    axs = ax_matrix[i]
    axs_dict = dict(zip(gain_counts.index.get_level_values(0).unique().tolist(), axs))
    
    for index, row in loss_counts.iterrows():
        ax = axs_dict.get(index[0])
        c = index[0]
        alpha = row['sum'] / len(samples)
        ax.barh(0,index[2]-index[1],left = index[1], color = 'blue', alpha = alpha)
        
    for index, row in gain_counts.iterrows():
        ax = axs_dict.get(index[0])
        c = index[0]
        alpha = row['sum'] / len(samples)
        ax.barh(1,index[2]-index[1],left = index[1], color = 'red', alpha = alpha)
        
    for c in gain_counts.index.get_level_values(0).unique().tolist():
        ax = axs_dict.get(c)
        ax.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False)
        ax.spines['left'].set_visible(False)
        ax.set_xlim(0,chr_size.at[c,'Total length (bp)'])
        if pt == 'ID43':
            ax.set_xlabel(c, rotation = 90, fontsize = 6)
        if c == 'chr1':
            ax.set_ylabel(pt, rotation = 0, fontsize = 6, ha = 'right')
    
    loss_tmp = loss_counts.copy()
    loss_tmp = loss_tmp.reset_index()
    loss_tmp = loss_tmp[loss_tmp['sum'] == len(samples)]
    loss_tmp = loss_tmp[[col for col in loss_tmp.columns.tolist() if col in ['CHR','START','END']]]
    loss_tmp['Patient ID'] = pt
    truncal_losses.append(loss_tmp)
    
    gain_tmp = gain_counts.copy()
    gain_tmp = gain_tmp.reset_index()
    gain_tmp = gain_tmp[gain_tmp['sum'] == len(samples)]
    gain_tmp = gain_tmp[[col for col in gain_tmp.columns.tolist() if col in ['CHR','START','END']]]
    gain_tmp['Patient ID'] = pt
    truncal_gains.append(gain_tmp)

    loss_tmp = loss_counts.copy()
    loss_tmp = loss_tmp.reset_index()
    loss_tmp = loss_tmp[(loss_tmp['sum'] < len(samples))&(loss_tmp['sum'] > 0)]
    loss_tmp = loss_tmp[[col for col in loss_tmp.columns.tolist() if col in ['CHR','START','END']]]
    loss_tmp['Patient ID'] = pt
    branch_losses.append(loss_tmp)
    
    gain_tmp = gain_counts.copy()
    gain_tmp = gain_tmp.reset_index()
    gain_tmp = gain_tmp[(gain_tmp['sum'] < len(samples))&(gain_tmp['sum'] > 0)]
    gain_tmp = gain_tmp[[col for col in gain_tmp.columns.tolist() if col in ['CHR','START','END']]]
    gain_tmp['Patient ID'] = pt
    branch_gains.append(gain_tmp)
             # 
            
plt.tight_layout()
# plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/wes_cna_truncality.pdf')

truncal_losses_df = pd.concat(truncal_losses, ignore_index = True)
truncal_gains_df = pd.concat(truncal_gains, ignore_index = True)

tl_count_dict = {}
tg_count_dict = {}

fig, axs = plt.subplots(ncols = 22, gridspec_kw = {'width_ratios':chr_size['Total length (bp)']}, sharey = True)

for ax, (c, size) in zip(axs, chr_size['Total length (bp)'].iteritems()):
    c_tl = truncal_losses_df[truncal_losses_df['CHR'] == c].copy()
    c_tg = truncal_gains_df[truncal_gains_df['CHR'] == c].copy()
    tl_count_dict[c] = {}
    tg_count_dict[c] = {}
    for start in np.arange(0,size,10**6):
        end = start + 10**6 -1
        
        # =============================================================================
        # Count truncal losses for each bin        
        # =============================================================================
        bin_loss = len(c_tl[((c_tl['START'] >= start)&(c_tl['START'] <= end))|((c_tl['END'] >= start)&(c_tl['END'] <= end))|((c_tl['START'] <= start)&(c_tl['END'] >= end))])
        tmp = tl_count_dict.get(c)
        tmp[(start+end)//2] = bin_loss*-1
        tl_count_dict[c]= tmp
        if len(c_tl) > 0:
            c_tl = c_tl
        
        # =============================================================================
        # Count truncal gains for each bin        
        # =============================================================================
        bin_gain = len(c_tg[((c_tg['START'] >= start)&(c_tg['START'] <= end))|\
                    ((c_tg['END'] >= start)&(c_tg['END'] <= end))|\
                        ((c_tg['START'] <= start)&(c_tg['END'] >= end))])
        tmp = tg_count_dict.get(c)
        tmp[(start+end)//2] = bin_gain
        tg_count_dict[c]= tmp
        
    
    ax.plot(tl_count_dict[c].keys(), tl_count_dict[c].values(), color = 'blue', marker = None, zorder = 100)
    ax.plot(tg_count_dict[c].keys(), tg_count_dict[c].values(), color = 'red', marker = None, zorder = 100)
    ax.plot([0,size], [0,0], ls = 'dashed', color = 'grey', zorder = 0)
    ax.set_xlim(0, size)
    ax.tick_params(bottom = False, labelbottom = False)
    ax.set_xlabel(c, rotation = 90, fontsize = 6)
    if c != 'chr1':
        ax.spines['left'].set_visible(False)
        ax.tick_params(left = False)
        
del c_tl, c_tg, tl_count_dict, tg_count_dict
# =============================================================================
#         
# =============================================================================
        
branch_losses_df = pd.concat(branch_losses, ignore_index = True)
branch_gains_df = pd.concat(branch_gains, ignore_index = True)

bl_count_dict = {}
bg_count_dict = {}

        
fig, axs = plt.subplots(ncols = 22, gridspec_kw = {'width_ratios':chr_size['Total length (bp)']}, sharey = True)

for ax, (c, size) in zip(axs, chr_size['Total length (bp)'].iteritems()):
    c_bl = branch_losses_df[branch_losses_df['CHR'] == c].copy()
    c_bg = branch_gains_df[branch_gains_df['CHR'] == c].copy()
    bl_count_dict[c] = {}
    bg_count_dict[c] = {}
    for start in np.arange(0,size,10**6):
        end = start + 10**6 -1
        
        # =============================================================================
        # Count truncal losses for each bin        
        # =============================================================================
        bin_loss = len(c_bl[((c_bl['START'] >= start)&(c_bl['START'] <= end))|((c_bl['END'] >= start)&(c_bl['END'] <= end))|((c_bl['START'] <= start)&(c_bl['END'] >= end))])
        tmp = bl_count_dict.get(c)
        tmp[(start+end)//2] = bin_loss*-1
        bl_count_dict[c]= tmp
        if len(c_bl) > 0:
            c_bl = c_bl
        
        # =============================================================================
        # Count truncal gains for each bin        
        # =============================================================================
        bin_gain = len(c_bg[((c_bg['START'] >= start)&(c_bg['START'] <= end))|\
                    ((c_bg['END'] >= start)&(c_bg['END'] <= end))|\
                        ((c_bg['START'] <= start)&(c_bg['END'] >= end))])
        tmp = bg_count_dict.get(c)
        tmp[(start+end)//2] = bin_gain
        bg_count_dict[c]= tmp
        
    
    ax.plot(bl_count_dict[c].keys(), bl_count_dict[c].values(), color = 'blue', marker = None, zorder = 100)
    ax.plot(bg_count_dict[c].keys(), bg_count_dict[c].values(), color = 'red', marker = None, zorder = 100)
    ax.plot([0,size], [0,0], ls = 'dashed', color = 'grey', zorder = 0)
    ax.set_xlim(0, size)
    ax.tick_params(bottom = False, labelbottom = False)
    ax.set_xlabel(c, rotation = 90, fontsize = 6)
    if c != 'chr1':
        ax.spines['left'].set_visible(False)
        ax.tick_params(left = False)
        