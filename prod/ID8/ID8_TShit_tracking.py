# -*- coding: utf-8 -*-
"""
Created on Tue May  4 11:43:09 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string
import scipy.stats as stats

# =============================================================================
# Constants
# =============================================================================

def get_snp_estimation(tc, cat):
    tc = tc[tc['Sample Category'] == cat]
    est_vaf = (1/(2-tc['Final tNGS_TC']))-0.5
    return est_vaf;

sample_cat = {'RP':'Primary',
              'PB':'Primary',
              'MLN': 'Met',
              'MB':'Met',
              'cfDNA': 'cfDNA'}

color_dict = {'Synonymous':'grey',
              "3'-UTR":'grey',
              'Missense':'green',
              'Intergenic':'grey',
              'Intronic':'grey',
              'Frameshift':'yellow',
              'Stopgain':'yellow',
              "5'-UTR":'grey'}
fprops = {'marker':'o','markeredgewidth':0,'markersize':2,'markerfacecolor':'k'}
# order = {
#  'M1RP_ID10_PB5':0,
#  'M1RP_ID10_PB8':1,
#  'M1RP_ID10_PB7':1,
#  'M1RP_ID10_RP2':0,
#  'M1RP_ID10_PB2':2,
#  'M1RP_ID10_RP5':0,
#  'M1RP_ID10_PB1':0,
#  'M1RP_ID10_PB3':1,
#  'M1RP_ID10_PB4':1,
#  'M1RP_ID10_RP1':1,
#  'M1RP_ID10_RP3':3,
#  'M1RP_ID10_RP4':1,
#  'M1RP_ID10_MB1':1,
#  'M1RP_ID10_cfDNA_2018Dec03':1}

# =============================================================================
# Import samples, cn data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
snp = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/M1RP_snps_melted.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] >= 0.20]

tc['Sample Category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits).apply(lambda x: sample_cat.get(x))
# tc = tc[tc['Sample Category'] != 'cfDNA']

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Get ID2
# =============================================================================

pt = 'ID8'

tc = tc[tc['Patient ID'] == pt]
cn = cn[cn['Patient ID'] == pt]
snp = snp[snp['Sample ID'].isin(tc['Sample ID'].tolist())]

# tc['order'] = tc['Sample ID'].map(order)
# tc = tc.sort_values(['order','Final tNGS_TC'], ascending = True)
tc = tc.sort_values(['Final tNGS_TC'], ascending = True)

# =============================================================================
# get snp deviation
# =============================================================================

snp['deviation'] = (snp['VAF']-0.5).abs()

# =============================================================================
# Get sample counts
# =============================================================================

n_pri = len(tc[tc['Sample Category'] == 'Primary'])
n_met = len(tc[tc['Sample Category'] == 'Met'])
n_cf = len(tc[tc['Sample Category'] == 'cfDNA'])

# =============================================================================
# Set up mutations
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')

muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
muts['GENE'] = muts['GENE'].replace({'BIVM-ERCC5;ERCC5':'ERCC5'}).fillna('Intergenic')

genes = muts['GENE'].unique().tolist()
muts['y'] = muts['GENE'].map(dict(zip(genes,np.arange(len(genes)))))

muts['color'] = muts['EFFECT'].str.split(' ').str[0].map(color_dict)

p_muts = muts[muts['Sample ID'].isin(tc[tc['Sample Category'] == 'Primary']['Sample ID'].tolist())].copy()
m_muts = muts[muts['Sample ID'].isin(tc[tc['Sample Category'] == 'Met']['Sample ID'].tolist())].copy()
cf_muts = muts[muts['Sample ID'].isin(tc[tc['Sample Category'] == 'cfDNA']['Sample ID'].tolist())].copy()

# =============================================================================
# Loop over genes and samples. Plot estimated copy number and SNP deviation for
# each sample
# =============================================================================

goi = ['TP53','PTEN','BRCA2','MSH2','MSH6']

# =============================================================================
# Get SNP counts per gene
# =============================================================================

snp_count = snp.copy()
snp_count = snp_count.drop_duplicates(['CHROM','POSITION','REF','ALT'])
snp_count = snp_count.groupby(['GENE']).count().reset_index()

snp_count = pd.DataFrame({'GENE':goi}).merge(snp_count, on = 'GENE', how = 'left').fillna(0)

snp_count = snp_count[['GENE','CHROM']]
snp_count.columns = ['GENE','Count']
no_snp = len(snp_count[snp_count['Count'] <= 1])
snp_count = snp_count.set_index('GENE')

# =============================================================================
# Plots
# =============================================================================


fig,axs = plt.subplots(ncols = 3, nrows = 2*len(goi)+1-no_snp, gridspec_kw={'width_ratios': [n_pri,n_met,n_cf]}, sharey = 'row', sharex = 'col', figsize = (6,7))

a = 0
first = True
for gene in goi:
    primary = []
    p_order = []
    p_snps = []
    
    metastatic = []
    m_order = []
    m_snps = []
    
    cfdna = []
    c_order = []
    c_snps = []
    
    for index,row in tc.iterrows():
        s = []
        sg_snps = snp[(snp['Sample ID'] == row['Sample ID'])&(snp['GENE'] == gene)].copy()
        
        for i in range(1000):
            s_tc = sorted([0.0001,row['Final tNGS_TC']*2**np.random.normal(0,0.1)+ np.random.normal(0,0.02),1])[1]
            lr = cn.at[(row['Sample ID'],gene), 'Log_ratio']
            s.append(2*(s_tc+2**lr-1)/s_tc)
        if row['Sample Category'] == 'Primary':
            primary.append(s)
            p_order.append(row['Sample ID'])
            p_snps.append(sg_snps['deviation'])
        elif row['Sample Category'] == 'Met':
            metastatic.append(s)
            m_order.append(row['Sample ID'])
            m_snps.append(sg_snps['deviation'])
        else:
            cfdna.append(s)
            c_order.append(row['Sample ID'])
            c_snps.append(sg_snps['deviation'])
         
    n_met = len(tc[tc['Sample Category'] == 'Met'])
    n_pri = len(tc[tc['Sample Category'] == 'Primary'])
    n_cf = len(tc[tc['Sample Category'] == 'cfDNA'])
    
    # =============================================================================
    # Boxplots for estimated copy number
    # =============================================================================
    
    ax1 = axs[a][0]
    ax2 = axs[a][1]
    ax3 = axs[a][2]
    ax4 = axs[a+1][0]
    ax5 = axs[a+1][1]
    ax6 = axs[a+1][2]
    
       
    ax1.boxplot(primary, flierprops = fprops)
    ax2.boxplot(metastatic, flierprops = fprops)
    # ax3.boxplot(cfdna, flierprops = fprops)
    
    ax1.plot([-1,n_pri],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
    ax2.plot([-1,n_met],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
    # ax3.plot([-1,n_cf],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
    
    ax1.plot([-1,n_pri],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
    ax2.plot([-1,n_met],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
    # ax3.plot([-1,n_cf],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
    
    ax1.set_ylim(0,3)
    ax1.set_yticks(np.arange(0,4,1))
    ax1.tick_params(labelsize = 6)
    ax1.set_ylabel('%s\ncopy number' % gene, fontsize = 6)
    
    # =============================================================================
    # Boxplots and scatter for SNPs  
    # =============================================================================
    
    # ax3.boxplot(p_snps, showfliers = False, zorder = 10)
    a+=1
    if snp_count.at[gene,'Count'] > 1:
        a+=1
        for i, sg_snp in enumerate(p_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax4.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100)
            ax4.scatter(np.arange(1,n_pri+1,1),get_snp_estimation(tc,'Primary'),c='blue',marker='_',lw=1)
        for i, sg_snp in enumerate(m_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax5.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100)
            ax5.scatter(np.arange(1,n_met+1,1),get_snp_estimation(tc,'Met'), c = 'blue',marker='_',lw=1)
        for i, sg_snp in enumerate(c_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax6.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100)
            ax6.scatter(np.arange(1,n_cf+1,1),get_snp_estimation(tc,'cfDNA'), c = 'blue',marker='_',lw=1)
        ax4.set_ylim(-0.02, 0.5)
        ax4.set_yticks(np.arange(0,0.6,0.1))
        ax4.tick_params(labelsize = 6)
        ax4.set_ylabel('%s\nSNP deviation' % gene, fontsize = 6)
        if first == True:
            ax4.legend(labels = ['SNP deviation','Est. hetz. loss SNP deviation'], fontsize = 6, loc = 'upper left', bbox_to_anchor=(-.02,1.1), frameon=False)
            first = False
        
    
# =============================================================================
# Mutation oncoprint    
# =============================================================================
   
ax1 = axs[-1][0]
ax2 = axs[-1][1]
# ax3 = axs[-1][2] 
   
p_muts['x'] = p_muts['Sample ID'].map(dict(zip(p_order,np.arange(1,len(p_order)+1,1))))
m_muts['x'] = m_muts['Sample ID'].map(dict(zip(m_order,np.arange(1,len(m_order)+1,1))))
# cf_muts['x'] = cf_muts['Sample ID'].map(dict(zip(c_order,np.arange(1,len(c_order)+1,1))))

ax1.scatter(p_muts['x'],p_muts['y'],c=p_muts['color'],marker = 's', s = 6)
ax2.scatter(m_muts['x'],m_muts['y'],c=m_muts['color'],marker = 's', s = 6)
# ax3.scatter(cf_muts['x'],cf_muts['y'],c=cf_muts['color'],marker = 's', s = 6)

ax1.set_yticks(np.arange(0,len(genes),1))
ax1.set_yticklabels(genes, fontsize = 6)
ax1.set_ylim(-0.4,len(genes)-0.6)
    
# =============================================================================
# Axis adjustments    
# =============================================================================
    
ax1.set_xticks(np.arange(1,n_pri+1,1))
ax2.set_xticks(np.arange(1,n_met+1,1))
# ax3.set_xticks(np.arange(1,n_cf+1,1))
ax1.set_xticklabels(pd.Series(p_order).str.split('_').str[2], rotation = 90)
ax2.set_xticklabels(pd.Series(m_order).str.split('_').str[2], rotation = 90)
# ax3.set_xticklabels(pd.Series(c_order).str.split('_').str[2], rotation = 90)

ax1.set_xlim(0.4, n_pri+0.6)
ax2.set_xlim(0.4, n_met+0.6)
# ax3.set_xlim(0.4, n_cf+0.6)   
    

axs[-1][0].tick_params(labelsize = 6)
axs[-1][1].tick_params(labelsize = 6)
# axs[-1][2].tick_params(labelsize = 6)
    
fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Patient specific/%s/%s_TStracking.pdf' % (pt, pt))
