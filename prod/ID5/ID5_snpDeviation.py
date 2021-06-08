# -*- coding: utf-8 -*-
"""
Created on Fri May  7 14:16:30 2021

@author: amurtha
"""

from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

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
              'cfDNA': 'cfDNA'}

color_dict = {'Synonymous':'grey', "3'-UTR":'grey', 'Missense':'green', 'Intergenic':'grey'}
fprops = {'marker':'o','markeredgewidth':0,'markersize':2,'markerfacecolor':'k'}


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

tc = tc[tc['Patient ID'] == 'ID5']
cn = cn[cn['Patient ID'] == 'ID5']
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

goi = ['TP53','RB1','PTEN']

# =============================================================================
# Get SNP counts per gene
# =============================================================================

snp_count = snp.copy()
snp_count = snp_count.drop_duplicates(['CHROM','POSITION','REF','ALT'])
snp_count = snp_count.groupby(['GENE']).count().reset_index()

snp_count = pd.DataFrame({'GENE':goi}).merge(snp_count, on = 'GENE', how = 'left').fillna(0)

snp_count = snp_count[['GENE','CHROM']]
snp_count.columns = ['GENE','Count']
no_snp = len(snp_count[snp_count['Count'] == 0])
snp_count = snp_count.set_index('GENE')

# =============================================================================
# Plots
# =============================================================================


fig,axs = plt.subplots(ncols = 3, nrows = 3, gridspec_kw={'width_ratios': [n_pri,n_met,n_cf]}, sharey = 'row', sharex = 'col', figsize = (3,2.5))

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
    
    ax1.set_ylim(0,3)
    ax1.set_yticks(np.arange(0,4,1))
    ax1.tick_params(labelsize = 6)
    ax1.set_ylabel('%s\ncopy number' % gene, fontsize = 6)
    
    # =============================================================================
    # Boxplots and scatter for SNPs  
    # =============================================================================
    
    # ax3.boxplot(p_snps, showfliers = False, zorder = 10)
    if snp_count.at[gene,'Count'] > 0 and gene != 'PTEN':
        a+=1
        for i, sg_snp in enumerate(p_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax1.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100, marker = 'o', label = 'SNP deviation', clip_on = False)
            ax1.scatter(np.arange(1,n_pri+1,1),get_snp_estimation(tc,'Primary'),c='blue',marker='_',lw=0.8)
        for i, sg_snp in enumerate(m_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax2.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100, clip_on = False)
            ax2.scatter(np.arange(1,n_met+1,1),get_snp_estimation(tc,'Met'), c = 'blue',marker='_',lw=0.8)
        for i, sg_snp in enumerate(c_snps):
            xs = pd.Series([i+1]*len(sg_snps), dtype = np.float64).apply(lambda x: x + np.random.uniform(low=-0.2, high = 0.2))
            ax3.scatter(xs,sg_snp, lw = 0, color = 'k', alpha = 0.6, s = 8, zorder = 100, clip_on = False)
            ax3.scatter(np.arange(1,n_cf+1,1),get_snp_estimation(tc,'cfDNA'),c = 'blue',marker='_',lw=0.8)
        ax1.set_ylim(-0.02, 0.5)
        ax1.set_yticks(np.arange(0,0.6,0.1))
        ax1.tick_params(labelsize = 6, pad = 1)
        ax1.set_ylabel('%s\nSNP deviation' % gene, fontsize = 6)
    elif gene == 'PTEN':
        ax1.boxplot(primary, flierprops = fprops)
        ax2.boxplot(metastatic, flierprops = fprops)
        ax3.boxplot(cfdna, flierprops = fprops)
        
        ax1.plot([-1,n_pri],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
        ax2.plot([-1,n_met],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
        ax3.plot([-1,n_cf],[2,2],lw = 0.5, ls = 'dashed', color = 'grey')
        
        ax1.plot([-1,n_pri],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
        ax2.plot([-1,n_met],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
        ax3.plot([-1,n_cf],[1,1],lw = 0.5, ls = 'dashed', color = 'blue')
        

handles = [Line2D([0,0],[0,0],lw = 0, marker = 'o', markerfacecolor='k',alpha=0.6, markersize = 3, markeredgewidth = 0),
           Line2D([0,0],[0,0],lw = 0.8, color = 'blue', marker = None)]
labels = ['SNP deviaiton','Est. hetz. loss deviation']

axs[0][0].legend(handles, labels, fontsize = 6, loc = 'upper left', bbox_to_anchor=(-.03,1.17), frameon=False, handlelength = 0.8, labelspacing = 0.35)

ax1.set_xticks(np.arange(1,n_pri+1,1))
ax2.set_xticks(np.arange(1,n_met+1,1))
ax3.set_xticks(np.arange(1,n_cf+1,1))
ax1.set_xticklabels(pd.Series(p_order).str.split('_').str[2], rotation = 90)
ax2.set_xticklabels(pd.Series(m_order).str.split('_').str[2], rotation = 90)
ax3.set_xticklabels(pd.Series(c_order).str.split('_').str[2], rotation = 90)

ax1.set_xlim(0.4, n_pri+0.6)
ax2.set_xlim(0.4, n_met+0.6)
ax3.set_xlim(0.4, n_cf+0.6)   
    

axs[-1][0].tick_params(labelsize = 6, pad = 1)
axs[-1][1].tick_params(labelsize = 6, pad = 1)
axs[-1][2].tick_params(labelsize = 6, pad = 1)


        
fig.tight_layout()
fig.subplots_adjust(left = 0.15,right = 0.96, top = 0.94, bottom = 0.15, hspace = 0.3)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Patient specific/ID5/TS_snpDeviation.pdf')