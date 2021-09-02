# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 12:50:27 2021

@author: amurtha
"""


import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from statsmodels.stats import proportion
from matplotlib import gridspec
import seaborn as sns
from scipy.stats import linregress
import matplotlib.patches as mpatches
import numpy as np

# Set some default fonts and style for matplotlib etc
mpl.rc('hatch', color = '#FFFFFF', linewidth = 0.2)
mpl.rcParams['font.size'] = 6
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['legend.handletextpad'] = '0.7'
mpl.rcParams['legend.labelspacing'] = '0.2'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['legend.handlelength'] = 0.7
plt.rcParams['legend.handleheight'] = 0.70
plt.rcParams["font.family"] = "Arial"



# =============================================================================
# Constants
# =============================================================================

#Genes to include
genes = ['SPOP', 'FOXA1', 'TP53', 'PTEN', 'RB1', 'AR', 'BRCA2', 'CDK12']

#Tumour suppressors to filter out unexpected changes (like gains in TS)
ts= ["TP53", "RB1", "PTEN", "BRCA2", "SPOP"]

#Genes to include for copy number changes
cn_genes= ["TP53", "RB1", "PTEN", "BRCA2", "AR"]

# Dictionaries to order the cohorts and genes on the barchart
cohort_dict = {"Primary TCGA PanCancer Atlas": 1, "SU2C mCRPC": 3, "MSK mCSPC":2, "M1RP de novo CSPC": 0}

gene_dict = {"SPOP\nmut": 3, "FOXA1\nmut": 1, "TP53\nmut": 0,  "PTEN\nmut": 2,
             "RB1\nmut": 4, "BRCA2\nmut": 5, "AR\nmut": 7, 'CDK12\nmut':6,
             "TP53\ndeep del": 10, "PTEN\ndeep del":8, "RB1\ndeep del": 9,"BRCA2\ndeep del": 12, "AR\namp":11 }

#Types of mutations to include
clinical_muts =["Frameshift", "Stopgain", "Splice", "indel", "Missense"]


# =============================================================================
# Import data
# =============================================================================
#Import files

# Paths for published cohorts
tcga_mut_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/prostate_adenocarcinoma_tcga_pancancer_atlas_mutations.txt"
tcga_cn_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/prostate_adenocarcinoma_tcga_pancancer_atlas_copy_number.txt"

msk_mut_path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/msk_mcspc_mutations.txt'
msk_cn_path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/msk_mcspc_copy_number.txt'

su2c_mut_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/SU2C_2019_mutations.txt"
su2c_cn_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/cBioportal/SU2C_2019_copy_number.txt"

# paths for M1RP cohort
snv_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv"

#To exclude the dependent mutations, use below:
#snv_path = "C:/Users/Sarah/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations.tsv"

cnv_path = "C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv"


#Load m1 data
m1_snvs = pd.read_csv(snv_path, sep='\t', header=[0])
m1_cnvs = pd.read_csv(cnv_path, sep='\t', header=[0])

# Load cbioportal data
tcga_snv = pd.read_csv(tcga_mut_path, sep='\t', header=[0])
tcga_cn = pd.read_csv(tcga_cn_path, sep='\t', header=[0])
su_snv = pd.read_csv(su2c_mut_path, sep='\t', header=[0])
su_cn = pd.read_csv(su2c_cn_path, sep='\t', header=[0])
msk_snv = pd.read_csv(msk_mut_path, sep = '\t', header=[0])
msk_cn = pd.read_csv(msk_cn_path, sep='\t', header=[0])

#Delete path variables
del tcga_mut_path, tcga_cn_path, su2c_mut_path, su2c_cn_path, snv_path, cnv_path, msk_cn_path, msk_mut_path

# Get total number of patients in M1RP
pt_num = len(m1_cnvs['Patient ID'].unique())
# =============================================================================
# Format data
# =============================================================================
# Label the cohorts
tcga_snv["Cohort"] = "Primary TCGA PanCancer Atlas"
tcga_cn["Cohort"] = "Primary TCGA PanCancer Atlas"

su_snv["Cohort"] = "SU2C mCRPC"
su_cn["Cohort"] = "SU2C mCRPC"

msk_snv["Cohort"] = 'MSK mCSPC'
msk_cn["Cohort"] = 'MSK mCSPC'

m1_snvs["Cohort"] = "M1RP de novo CSPC"
m1_cnvs["Cohort"] = "M1RP de novo CSPC"


# Combine cbioportal dataframes into one
snvs = pd.concat([tcga_snv, su_snv, msk_snv])
cnas =  pd.concat([tcga_cn, su_cn, msk_cn])

# Filter for desired genes
snvs = snvs[snvs['Gene'].isin(genes)]
cnas = cnas[cnas['Gene'].isin(cn_genes)]

# Get rid of percentage signs
snvs['Freq'] = snvs['Freq'].str.split("%").str[0].astype(float)
cnas['Freq'] = cnas['Freq'].str.split("%").str[0].astype(float)

m1_snvs = m1_snvs[m1_snvs['GENE'].isin(genes)]
m1_cnvs = m1_cnvs[m1_cnvs['GENE'].isin(cn_genes)]

# Filter for changes of interest
m1_cnvs = m1_cnvs[(m1_cnvs['Copy_num'] == 2) | (m1_cnvs['Copy_num'] == -2)]
m1_snvs = m1_snvs[m1_snvs["EFFECT"].str.contains("Frameshift|Stopgain|Missense|Splice|indel")]

# Drop duplicates
m1_snvs = m1_snvs.drop_duplicates(subset=["Patient ID", "GENE"])
m1_snvs["mutation_count"] = m1_snvs.groupby(["GENE"])["GENE"].transform("count")
m1_snvs = m1_snvs.drop_duplicates(subset=["GENE"] )

# Get frequency
m1_snvs['Mut_perc'] = m1_snvs['mutation_count']/pt_num * 100


# CNV formatting
m1_cnvs.loc[m1_cnvs["Copy_num"]==-2, "CN_type" ] = "HOMDEL"
m1_cnvs.loc[m1_cnvs["Copy_num"]== 2, "CN_type" ] = "AMP"

m1_cnvs = m1_cnvs.drop_duplicates(subset=["Patient ID", "GENE", "CN_type"])

# Get frequencies
m1_cnvs["cnv_type_count"] = m1_cnvs.groupby(["GENE", "CN_type"])["CN_type"].transform("count")
m1_cnvs = m1_cnvs.drop_duplicates(subset=["GENE", "CN_type"] )

m1_cnvs["perc"] = m1_cnvs["cnv_type_count"]/pt_num * 100

# Format to allow merge
m1_cnvs = m1_cnvs[['GENE', 'CN_type', 'perc', 'Cohort', 'cnv_type_count']]
m1_cnvs.columns = ['Gene', 'CN_type', 'Freq', 'Cohort', 'Count']
m1_cnvs['Profiled Samples'] = pt_num

m1_snvs = m1_snvs[['GENE', 'Mut_perc', 'Cohort','mutation_count']]
m1_snvs.columns = ['Gene', 'Freq', 'Cohort', 'Count']
m1_snvs['Profiled Samples'] = pt_num

m1_snvs['CN_type'] = ''

snvs = snvs[['Gene', 'Freq', 'Cohort','#','Profiled Samples']]
snvs.columns = ['Gene', 'Freq', 'Cohort', 'Count', 'Profiled Samples']
snvs['CN_type'] = ''

cnas = cnas[['Gene','CNA', 'Freq', 'Cohort','Profiled Samples', '#']]
cnas.columns = ['Gene','CN_type', 'Freq', 'Cohort','Count','Profiled Samples']

# Concat cohorts with m1 data
snvs = pd.concat([snvs, m1_snvs], sort=False)
cnas = pd.concat([cnas, m1_cnvs], sort=False)

# Get rid of changes in unexpected directions
cnas = cnas[~((cnas["Gene"].isin(ts)) & (cnas["CN_type"] == "AMP"))]


# Change gene names to reflect mut vs CN
snvs['Gene'] = snvs['Gene'] + '\nmut'
cnas['Gene'] = cnas['Gene'].apply(lambda x: x+'\ndeep del' if x in ts else x+'\namp')

# Everything into one dataframe, both CN and SNVs
final = pd.concat([snvs, cnas], sort=False)

# Sort and order
final["cohort_order"] = final["Cohort"].map(cohort_dict)
final["gene_order"] = final["Gene"].map(gene_dict)

final = final.sort_values(by=["gene_order", "cohort_order"], ascending=True)
final = final.drop_duplicates()

# =============================================================================
# Figure
# =============================================================================
# fig = plt.figure(figsize=(5.5,2.75))
# ax1 = fig.add_subplot()

color_pal = {'M1RP de novo CSPC':'#333399', 'SU2C mCRPC':'#009933', 'MSK mCSPC':'#66CC66', 'Primary TCGA PanCancer Atlas':'#CCCCCC'}

fg = sns.catplot(data=final, kind="bar", x='Gene', y="Freq", hue="Cohort", legend=True, legend_out = False, palette=color_pal)

fig = fg.fig
fig.set_size_inches(5.5,2)
ax1 = fg.ax

# =============================================================================
# Set up statistics
# =============================================================================

alts = final['Gene'].unique().tolist()

gene_x = dict(zip(alts,np.arange(0,len(alts),1)))
cohort_x = {'M1RP de novo CSPC':-0.3, 'Primary TCGA PanCancer Atlas':-0.1, 'MSK mCSPC':0.1, 'SU2C mCRPC':0.3}

final['x'] = final.apply(lambda x: gene_x.get(x['Gene'])+cohort_x.get(x['Cohort']), axis = 1)

for index, row in final.iterrows():
    ci_low, ci_upp = proportion.proportion_confint(row['Count'], row['Profiled Samples'])
    err = row['Freq'] - ci_low*100
    ax1.errorbar(x = row['x'], y = row['Freq'], yerr = err, lw = 0.8, elinewidth = 0.5, ecolor = 'k')
    

# =============================================================================
# Aethetics 
# =============================================================================

ax1.set_ylabel("Proportion with alteration", size=6)
ax1.set_xlabel(None)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment="center")


for ax in [ax1]:
    ax.margins(x=0.023)
    ax.legend(loc='upper right')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    for side in ax.spines.keys():
        ax.spines[side].set_linewidth(0.5)
        
plt.legend(loc=(0.06,0.74))
        
fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/cohort_comparison.pdf',  dpi = 600)
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/cohort_comparison.png',  dpi = 600)

