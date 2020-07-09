"""
Created on Wed Jul  8 9:49:03 2020

@author: echen
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
from matplotlib import gridspec
import matplotlib.lines as mlines

pb_avil = ['M1RP_ID1','M1RP_ID10','M1RP_ID12','M1RP_ID14','M1RP_ID15','M1RP_ID17','M1RP_ID18','M1RP_ID19','M1RP_ID21','M1RP_ID6']
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']
hypermutated = ['M1RP_ID8']

tc = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022")
sample_summary = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0")
mutations = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1RP_mutations.tsv", sep='\t')
mutations_m1b = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1B_mutations.tsv", sep='\t')
sample_summary.columns = sample_summary.iloc[0]
sample_summary = sample_summary.drop(sample_summary.index[0])
sample_summary = sample_summary[['Cohort','Sample ID','GGG','Targeted sequencing']]
sample_summary = sample_summary[sample_summary['Targeted sequencing'].str.contains("Completed",na=False)]
sample_summary = sample_summary.drop('Targeted sequencing',axis=1)

genes = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\copy number\final melted cna files\M1RP_cfdna_cna.tsv", sep='\t')
genes = genes.drop_duplicates(['GENE'],keep='first')
genes = genes.reset_index(drop=True)
gene_order = genes['GENE'].tolist()

def meltBet(path, cohort):
    df = pd.read_excel(path, index_col=None)
    df = pd.melt(df, id_vars=['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES'])
    df.rename(columns={'value': 'Allele_frequency'}, inplace=True)
    df.rename(columns={'variable': 'Sample ID'}, inplace=True)
    df['Read_depth'] = df['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[1]
    df = df[~df['Read_depth'].str.contains("*", regex=False)]
    df['Read_depth'] = df['Read_depth'].replace('\(','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\)','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\*','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace('\*','', regex=True)
    df['Read_depth'] = df['Read_depth'].replace("\[[^]]*\]",'', regex=True)
    df['Allele_frequency'] = df['Allele_frequency'].str.split(pat='%', n=-1, expand=False).str[0]
    df['Cohort'] = cohort
    df['Patient ID'] = df['Sample ID'].str.split('_').str[1]
    df = df[['Cohort','Patient ID','Sample ID', 'CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'Allele_frequency', 'Read_depth','NOTES']]
    df[['Read_depth','Allele_frequency']] = df[['Read_depth', 'Allele_frequency']].apply(pd.to_numeric)
    df['Allele_frequency'] = df['Allele_frequency'] / 100.
    return df;

"""
Getting uncalled mutations in M1RP and append them to called mutations
"""
muts_uncalled = meltBet('C:/Users/Emilia Chen/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/M1RP_betastasis_all.xlsx','M1RP')
muts_uncalled = muts_uncalled[muts_uncalled['Allele_frequency'] >= 0.01]
muts_uncalled['Mutant reads'] = muts_uncalled['Allele_frequency'] * muts_uncalled['Read_depth']
muts_uncalled = muts_uncalled[muts_uncalled['Mutant reads'] >= 3]
mutations = mutations.set_index(unique_mutation_columns)
muts_uncalled = muts_uncalled.set_index(unique_mutation_columns)
muts_uncalled = muts_uncalled[muts_uncalled.index.isin(mutations.index)]

mutations = mutations.reset_index()
muts_uncalled = muts_uncalled.reset_index()

muts_uncalled = muts_uncalled[['Patient ID','GENE','EFFECT','POSITION','Sample ID','Allele_frequency']]
mutations = mutations[['Patient ID','GENE','EFFECT','POSITION','Sample ID','Allele_frequency']]
mutation_combined = mutations.append(muts_uncalled,ignore_index=True)
mutation_combined['Cohort'] = mutation_combined['Sample ID'].str.split('_').str[0]

"""
Getting uncalled mutations in M1B and append them to called mutations
"""
muts_uncalled_m1b = meltBet('C:/Users/Emilia Chen/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/M1B_betastasis_all.xlsx','M1B')
muts_uncalled_m1b = muts_uncalled_m1b[muts_uncalled_m1b['Allele_frequency'] >= 0.01]
muts_uncalled_m1b['Mutant reads'] = muts_uncalled_m1b['Allele_frequency'] * muts_uncalled_m1b['Read_depth']
muts_uncalled_m1b = muts_uncalled_m1b[muts_uncalled_m1b['Mutant reads'] >= 3]

mutations_m1b = mutations_m1b.set_index(unique_mutation_columns)
muts_uncalled_m1b = muts_uncalled_m1b.set_index(unique_mutation_columns)
muts_uncalled_m1b = muts_uncalled_m1b[muts_uncalled_m1b.index.isin(mutations_m1b.index)]
mutations_m1b = mutations_m1b.reset_index()
muts_uncalled_m1b = muts_uncalled_m1b.reset_index()

muts_uncalled_m1b = muts_uncalled_m1b[['Patient ID','GENE','EFFECT','POSITION','Sample ID','Allele_frequency']]
mutations_m1b = mutations_m1b[['Patient ID','GENE','EFFECT','POSITION','Sample ID','Allele_frequency']]
mutation_combined_m1b = mutations_m1b.append(muts_uncalled_m1b,ignore_index=True)
mutation_combined_m1b['Cohort'] = mutation_combined_m1b['Sample ID'].str.split('_').str[0]

"""
Combining M1RP mutations and M1B mutations
"""
mutation_combined = mutation_combined.append(mutation_combined_m1b,ignore_index=True)
"""
Filter out hypermutated patients and silent mutations
"""
mutation_combined['Cohort'] = mutation_combined['Cohort'] + "_" + mutation_combined['Patient ID']
mutation_combined = mutation_combined[~mutation_combined['Cohort'].isin(hypermutated)]
mutation_combined = mutation_combined[mutation_combined['GENE'].isin(gene_order)]
silent = ["Intronic", "3'-UTR", "5'-UTR", "Upstream", "Intergenic", "Downstream"]
for index, row in mutation_combined.iterrows():
    mutation = mutation_combined.at[index, 'EFFECT']
    if mutation in silent or "Synonymous" in mutation:
        mutation_combined = mutation_combined.drop(index)
mutation_combined = mutation_combined.reset_index(drop=True)
mutation_combined = mutation_combined.drop('Patient ID',axis = 1)
"""
Restricting mutations to those from primary samples
"""
for index, row in mutation_combined.iterrows():
    sample = mutation_combined.at[index, 'Sample ID']
    if 'cfDNA' in sample:
        mutation_combined.at[index, 'Sample Type'] = 'cfDNA'
    elif 'MB' in sample:
        mutation_combined.at[index, 'Sample Type'] = 'MB'
    elif 'PB' in sample:
        mutation_combined.at[index, 'Sample Type'] = 'PB'
    elif 'MLN' in sample:
        mutation_combined.at[index, 'Sample Type'] = 'MLN'
    elif 'NL' in sample:
        mutation_combined.at[index, 'Sample Type'] = 'NL'
    else:
        mutation_combined.at[index, 'Sample Type'] = 'RP'
    if 'merged' in sample:
        names = mutation_combined.at[index, 'Sample ID'].split("_")
        mutation_combined.at[index, 'Sample ID'] =  names[0]+"_"+names[1]+"_"+names[2]
mutation_combined = mutation_combined[mutation_combined['Sample Type'].isin(['RP','PB'])]
"""
Reading GGG from sample summaries and restricting to primary samples
"""
sample_summary['Sample ID'] = sample_summary['Sample ID'].fillna('non')
for index, row in sample_summary.iterrows():
    sample = sample_summary.at[index, 'Sample ID']
    if 'cfDNA' in sample:
        sample_summary.at[index, 'Sample Type'] = 'cfDNA'
    elif 'MB' in sample:
        sample_summary.at[index, 'Sample Type'] = 'MB'
    elif 'PB' in sample:
        sample_summary.at[index, 'Sample Type'] = 'PB'
    elif 'MLN' in sample:
        sample_summary.at[index, 'Sample Type'] = 'MLN'
    elif 'NL' in sample:
        sample_summary.at[index, 'Sample Type'] = 'NL'
    elif 'gDNA' in sample:
        sample_summary.at[index, 'Sample Type'] = 'gDNA'
    elif 'non' in sample:
        sample_summary = sample_summary.drop(index)
    else:
        sample_summary.at[index, 'Sample Type'] = 'RP'
sample_summary = sample_summary[sample_summary['Sample Type'].isin(['PB','RP'])]
sample_summary['GGG'] = sample_summary['GGG'].fillna(0)
sample_summary['GGG'] = sample_summary['GGG'].astype(int)
sample_summary = sample_summary[sample_summary['GGG']>0]
sample_summary = sample_summary.reset_index(drop=True)
sample_summary['Cohort'] = sample_summary['Cohort'] + "_" + sample_summary['Sample ID'].str.split("_").str[1] 

"""
Merging mutations to sample sammaries
"""

merged = sample_summary.merge(mutation_combined,how='outer',on=['Sample ID','Sample Type','Cohort'])

"""
Reading TC and merge it to mutation_sample_summary table
"""
tc['Cohort'] = tc['Cohort']+"_"+tc['Patient ID']
for index, row in tc.iterrows():
    sample = tc.at[index, 'Sample ID']
    if 'cfDNA' in sample:
        tc.at[index, 'Sample Type'] = 'cfDNA'
    elif 'MB' in sample:
        tc.at[index, 'Sample Type'] = 'MB'
    elif 'PB' in sample:
        tc.at[index, 'Sample Type'] = 'PB'
    elif 'MLN' in sample:
        tc.at[index, 'Sample Type'] = 'MLN'
    elif 'NL' in sample:
        tc.at[index, 'Sample Type'] = 'NL'
    else:
        tc.at[index, 'Sample Type'] = 'RP'
    if 'merged' in sample:
        names = sample.split("_")
        tc.at[index, 'Sample ID'] = names[0] +"_" + names[1] + "_" + names[2]
tc = tc[tc['Sample Type'].isin(['RP','PB'])]
merged_tc = merged.merge(tc,how='left',on=['Cohort','Sample ID','Sample Type'])

merged_tc['Final tNGS_TC'] = merged_tc['Final tNGS_TC'].fillna(0)
merged_tc['Final tNGS_TC'] = merged_tc['Final tNGS_TC'].astype(float)
merged_tc = merged_tc[merged_tc['Final tNGS_TC']>=1]
merged_tc = merged_tc.sort_values(by=['GGG','Final tNGS_TC'], ascending=[False, False])
merged_tc = merged_tc.reset_index(drop = True)

"""
Set up new table for mutation frequencies
"""
patients = set(merged_tc['Cohort'])
patients = list(patients)
df_patient = pd.DataFrame(patients,columns=['Patient ID']) 
df_patient = df_patient.set_index(['Patient ID'])
df_patient_other = pd.DataFrame([],columns=['Patient ID'])
df_patient_other = df_patient_other.set_index(['Patient ID'])
sample_dict = {}
for p in patients:
    df_one_tc = tc[tc['Cohort'].isin([p])]
    df_one_tc['Final tNGS_TC'] = df_one_tc['Final tNGS_TC'].fillna(0)
    df_one_tc['Final tNGS_TC'] = df_one_tc['Final tNGS_TC'].astype(float)
    df_one_tc = df_one_tc[df_one_tc['Final tNGS_TC']>0]
    sample_dict[p] = len(set(df_one_tc['Sample ID']))
    
def mutation_frequency(df, p, gene):
    df_one_p = df[df['Cohort'].isin([p])]
    df_one_p = df_one_p.sort_values(by=['GGG','Final tNGS_TC'], ascending=[False,False])
    rep = ''
    count = len(df_one_p[df_one_p['GENE'].isin([gene])].drop_duplicates(subset=['Sample ID','GENE'],keep='first'))
    for index, row in df_one_p.iterrows():
        if df_one_p.at[index,'Sample Type'] == 'PB' and rep == '':
            rep = df_one_p.at[index,'Sample ID']
        if df_one_p.at[index,'Cohort'] not in pb_avil and rep == '':
            rep = df_one_p.at[index,'Sample ID']
        if df_one_p.at[index,'Sample ID'] == rep and df_one_p.at[index,'GENE'] == gene:
            p_id = rep.split('_')[0]+ "_" +rep.split('_')[1]
            df_patient.at[p_id,'Rep has ' + gene + ' mut'] = 'Yes'
    per = count/sample_dict[p] * 100
    if sample_dict[p] > 1:
        p_id = rep.split('_')[0]+ "_" + rep.split('_')[1]
        print(p_id + ":" + rep + ":" + str(per))
        df_patient.at[p_id,'% '+ gene +' mutation'] = per
        df_patient.at[p_id, gene +' rep'] = rep

for p in patients:
    mutation_frequency(merged_tc, p, 'TP53')
    mutation_frequency(merged_tc, p, 'FOXA1')
    mutation_frequency(merged_tc, p, 'PTEN')
    mutation_frequency(merged_tc, p, 'SPOP')
    mutation_frequency(merged_tc, p, 'APC')

checked = ['TP53','PTEN','FOXA1','SPOP','APC']
for p in patients:
    df_one_p = merged_tc[merged_tc['Cohort'].isin([p])]
    df_one_tc = tc[tc['Cohort'].isin([p])]
    df_one_tc['Final tNGS_TC'] = df_one_tc['Final tNGS_TC'].fillna(0)
    df_one_tc['Final tNGS_TC'] = df_one_tc['Final tNGS_TC'].astype(float)
    df_one_tc = df_one_tc[df_one_tc['Final tNGS_TC']>0]
    df_one_p = df_one_p.sort_values(by=['GGG','Final tNGS_TC'], ascending=[False,False])
    num_sample = len(set(df_one_tc['Sample ID']))
    rep = ''
    for index, row in df_one_p.iterrows():
        if df_one_p.at[index,'Sample Type'] == 'PB' and rep == '':
            rep = df_one_p.at[index,'Sample ID']
        if df_one_p.at[index,'Cohort'] not in pb_avil and rep == '':
            rep = df_one_p.at[index,'Sample ID']
    mut_rep = df_one_p[df_one_p['Sample ID'].isin([rep])]['GENE'].tolist()
    for mut in mut_rep:
        if mut not in checked and mut in gene_order:
            count = len(df_one_p[df_one_p['GENE'].isin([mut])].drop_duplicates(subset=['Sample ID','GENE'],keep='first'))
            per = count/num_sample * 100
            if num_sample > 1:
                p_id = rep.split('_')[0]+ "_" + rep.split('_')[1] + str(mut)
                print(p_id + ":" + rep + ":" + str(per))
                df_patient_other.at[p_id,'Rep has other mut'] = 'Yes'
                df_patient_other.at[p_id,'% other mutation'] = per
    allmut = list(set(df_one_p['GENE']))
    for mut in allmut:
        if mut not in checked and mut in gene_order and mut not in mut_rep:
            count = len(df_one_p[df_one_p['GENE'].isin([mut])].drop_duplicates(subset=['Sample ID','GENE'],keep='first'))
            per = count/num_sample * 100
            p_id = rep.split('_')[0]+ "_" + rep.split('_')[1] + str(mut)
            print(p_id + ":" + rep + ":" + str(per))
            df_patient_other.at[p_id,'Rep has other mut'] = 'No'
            df_patient_other.at[p_id,'% other mutation'] = per

df_patient['Rep has TP53 mut'] = df_patient['Rep has TP53 mut'].fillna('No')
df_patient['Rep has FOXA1 mut'] = df_patient['Rep has FOXA1 mut'].fillna('No')
df_patient['Rep has PTEN mut'] = df_patient['Rep has PTEN mut'].fillna('No')
df_patient['Rep has SPOP mut'] = df_patient['Rep has SPOP mut'].fillna('No')
df_patient['Rep has APC mut'] = df_patient['Rep has APC mut'].fillna('No')
df_patient_other['Rep has other mut'] = df_patient_other['Rep has other mut'].fillna('No')

df_patient = df_patient.reset_index()
df_patient_other = df_patient_other.reset_index()

df_patient['key'] = df_patient['Patient ID'].str.split('_').str[0]
df_patient_other['key'] = df_patient_other['Patient ID'].str.split('_').str[0]

df_patient_tp53 = df_patient.copy()
df_patient_tp53 = df_patient_tp53[df_patient_tp53['% TP53 mutation']>0]
df_patient_foxa1 = df_patient.copy()
df_patient_foxa1 = df_patient_foxa1[df_patient_foxa1['% FOXA1 mutation']>0]
df_patient_pten = df_patient.copy()
df_patient_pten = df_patient_pten[df_patient_pten['% PTEN mutation']>0]
df_patient_spop = df_patient.copy()
df_patient_spop = df_patient_spop[df_patient_spop['% SPOP mutation']>0]
df_patient_apc = df_patient.copy()
df_patient_apc = df_patient_apc[df_patient_apc['% APC mutation']>0]

"""
Plotting
"""
fig = plt.figure(figsize=(14, 4))
gs = gridspec.GridSpec(1,7,width_ratios=[1,1,1,1,1,1,0.5])
gs.update(wspace=0, hspace=0.1, bottom = 0.25, left = 0.15, right = 0.99, top=0.95)
ax1 =  fig.add_subplot(gs[0])
ax2 =  fig.add_subplot(gs[1])
ax3 =  fig.add_subplot(gs[2])
ax4 =  fig.add_subplot(gs[3])
ax5 =  fig.add_subplot(gs[4])
ax7 =  fig.add_subplot(gs[5])

l = 8
sns.swarmplot(x = 'key', y = '% TP53 mutation', data = df_patient_tp53, ax = ax1, size = 5, zorder = 0, hue="Rep has TP53 mut",hue_order=["Yes","No"],order=['M1RP','M1B'], palette=['r','b'])
n_m1rp = len(df_patient_tp53[df_patient_tp53['key'].isin(['M1RP'])])
n_m1b = len(df_patient_tp53[df_patient_tp53['key'].isin(['M1B'])])
ax1.set_xticklabels(['M1RP\nn='+str(n_m1rp),'M1B\nn='+str(n_m1b)])
ax1.set_xlabel('TP53', fontsize=14,labelpad=l)
ax1.legend_.remove()

sns.swarmplot(x = 'key', y = '% FOXA1 mutation', data = df_patient_foxa1, ax = ax2, size = 5, zorder = 0, hue="Rep has FOXA1 mut",hue_order=["Yes","No"],order=['M1RP','M1B'],palette=['r','b'])

n_m1rp = len(df_patient_foxa1[df_patient_foxa1['key'].isin(['M1RP'])])
n_m1b = len(df_patient_foxa1[df_patient_foxa1['key'].isin(['M1B'])])
ax2.set_xticklabels(['M1RP\nn='+str(n_m1rp),'M1B\nn='+str(n_m1b)])
ax2.set_xlabel('FOXA1', fontsize=14,labelpad=l)
ax2.legend_.remove()

sns.swarmplot(x = 'key', y = '% PTEN mutation', data = df_patient_pten, ax = ax3, size = 5, zorder = 0, hue="Rep has PTEN mut",hue_order=["Yes","No"],order=['M1RP','M1B'],palette=['r','b'])

n_m1rp = len(df_patient_pten[df_patient_pten['key'].isin(['M1RP'])])
n_m1b = len(df_patient_pten[df_patient_pten['key'].isin(['M1B'])])
ax3.set_xticklabels(['M1RP\nn='+str(n_m1rp),'M1B\nn='+str(n_m1b)])
ax3.set_xlabel('PTEN', fontsize=14,labelpad=l)
ax3.legend_.remove()

sns.swarmplot(x = 'key', y = '% SPOP mutation', data = df_patient_spop, ax = ax4, size = 5, zorder = 0, hue="Rep has SPOP mut",hue_order=["Yes","No"],order=['M1RP','M1B'],palette=['r','b'])

n_m1rp = len(df_patient_spop[df_patient_spop['key'].isin(['M1RP'])])
n_m1b = len(df_patient_spop[df_patient_spop['key'].isin(['M1B'])])
ax4.set_xticklabels(['M1RP\nn='+str(n_m1rp),'M1B\nn='+str(n_m1b)])
ax4.set_xlabel('SPOP', fontsize=14,labelpad=l)
ax4.legend_.remove()

sns.swarmplot(x = 'key', y = '% APC mutation', data = df_patient_apc, ax = ax5, size = 5, zorder = 0, hue="Rep has APC mut",hue_order=["Yes","No"],order=['M1RP','M1B'],palette=['r','b'])

n_m1rp = len(df_patient_apc[df_patient_apc['key'].isin(['M1RP'])])
n_m1b = len(df_patient_apc[df_patient_apc['key'].isin(['M1B'])])
ax5.set_xticklabels(['M1RP\nn='+str(n_m1rp),'M1B\nn='+str(n_m1b)])
ax5.set_xlabel('APC', fontsize=14,labelpad=l)
ax5.legend_.remove()

sns.swarmplot(x = 'key', y = '% other mutation', data = df_patient_other, ax = ax7, size = 5, zorder = 0, hue="Rep has other mut",hue_order=["Yes","No"],order=['M1RP','M1B'],palette=['r','b'])

ax7.set_xticklabels(['M1RP\n','M1B\n'])
ax7.set_xlabel('Others', fontsize=14,labelpad=l)
ax7.legend_.remove()

def turnOffSpines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis = 'y',labelleft = False, left = False)
    ax.set_ylabel('')
    
for i,ax in enumerate([ax1,ax2,ax3,ax4,ax5,ax7]):
    if i != 0:
        turnOffSpines(ax)
    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('% primary samples per patient\nhaving the mutation',fontsize=12)
    if i%2 == 0:
        ax.set_facecolor('#f0f0f0')

    ax.set_ylim(-2,102)
    ax.set_yticks([0,20,40,60,80,100])
    ax.tick_params(axis = 'x',labelsize = 12, length=0, pad=12)
    ax.tick_params(axis = 'y',labelsize = 12, length=1, pad=6, width=0.5)
    

ax6 =  fig.add_subplot(gs[6])
labels = ['Present','Absent']
handles = [mlines.Line2D([], [], color='r', marker='o', lw=0, markersize=5),
           mlines.Line2D([], [], color='b', marker='o', lw=0, markersize=5)]
ax6.legend(handles, labels, fontsize =12, title_fontsize=12,title='Mutation status in the\nrepresentative sample', loc = 'center left',handlelength = 0.8, frameon = False)
ax6.spines['bottom'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)
fig.savefig('mutation analysis.pdf',bbox_extra_artists=(ax6,), bbox_inches='tight')