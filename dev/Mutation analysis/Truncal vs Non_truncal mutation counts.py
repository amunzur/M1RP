"""
Created on Thur Jul 9 11:09:03 2020

@author: echen
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
from matplotlib import gridspec

unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']
hypermutated = ['M1RP_ID8']

tc = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022")
mutations = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1RP_mutations.tsv", sep='\t')
mutations_m1b = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1B_mutations.tsv", sep='\t')
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
Restricting mutations to those from primary samples (or comment out the last line to analyze all samples)
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
mutation_combined = mutation_combined[mutation_combined['Sample Type'].isin(['PB','RP'])]
"""
Restricting TC to those from primary samples (or comment out the last line to analyze all samples)
"""
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[['Cohort','Patient ID', 'Sample ID', 'Final tNGS_TC']]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)
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
tc = tc[tc['Sample Type'].isin(['PB','RP'])]

"""
Merge TC to mutation table
"""
merged_tc = mutation_combined.merge(tc,how='outer',on=['Cohort','Sample ID','Sample Type'])

merged_tc['Final tNGS_TC'] = merged_tc['Final tNGS_TC'].fillna(0)
merged_tc['Final tNGS_TC'] = merged_tc['Final tNGS_TC'].astype(float)
merged_tc = merged_tc[merged_tc['Final tNGS_TC']>0.01]

"""
Create a new table to summarize mutation frequency
"""
patients = set(merged_tc['Cohort'])
patients = list(patients)
df_patient = pd.DataFrame(patients,columns=['Patient ID']) 
df_patient['Cohort'] = df_patient['Patient ID'].str.split("_").str[0]
df_patient = df_patient.set_index(['Patient ID'])

"""
Calculating mutation frequency (per mutation per patient)
"""
df_genes = pd.DataFrame(gene_order,columns=['GENE'])
df_genes = df_genes.set_index(['GENE'])
df_genes['M1RP Truncal'] = 0
df_genes['M1RP Non-truncal'] = 0
df_genes['M1B Truncal'] = 0
df_genes['M1B Non-truncal'] = 0

for p in patients:
    if True:
        cohort = p.split("_")[0]
        df_one_p = merged_tc[merged_tc['Cohort'].isin([p])]
        df_one_p['Identifier'] = df_one_p['GENE']+"+"+df_one_p['EFFECT']
        all_mut = set(df_one_p['Identifier'].dropna())
        for mut in all_mut:
            gene = mut.split("+")[0]
            effect = mut.split("+")[1]
            count = len(df_one_p.query('GENE == "' + gene +'" and EFFECT== "' + effect + '"').drop_duplicates(subset=['Sample ID'],keep='first'))
            num_sumple = len(df_one_p.drop_duplicates(subset=['Sample ID'],keep='first'))
            mut_freq = count/num_sumple
#             print(gene+" "+str(count))
            if mut_freq >= 0.8 and num_sumple > 1:
#                 print("T " + p + " " + str(mut_freq))
                df_genes.at[gene, cohort + ' Truncal'] = df_genes.at[gene, cohort + ' Truncal'] + 1
            elif mut_freq > 0 and num_sumple > 1:
#                 print("NT " + p + " " + str(mut_freq))
                df_genes.at[gene, cohort + ' Non-truncal'] = df_genes.at[gene, cohort + ' Non-truncal'] + 1
df_genes['Sum'] = df_genes['M1RP Truncal'] + df_genes['M1RP Non-truncal'] + df_genes['M1B Truncal'] + df_genes['M1B Non-truncal']
df_genes = df_genes.sort_values(by=['Sum'],ascending=False)
df_genes_filtered = df_genes[df_genes['Sum']>0]

"""
Plotting
"""
genes = df_genes_filtered.index.values
m1rp_trunk = np.array(df_genes_filtered['M1RP Truncal'])
m1rp_non_trunk = np.array(df_genes_filtered['M1RP Non-truncal'])
m1b_trunk = np.array(df_genes_filtered['M1B Truncal'])
m1b_non_trunk = np.array(df_genes_filtered['M1B Non-truncal'])

barWidth = 0.25
m1rp_pos = [x for x, _ in enumerate(genes)]
m1b_pos = [x + barWidth+0.1 for x in m1rp_pos]
tick_pos = [x + 0.18 for x in m1rp_pos]
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams.update({'hatch.color': 'white'})

fig = plt.figure(figsize=(16, 6))
ax = fig.add_subplot(111)
plt.bar(m1rp_pos, m1rp_trunk, width=barWidth, label='M1RP Truncal', color='#9C133A')
plt.bar(m1rp_pos, m1rp_non_trunk, width=barWidth, label='M1RP Non-truncal', color='#3b50b2', bottom=m1rp_trunk)
plt.bar(m1b_pos, m1b_trunk, width=barWidth, color='#9C133A',hatch='////', label='M1B Truncal')
plt.bar(m1b_pos, m1b_non_trunk, width=barWidth, color='#3b50b2', hatch='////', label='M1B Non-truncal',bottom=m1b_trunk)
ax.tick_params(bottom = False)
gene_label = [r'$\it{%s}$' % label for label in genes]
plt.xticks(m1rp_pos, gene_label,rotation = 90,fontsize=16)
ax.set_ylabel('Mutation count',fontsize=16)
ax.set_title('Mutation count in all samples',fontsize=16)
# plt.annotate('N (M1RP) = 353\nN (M1B) = 181', xy=(0, 1), xytext=(700, -112), fontsize=16, va='top', xycoords='axes fraction', textcoords='offset points')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(fontsize=16)

plt.show()
fig.savefig('Primary sample mutation count.pdf',bbox_inches = "tight")