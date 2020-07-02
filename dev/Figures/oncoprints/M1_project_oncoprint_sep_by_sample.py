"""
Created on Wed Jun 18 20:00:28 2020

M1 project oncorpint, separated by sample, include unable to call mutation

@author: echen
"""
#TODO: figure out how to add dahsed line on squares

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from datetime import datetime as dt

"""
Setting parameters
"""
# I. Formating
# 1. Font & Size
axis_label_font_size = 12
tick_label_font_size = 12
legend_font_size = 12

# 2. Axis names

# II. Artistics
# 1. cna
cnv_color_dict = {
    '-2' : '#3F60AC',
    '-1' : '#9CC5E9',
    '0'  : '#E6E7E8',
    '1'  : '#F59496',
    '2'  : '#EE2D24',
    'unable': 'unable_#E6E7E8'
}

# 2. Somatic mutations
mutation_color_dict = {
    'missense': '#79B443',
    'non-frameshift indel': '#a9a9a9',
    'non-frameshift': '#a9a9a9',
    'frameshift': '#FFC907',
    'stopgain': '#FFC907',
    'splice site': '#FFC907',
    'splice':'#FFC907',
    'intronic': '#605857',
    'upstream': '#605857',
    'downstream': '#605857',
    "3'-utr": '#605857',
    "5'-utr": '#605857',
    'synonymous': '#605857',
    'unable': 'unable_',
    'nan':'nan'
}


# 3. Figure size
width = 10
height = 16

# III. Data
cohort = 'M1RP'
patientID = 'ID22'
single_loss_threshold = 37.5
double_loss_threshold = 19.5
single_gain_threshold = 46
amplification_threshold = 23.5
genes = []
gene_order = []
num_genes = 73
num_genes = num_genes + 1 # one for making it into column index
 
tc = pd.read_csv("https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022")
mutations = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1RP_mutations.tsv", sep='\t')
cfDNA_cna = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\copy number\final melted cna files\M1RP_cfdna_cna.tsv", sep='\t')
ffpe_cna = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\copy number\final melted cna files\M1RP_FFPE_cna.tsv", sep='\t')
date = pd.read_excel(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\Clinical data tables\PB and RP dates.xlsx",index_col=0)
betastasis = pd.read_excel(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\betastasis\M1RP_betastasis_all.xlsx",index_col=0)
genes = cfDNA_cna.copy()
mutations_copy = mutations.copy() # need it later, declare here first
"""
Helper functions
"""

# Helper for scatter plots
def plot_scatter(ax, markstyle, ul, lr, cell):
    if 'unable' in cell:
        cellcol = '#FFFFFF'
        edgecol = 'k'
    elif 'non-independent' in cell:
        cellcol = cell.split("_")[1]
        print(cell)
        ax.scatter(ul, lr, marker = 'o', s = 35, color = cellcol, edgecolors = 'none', zorder = 20, linewidth = 0.4) 
        return
    else:
        cellcol = cell
        edgecol = 'none'
    ax.scatter(ul, lr, marker = markstyle, s = 35, color = cellcol, edgecolors = edgecol, zorder = 20, linewidth = 0.4)

# Helper for plotting oncoprint
def ploting_oncoprint(tc_cna_table, mut_table, first_plot, tc_ax, onco_axis):
    tc_cna_table['Final tNGS_TC'] = tc_cna_table['Final tNGS_TC'].astype(float)
    if (first_plot):
        tc_ax.set_ylabel('TC%', fontsize=axis_label_font_size)
        tc_ax.set_yticks([0, 20, 40, 60, 80, 100])
        tc_ax.set_yticklabels([0, 20, 40, 60, 80, 100],fontsize=tick_label_font_size)
    else:
        tc_ax.yaxis.set_ticklabels([])
        tc_ax.yaxis.set_ticks_position('none')
        tc_ax.spines['left'].set_visible(False)
    
    tc_ax.bar(tc_cna_table['Sample ID'], tc_cna_table['Final tNGS_TC'], zorder = 10, color = '#909396')
    tc_ax.set_ylim(0,100)
    tc_ax.tick_params(bottom = False, labelbottom = False)
    tc_ax.spines['bottom'].set_zorder(100)
    tc_ax.spines['top'].set_visible(False)
    tc_ax.spines['right'].set_visible(False)
#     tc_ax.grid(axis = 'y', linestyle = 'dashed', linewidth = 0.5, color = '0.7', zorder = 0)
    if len(tc_cna_table)>0 or first_plot:
        tc_ax.axhline(y=single_loss_threshold,xmin=0,xmax=1,c="#9CC5E9",linewidth=1,zorder=20)
        tc_ax.axhline(y=double_loss_threshold,xmin=0,xmax=1,c="#3F60AC",linewidth=1,zorder=20)
        tc_ax.axhline(y=single_gain_threshold,xmin=0,xmax=1,c="#F59496",linewidth=1,zorder=20)
        tc_ax.axhline(y=amplification_threshold,xmin=0,xmax=1,c="#EE2D24",linewidth=1,zorder=20)

    for index, row in tc_cna_table.iterrows():
        for i, col in enumerate(tc_cna_table.columns.tolist()[1:num_genes]):
            if 'unable' in tc_cna_table.at[index,col]:
                onco_axis.bar(row['Sample ID'], 0.8, bottom = i, color = 'none', edgecolor= row[col].split("_")[1], zorder = 10)
            else:
                onco_axis.bar(row['Sample ID'], 0.8, bottom = i, color = row[col], edgecolor= 'none', zorder = 10)
    
    for index, row in mut_table.iterrows():
        for i, col in enumerate(mut_table.columns.tolist()[1:num_genes]):
            mut_table[col] = mut_table[col].astype(str)
            cell = mut_table.loc[index, col]
            jitter = 0
            more_than_two_mutations = 0
            # Assumption: 
            #   if there is one mut in a sample that is unable to be called due to low TC or read depth
            #   then all potential muts (>2 muts) in that sample are more likely to be non-callable
            if ':' in cell:
                cell1 = cell.split(":")[1]
                cell2 = cell.split(":")[2]
                if len(cell.split(":")) > 3:
                    jitter = 0
                    more_than_two_mutations = 1
                else:
                    jitter = 0.1
            else:
                jitter = 0
            if 'nan' == cell: continue;
            if 'no' == cell: continue;
                
            if jitter > 0:
                plot_scatter(ax=onco_axis, markstyle='s', ul=index+jitter, lr=i+0.4+jitter, cell=cell1)
                plot_scatter(ax=onco_axis, markstyle='s', ul=index-jitter, lr=i+0.4-jitter, cell=cell2)
            elif more_than_two_mutations != 0:
                plot_scatter(ax=onco_axis, markstyle='^', ul=index, lr=i+0.4, cell=cell1)
            else:
                plot_scatter(ax=onco_axis, markstyle='s', ul=index, lr=i+0.4, cell=cell)
    
    onco_axis.set_xlim(-0.6,len(mut_table)-0.4)
    onco_axis.set_ylim(-0.1, len(mut_table.columns.tolist()[1:num_genes]))
    onco_axis.set_xticklabels(mut_table['Sample ID'], rotation = 90,fontsize=tick_label_font_size)
    onco_axis.tick_params(left = False, bottom = False)
    onco_axis.spines['bottom'].set_visible(False)
    onco_axis.spines['left'].set_visible(False)
    onco_axis.spines['right'].set_visible(False)
    onco_axis.spines['top'].set_visible(False)
    onco_axis.invert_yaxis()
    if (first_plot):
        yticklabels = [r'$\it{%s}$' % label for label in mut_table.columns.tolist()[1:num_genes]]
        onco_axis.set_yticks(np.arange(0.4,len(mut_table.columns.tolist()[1:num_genes]), 1))
        onco_axis.set_yticklabels(yticklabels,fontsize=tick_label_font_size)
    else:
        onco_axis.axes.yaxis.set_ticks([])


"""
Reorganizing individual source table
"""
# I. Target genes 
# Get all genes from cfDNA_cna table and sort by chrom and start pos in ascending order
genes = genes[['GENE']]
genes = genes.drop_duplicates(['GENE'],keep='first')
genes = genes.reset_index(drop=True)
gene_order = genes['GENE'].tolist()

# II. Tumor content estimates
# Remove the first row due to the formatting of the source table
# Extract relevant columns and rows with cohort of interest
# Add a column to keep track of sample type
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[['Cohort','Patient ID', 'Sample ID', 'Final tNGS_TC']]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(str)
tc_cohort = tc[tc['Cohort'].isin([cohort])]
tc_cohort['Sample Type'] = None
for index, row in tc_cohort.iterrows():
    sample = tc_cohort.at[index, 'Sample ID']
    if 'cfDNA' in sample:
        tc_cohort.at[index, 'Sample Type'] = 'cfDNA'
    elif 'MB' in sample:
        tc_cohort.at[index, 'Sample Type'] = 'MB'
    elif 'PB' in sample:
        tc_cohort.at[index, 'Sample Type'] = 'PB'
    elif 'MLN' in sample:
        tc_cohort.at[index, 'Sample Type'] = 'MLN'
    elif 'NL' in sample:
        tc_cohort.at[index, 'Sample Type'] = 'NL'
    else:
        tc_cohort.at[index, 'Sample Type'] = 'RP'

# III. Dates
# Spliting the date table to PB and RP subtable and attach one to the back to the other, for the ease of sorting
date_PB = date[['ID ', 'Date PB']]
date_RP = date[['ID ', 'Date RP']]
date_PB['Sample Type'] = 'PB'
date_RP['Sample Type'] = 'RP'
date_PB = date_PB.rename(columns = {'Date PB':'Date'})
date_RP = date_RP.rename(columns = {'Date RP':'Date'})
date_all = date_PB.append(date_RP,ignore_index=True)
date_all['Cohort'] = date_all['ID '].str.split("_").str[0]
date_all = date_all[date_all['Cohort'].isin([cohort])]
date_all = date_all.drop('Cohort', axis=1)
date_all['Patient ID'] = date_all['ID '].str.split("_").str[1]

# IV. Mutations
# 1) Extract relevant columns,
# 2) Remove genes not in the gene_order list,
# 3) Curate gene name if multiple genes are called in one record
mutations = mutations[['Patient ID', 'Sample ID', 'GENE', 'EFFECT']]
mutations['EFFECT'] = mutations['EFFECT'].str.split(" ").str[0]
mutations['GENE'] = mutations['GENE'].astype(str)
for index, row in mutations.iterrows():
    if len(mutations.at[index,'GENE'].split(';')) > 1:
        for gene in mutations.at[index,'GENE'].split(';'):
            if gene in gene_order:
                mutations.at[index,'GENE'] = gene
    if mutations.at[index, 'GENE'] not in gene_order:
        mutations = mutations.drop(index)
mutations = mutations.reset_index(drop=True)
# 4) Merge mutations at the same gene into one record
mutations['EFFECT'] = mutations['EFFECT'].str.lower()
mutations['EFFECT'] = mutations.groupby(['Sample ID','GENE'])['EFFECT'].transform(lambda x: ':'.join(x))
mutations = mutations.drop_duplicates(subset=['Sample ID','GENE'], keep='first')
mutations = mutations.reset_index(drop=True)

# V. cna
# Extracting relevant columns and append ffpe_cna to the back of cfDNA_cna
cfDNA_cna = cfDNA_cna[['Patient ID', 'Sample ID', 'GENE', 'Copy_num']]
ffpe_cna = ffpe_cna[['Patient ID', 'Sample ID', 'GENE', 'Copy_num']]
all_cna = cfDNA_cna.append(ffpe_cna, ignore_index = True)
all_cna['Copy_num'] = all_cna['Copy_num'].astype(str)

# VI. betastasis
# 1) Curate gene names
# 2) Extract read depth from each sample for each mutation, called mutations are flagged by "*"
betastasis = betastasis.reset_index(drop=True)
betastasis['EFFECT'] = betastasis['EFFECT'].str.split(" ").str[0].str.lower()
for index,row in betastasis.iterrows():
    cell = str(betastasis.at[index,'GENE'])
    if ";" in cell:
        for gene in cell.split(";"):
            if gene in gene_order:
                betastasis.at[index,'GENE'] = gene
                cell = gene
    if cell not in gene_order:
        betastasis = betastasis.drop(index)
betastasis = betastasis.reset_index(drop=True)           
for index,row in betastasis.iterrows():
    for col in betastasis.columns.tolist()[6:]:
        cell = str(betastasis.at[index,col])
        vaf = ":" + cell.split("%")[0]
        if "*" not in cell:
            betastasis.at[index,col] = cell[cell.find("(")+1:cell.find(")")] + vaf
        else:
            betastasis.at[index,col] = cell[cell.find("(")+1:cell.find(")")] + vaf + "*"


"""
Merging date, tc, cna, mutation
"""
# I. Merging date to tumor content --> tc_cohort_merge
# If cfDNA sample, date is from Sample ID
# TODO: MLN, NL don't have date info
tc_cohort_merge = tc_cohort.merge(date_all, how='left', on=['Patient ID', 'Sample Type'])
tc_cohort_merge = tc_cohort_merge.drop('ID ',axis = 1)
for index,row in tc_cohort_merge.iterrows():
    sample = tc_cohort_merge.at[index,'Sample ID'].lower()
    if 'cfdna' in sample:
        date = sample.split("_")[3].lower()
        tc_cohort_merge.at[index, 'Date'] = dt.strptime(date, '%Y%b%d')
    elif 'MB' in sample:
        tc_cohort_merge.at[index, 'Date'] = dt.strptime('2018Oct', '%Y%b')


# II. Produece a template with only all patient and sample
#     Merge genes to the template --> template_genes
template = tc_cohort[['Cohort','Patient ID','Sample ID','Sample Type']]
genes['key'] = '0'
template['key'] = '0'
template_genes = template.merge(genes,how='left',on='key')
template_genes = template_genes.drop('key',axis=1)

# III. Merge template_genes with cna --> template_genes_cnv
template_genes_cnv = template_genes.merge(all_cna,how='left',on=['Patient ID','Sample ID','GENE'])

# IV. Merge template_genes_cnv with mutations --> template_genes_cnv_mutation
template_genes_cnv_mutation = template_genes_cnv.merge(mutations,how='left',on=['Patient ID','Sample ID','GENE'])

# V. Merge 'Copy_num' and 'EFFECT' column and they are separated by ';' 
template_genes_cnv_mutation['EFFECT'] = template_genes_cnv_mutation['EFFECT'].astype(str)
template_genes_cnv_mutation['Copy_num'] = template_genes_cnv_mutation['Copy_num'].astype(str)
for index, row in template_genes_cnv_mutation.iterrows():
    template_genes_cnv_mutation.at[index,'EFFECT'] = template_genes_cnv_mutation.at[index,'Copy_num']+ ";" + template_genes_cnv_mutation.at[index,'EFFECT']

# VI. Pivot the GENE column to the header row and rearrange the order of genes according to gene_order specified above
pivot_template_genes_cnv_mutation = template_genes_cnv_mutation.pivot(index='Sample ID', columns='GENE',values='EFFECT')
pivot_template_genes_cnv_mutation = pivot_template_genes_cnv_mutation.reset_index()
ordered = list()
ordered.append('Sample ID')
ordered.extend(gene_order)
pivot_template_genes_cnv_mutation = pivot_template_genes_cnv_mutation[ordered]

# VII. Merge the pivoted table from VI with tumor content and date --> pivot_template_genes_cnv_mutation_tc
pivot_template_genes_cnv_mutation_tc = pivot_template_genes_cnv_mutation.merge(tc_cohort_merge, how='left',on=['Sample ID'])

"""
Get individual patient info
"""
# I. Take info from pivot_template_genes_cnv_mutation_tc
#    and divide into tc_cna table and mutations table
one_patient_all_info = pivot_template_genes_cnv_mutation_tc[pivot_template_genes_cnv_mutation_tc['Patient ID'].isin([patientID])]
one_patient_all_info['Date'] =pd.to_datetime(one_patient_all_info.Date)
one_patient_all_info = one_patient_all_info.sort_values(by='Date')
one_patient_all_info = one_patient_all_info.reset_index(drop=True)

one_patient_tc_cna = one_patient_all_info.copy()
one_patient_mutations = one_patient_all_info.copy()
for col in one_patient_tc_cna.columns.tolist()[1:num_genes]:
    one_patient_tc_cna[col] = one_patient_tc_cna[col].str.split(';').str[0]

for col in one_patient_mutations.columns.tolist()[1:num_genes]:
    one_patient_mutations[col] = one_patient_mutations[col].str.split(';').str[1]

# II. Get the list of samples the patient has and the mutations called in the patient
patient_sample_list = set(one_patient_mutations['Sample ID'])

mutations_copy = mutations_copy[['Patient ID', 'Sample ID', 'GENE', 'EFFECT','POSITION']]
mutations_copy['EFFECT'] = mutations_copy['EFFECT'].str.split(" ").str[0].str.lower()
mutations_copy['POSITION'] = mutations_copy['POSITION'].astype(str)
mutations_copy['GENE'] = mutations_copy['GENE'].astype(str)
for index, row in mutations_copy.iterrows():
    if len(mutations_copy.at[index,'GENE'].split(';')) > 1:
        for gene in mutations_copy.at[index,'GENE'].split(';'):
            if gene in gene_order:
                mutations_copy.at[index,'GENE'] = gene
    if mutations_copy.at[index, 'GENE'] not in gene_order:
        mutations_copy = mutations_copy.drop(index)
mutations_copy = mutations_copy.reset_index(drop=True)
one_patient_from_mutations_copy = mutations_copy[mutations_copy['Patient ID'].isin([patientID])]
# Create a column called ' Mutation Identifier' to uniquely identify each mutation call
one_patient_from_mutations_copy['Mutation Identifier'] = one_patient_from_mutations_copy['GENE'] + ":" +one_patient_from_mutations_copy['EFFECT'] + ":" +one_patient_from_mutations_copy['POSITION']
patient_called_mutations = set(one_patient_from_mutations_copy['Mutation Identifier'])

# III. For each called mutation (in patient_called_mutations),
#        examine the read depth for each sample (from betastasis table)
#        and determine if we have enough TC to detect the mutation
tc_dict = tc_cohort.set_index('Sample ID')['Final tNGS_TC'].to_dict()
one_patient_mutations = one_patient_mutations.set_index('Sample ID', drop=False)

# We are confident with mutation calling if 
#  i)  read depth >= 100 and (read depth * TC) > 8 or ii) read depth < 100 and TC > 8% 
for gene_eff_pos in patient_called_mutations:
    gene = gene_eff_pos.split(":")[0]
    effect = gene_eff_pos.split(":")[1]
    pos = gene_eff_pos.split(":")[2]
    for sample in patient_sample_list:
        query_df = betastasis.query('GENE == "' + gene +'" and EFFECT== "' + effect +'" and POSITION== "' + pos + '"')
        query_df = query_df.reset_index(drop=True)
        if query_df.empty:
            print(sample + ":" + gene + ": " + effect + ' empty')
        elif len(query_df.index) > 1:
            print('more than one records')
        else:
            tc_sample = float(tc_dict[sample]) * 0.01
            read_depth = float(query_df.at[0, sample].split(":")[0])
            vaf = float(query_df.at[0, sample].split(":")[1].split("*")[0]) * 0.01
            mutant_read = read_depth * vaf
            eligible_read = read_depth * tc_sample
            cell = one_patient_mutations.at[sample,gene]
            if "*" not in query_df.at[0,sample] and vaf >= 0.01 and mutant_read >= 3:
                if cell == 'nan':
                    one_patient_mutations.at[sample,gene] = 'non-independent_' + effect
                    print(sample+":"+gene+":"+effect+" "+str(vaf)+" "+str(mutant_read))
                else:
                    one_patient_mutations.at[sample,gene] = cell + ":non-independent_" + effect
                    print(sample+":"+gene+":"+effect+" "+str(vaf)+" "+str(mutant_read))
            elif (tc_sample > 0) and ((read_depth < 100 and tc_sample < 0.08) or (read_depth >= 100 and eligible_read < 8)):
                if cell == 'nan':
                    one_patient_mutations.at[sample,gene] = 'unable_' + effect
                else:
                    one_patient_mutations.at[sample,gene] = cell + ":unable_" + effect

# IV. Map colors to mutation and cnv
for index,row in one_patient_mutations.iterrows():
    for col in one_patient_mutations.columns[1:num_genes]:
        cell = str(one_patient_mutations.at[index,col])
        color = ''
        if ':' in cell:
            for mut in cell.split(":"):
                if 'unable' in mut:
                    color = color + ":" + 'unable_' + mutation_color_dict[mut.split("_")[1]]
                elif 'non-independent' in mut:
                    color = color + ":" + 'non-independent_' + mutation_color_dict[mut.split("_")[1]]
                else:
                    color = color + ":" + mutation_color_dict[mut]
        elif 'unable' in cell:
            color = 'unable_' + mutation_color_dict[cell.split("_")[1]]
        elif 'non-independent' in cell:
            color = 'non-independent_' + mutation_color_dict[cell.split("_")[1]]
        else:
            color = mutation_color_dict[cell]
        one_patient_mutations.at[index,col] = color

for index,row in one_patient_tc_cna.iterrows():
    for col in one_patient_tc_cna.columns[1:num_genes]:
        cell = one_patient_tc_cna.at[index,col]
        sample = one_patient_tc_cna.at[index,'Sample ID']
        if tc_dict[sample] == '0':
            cell = 'unable'
        one_patient_tc_cna.at[index,col] = cnv_color_dict[cell]

# II. Count the number of PB, RP, MLN, BM, NL, cfDNA samples
# ** Since most cohort has PB, RP, MLM, cfDNA, the count is esablished to 0.5 **
grouped_by_sample_type = one_patient_all_info.groupby(['Sample Type','Cohort','Patient ID']).size().reset_index(name='Counts')

pb_count = 0.5
rp_count = 0.5
mln_count = 0
mb_count = 0
nl_count = 0
cfdna_count = 0

for index, row in grouped_by_sample_type.iterrows():
    sample_type = grouped_by_sample_type.at[index,'Sample Type']
    if sample_type == 'PB':
        pb_count = grouped_by_sample_type.at[index,'Counts']
    elif sample_type == 'RP':
        rp_count = grouped_by_sample_type.at[index,'Counts']
    elif sample_type == 'MLN':
        mln_count = grouped_by_sample_type.at[index,'Counts']
    elif sample_type == 'MB':
        mb_count = grouped_by_sample_type.at[index,'Counts']
    elif sample_type == 'NL':
        nl_count = grouped_by_sample_type.at[index,'Counts']
    else:
        cfdna_count = grouped_by_sample_type.at[index,'Counts']

# II. Subdivide the one_patient tables into each sample type
# TODO: adding all sample types
pb_one_patient_tc_cna = one_patient_tc_cna[one_patient_tc_cna['Sample Type'].isin(['PB'])]
pb_one_patient_tc_cna = pb_one_patient_tc_cna.reset_index(drop=True)
pb_one_patient_mutations = one_patient_mutations[one_patient_mutations['Sample Type'].isin(['PB'])]
pb_one_patient_mutations = pb_one_patient_mutations.reset_index(drop=True)

rp_one_patient_tc_cna = one_patient_tc_cna[one_patient_tc_cna['Sample Type'].isin(['RP'])]
rp_one_patient_tc_cna = rp_one_patient_tc_cna.reset_index(drop=True)
rp_one_patient_mutations = one_patient_mutations[one_patient_mutations['Sample Type'].isin(['RP'])]
rp_one_patient_mutations = rp_one_patient_mutations.reset_index(drop=True)

mln_one_patient_tc_cna = one_patient_tc_cna[one_patient_tc_cna['Sample Type'].isin(['MLN'])]
mln_one_patient_tc_cna = mln_one_patient_tc_cna.reset_index(drop=True)
mln_one_patient_mutations = one_patient_mutations[one_patient_mutations['Sample Type'].isin(['MLN'])]
mln_one_patient_mutations = mln_one_patient_mutations.reset_index(drop=True)

mb_one_patient_tc_cna = one_patient_tc_cna[one_patient_tc_cna['Sample Type'].isin(['MB'])]
mb_one_patient_tc_cna = mb_one_patient_tc_cna.reset_index(drop=True)
mb_one_patient_mutations = one_patient_mutations[one_patient_mutations['Sample Type'].isin(['MB'])]
mb_one_patient_mutations = mb_one_patient_mutations.reset_index(drop=True)

cfdna_one_patient_tc_cna = one_patient_tc_cna[one_patient_tc_cna['Sample Type'].isin(['cfDNA'])]
cfdna_one_patient_tc_cna = cfdna_one_patient_tc_cna.reset_index(drop=True)
cfdna_one_patient_mutations = one_patient_mutations[one_patient_mutations['Sample Type'].isin(['cfDNA'])]
cfdna_one_patient_mutations = cfdna_one_patient_mutations.reset_index(drop=True)


"""
Ploting oncoprints
"""

# I. Establishing basic structure
fig = plt.figure(figsize = (width,height))
gs = gridspec.GridSpec(ncols = 6, nrows = 2, height_ratios = [2.5,15], width_ratios = [pb_count,rp_count,mln_count,mb_count,cfdna_count, 6.5])
gs.update(hspace = 0.07, wspace = 0, right = 0.97, top = 0.97, bottom = 0.25, left = 0.17)

axpb = fig.add_subplot(gs[0,0])
axrp = fig.add_subplot(gs[0,1])
axmln = fig.add_subplot(gs[0,2])
axmb = fig.add_subplot(gs[0,3])
axcfdna = fig.add_subplot(gs[0,4])

axpb1 = fig.add_subplot(gs[1,0],sharex=axpb)
axrp1 = fig.add_subplot(gs[1,1], sharex=axrp)
axmln1 = fig.add_subplot(gs[1,2],sharex=axmln)
axmb1 = fig.add_subplot(gs[1,3], sharex=axmb)
axcfdna1 = fig.add_subplot(gs[1,4], sharex=axcfdna)

# II. Ploting PB, RP, MLN, MB, cfDNA
ploting_oncoprint(pb_one_patient_tc_cna, pb_one_patient_mutations, first_plot=True, tc_ax=axpb, onco_axis=axpb1)
ploting_oncoprint(rp_one_patient_tc_cna, rp_one_patient_mutations, first_plot=False, tc_ax=axrp, onco_axis=axrp1)
ploting_oncoprint(mln_one_patient_tc_cna, mln_one_patient_mutations, first_plot=False, tc_ax=axmln, onco_axis=axmln1)
ploting_oncoprint(mb_one_patient_tc_cna, mb_one_patient_mutations, first_plot=False, tc_ax=axmb, onco_axis=axmb1)
ploting_oncoprint(cfdna_one_patient_tc_cna, cfdna_one_patient_mutations,first_plot=False, tc_ax=axcfdna, onco_axis=axcfdna1)

ax5 = fig.add_subplot(gs[1,5])
onco_labels = ['Amplification','Gain', 'Deletion', 'Deep deletion', 'Unable to call CNV', 'Non-frameshift indel','Missense', 'Truncation', 'Silent','>2 mutations', 'Non-independent call', 'Insufficient']
onco_handles = [Patch(color = '#EE2D24', linewidth = 0), Patch(color = '#F59496', linewidth = 0),
               Patch(color = '#9CC5E9', linewidth = 0), Patch(color = '#3F60AC', linewidth = 0),
               mlines.Line2D([], [], color='none', markeredgecolor='#E6E7E8', marker='s', lw=0, markersize=8),
               Patch(color = '#a9a9a9', linewidth = 0), Patch(color = '#79B443', linewidth = 0),
               Patch(color = '#FFC907', linewidth = 0), Patch(color = '#605857', linewidth = 0),
               mlines.Line2D([], [], color='black', markeredgecolor='k', marker='^', lw=0, markersize=8),
               mlines.Line2D([], [], color='black', markeredgecolor='k', marker='o', lw=0, markersize=8),
               mlines.Line2D([], [], color='none', markeredgecolor='k', marker='s', lw=0, markersize=8)]
ax5.legend(onco_handles, onco_labels, loc = 'center left', handlelength = 0.8, frameon = False,fontsize=legend_font_size)
ax5.spines['bottom'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)

ax6 = fig.add_subplot(gs[0,5])
labels = ['Gain', 'Deletion','Amplification', 'Deep deletion']
handles = [mlines.Line2D([], [], color='#F59496', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#9CC5E9', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#EE2D24', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#3F60AC', marker='_', lw=0, markersize=8)]
ax6.legend(handles, labels, loc = 'center left', title='Detection threshold',handlelength = 0.8, frameon = False,fontsize=legend_font_size)
ax6.spines['bottom'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)

filename = 'new_sep_' + cohort + '_' + patientID + '.pdf'
fig.savefig(filename,bbox_extra_artists=(ax5,), bbox_inches='tight')