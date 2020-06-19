"""
Created on Wed Jun 18 20:00:28 2020

@author: echen
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import datetime

"""
Setting parameters
"""
# I. Formating
# 1. Font & Size
axis_label_font_size = 20
tick_label_font_size = 20
legend_font_size = 16

# 2. Axis names

# II. Artistics
# 1. cna
amplification_color = '#EE2D24'
gain_color ='#F59496'
deletion_color = '#9CC5E9'
deep_deletion_color = '#3F60AC'
no_change_color = '#E6E7E8'

# 2. Somatic mutations
missense_color = '#79B443'
truncation_color ='#FFC907'
silent_color = '#605857'
non_frameshift_indel_color = '#a9a9a9'

# 3. Figure size
width = 16
height = 28

# III. Data
cohort = 'M1RP'
patientID = 'ID10'
genes = []
gene_order = []
num_genes = 73
num_genes = num_genes + 1 + 1 # one for intergenic and another one for making it into column index

genes = pd.read_csv(r'C:\UBC\Coop\M1 project\baits_hg38.bed',names=['CHROM','START','END','GENE'],sep='\t')
tc = pd.read_csv(r"C:\UBC\Coop\M1 project\M1RP (Ghent) tables - Table SX - tNGS TC estimation.tsv", sep='\t')
mutations = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\mutations\final melted mutations\M1RP_mutations.tsv", sep='\t')
cfDNA_cna = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\copy number\final melted cna files\M1RP_cfdna_cna.tsv", sep='\t')
ffpe_cna = pd.read_csv(r"C:\Users\Emilia Chen\Dropbox\Ghent M1 2019\sandbox\copy number\final melted cna files\M1RP_FFPE_cna.tsv", sep='\t')
date = pd.read_excel(r"C:\UBC\Coop\M1 project\PB and RP dates.xlsx",index_col=0)

"""
Helper functions
"""
# Format the date into year-month-day as Date datatype for the ease of sorting
# Input = @date = in the format of YearMonthDay (e.g. 2020Jun19)
# Return = Year-Month-Day (e.g. 2020-06-19)
def translate_date(date):
    year = date[:4]
    mon_string = date[4:7].lower()
    day = date[7:]
    month_dict = {
        'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4, 'may': 5,'jun': 6, 
        'jul': 3, 'aug': 8, 'sep': 9, 'oct': 10, 'nov': 11, 'dec': 12
    }
    mon = month_dict[mon_string]
    return datetime.datetime(int(year), mon, int(day))

# Map somatic mutation to corresponding color
# Both input and output are of type 'str'
def translate_mutation_color(mutation):
    color = 'nan'
    mutation = mutation.lower()
    if 'missense' in mutation:
        color = missense_color
    elif 'non-frameshift indel' in mutation:
        color = non_frameshift_indel_color
    elif 'frameshift' in mutation:
        color = truncation_color
    elif 'stopgain' in mutation:
        color = truncation_color
    elif 'splice site' in mutation:
        color = truncation_color
    elif 'intronic' in mutation:
        color = silent_color
    elif 'intergenic'in mutation:
        color = silent_color
    elif 'upstream' in mutation:
        color = silent_color
    elif '3' in mutation:
        color = silent_color
    elif '5' in mutation:
        color = silent_color
    elif 'synonymous' in mutation:
        color = silent_color
    elif 'nan' == mutation:
        color = 'nan'
    else:
        color = 'nan'
    return color

# Map log_ratio to corresponding color
# @color is of type 'str'
def translate_cnv_color(lr):
    cell = float(lr)
    color = no_change_color
    if -1 < cell and cell <= -0.3:
        color = deletion_color
    elif  0.3 <= cell and cell < 1:
        color = gain_color
    elif cell >= 1:
        color = amplification_color
    elif cell <= -1:
        color = deep_deletion_color
    else:
        color = no_change_color
    return color   

"""
Reorganizing individual source table
"""
# I. Target genes 
# Add a column for chromosome number, if in X or Y chr, use 23 and 24 respectively, for the ease of ordering genes
genes['CHROM_NUM'] = None
for index,row in genes.iterrows():
    chrom = genes.loc[index,'CHROM']
    chrom_num = genes.loc[index,'CHROM'][3:]
    if chrom_num.lower() == 'x':
        chrom_num = '23'
    elif chrom_num.lower() == 'y':
        chrom_num = '24'
    genes.loc[index,'CHROM_NUM'] = chrom_num
genes['CHROM_NUM'] = genes['CHROM_NUM'].astype(int)
genes['START'] = genes['START'].astype(int)
# Sorting the target genes by CHROM_NUM and START position
genes = genes.sort_values(by=['CHROM_NUM','START'],ascending=[True,True])
genes = genes[['GENE']]
# Removing duplicates (as the bed file contains exon targets, so there are duplicated gene names)
genes = genes.drop_duplicates(['GENE'],keep='first')
genes = genes.reset_index(drop=True)
# Add 'intergenic'
intergenic = {'GENE':'intergenic'}
genes = genes.append(intergenic,ignore_index=True)
gene_order = genes['GENE'].tolist()

# II. Tumor content estimates
# Remove the first row due to the formatting of the source table
# Extract relevant columns and rows with cohort of interest
# Add a column to keep track of sample type
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[['Cohort','Patient ID', 'Sample ID', 'Final tNGS_TC']]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(str)
tc_M1RP = tc[tc['Cohort'].isin([cohort])]
tc_M1RP['Sample Type'] = None
for index, row in tc_M1RP.iterrows():
    sample = tc_M1RP.at[index, 'Sample ID']
    if 'cfDNA' in sample:
        tc_M1RP.at[index, 'Sample Type'] = 'cfDNA'
    elif 'MB' in sample:
        tc_M1RP.at[index, 'Sample Type'] = 'MB'
    elif 'PB' in sample:
        tc_M1RP.at[index, 'Sample Type'] = 'PB'
    elif 'MLN' in sample:
        tc_M1RP.at[index, 'Sample Type'] = 'MLN'
    elif 'NL' in sample:
        tc_M1RP.at[index, 'Sample Type'] = 'NL'
    else:
        tc_M1RP.at[index, 'Sample Type'] = 'RP'

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
# Extracting relevant columns and put 'intergenic' into GENE column
mutations = mutations[['Patient ID', 'Sample ID', 'GENE', 'EFFECT']]
mutations['EFFECT'] = mutations['EFFECT'].str.split(" ").str[0]
for index, row in mutations.iterrows():
    if mutations.at[index,'EFFECT'].lower() == 'intergenic':
        mutations.at[index, 'GENE'] = 'intergenic'
    if len(mutations.at[index,'GENE'].split(';')) > 1:
        for gene in mutations.at[index,'GENE'].split(';'):
            if gene in genes:
                mutations.at[index,'GENE'] = gene
# Sorting mutations so that entries with identical 'Sample ID', 'GENE', and 'EFFECT' are next to each other
# This is for latter when we concatenate entries with mutations at the same gene
mutations = mutations.sort_values(by=['Sample ID','GENE','EFFECT'])
mutations = mutations.reset_index(drop=True)
# Mark duplicates based on 'Sample ID' and 'GENE'
mutations['Duplicates'] = mutations.duplicated(subset=['Sample ID', 'GENE'], keep=False)
# Reorganize entries that have multiple mutations within one gene
# Since we sorted the mutation records, we know that the duplicates are next to one another
# Step 1: if the current row is in duplicate_row, then we drop the row
#         if the current row is not in duplicate_row, go to Step 2
# Step 2: a) translate mutation color
#         b) check to see if the current row has a duplicate
#         c) if there's a duplicate, then increment row number 
#                                     and append the duplicate effects to duplication_list 
#                                     and keep track of the row number
#            if there's no duplicate, the proceed to next row
num_row = len(mutations.index)
duplicate_row = list()
for index, row in mutations.iterrows():
    if index not in duplicate_row:
        effect_name = mutations.at[index, 'EFFECT'].lower()
        mutations.at[index, 'EFFECT'] = translate_mutation_color(effect_name)
        gene = mutations.at[index, 'GENE']
        next_row = index + 1
        duplication_list = ">1:" + translate_mutation_color(effect_name) 
        while next_row < num_row and mutations.at[next_row, 'Duplicates'] and mutations.at[next_row, 'GENE'] == gene:
            duplicate_row.append(next_row)
            another_effect = translate_mutation_color(mutations.at[next_row, 'EFFECT'])
            duplication_list = duplication_list +":" + another_effect
            next_row = next_row + 1
        if mutations.at[index, 'Duplicates']:
            mutations.at[index,'EFFECT'] = duplication_list
    else:
        mutations = mutations.drop(index)
        duplicate_row.remove(index)
mutations = mutations.drop('Duplicates',axis=1)

# V. cna
# Extracting relevant columns and append ffpe_cna to the back of cfDNA_cna
cfDNA_cna = cfDNA_cna[['Patient ID', 'Sample ID', 'GENE', 'Log_ratio']]
ffpe_cna = ffpe_cna[['Patient ID', 'Sample ID', 'GENE', 'Log_ratio']]
all_cna = cfDNA_cna.append(ffpe_cna, ignore_index = True)

"""
Merging date, tc, cna, mutation
"""
# I. Merging date to tumor content --> tc_M1RP_merge
# If cfDNA sample, date is from Sample ID
# TODO: MLN, MB, NL don't have date info
tc_M1RP_merge = tc_M1RP.merge(date_all, how='left', on=['Patient ID', 'Sample Type'])
tc_M1RP_merge = tc_M1RP_merge.drop('ID ',axis = 1)
for index,row in tc_M1RP_merge.iterrows():
    sample = tc_M1RP_merge.at[index,'Sample ID'].lower()
    if 'cfdna' in sample:
        date = sample.split("_")[3].lower()
        newdate = translate_date(date)
        tc_M1RP_merge.at[index, 'Date'] = newdate


# II. Produece a template with only all patient and sample
#     Merge genes to the template --> template_genes
template = tc_M1RP[['Cohort','Patient ID','Sample ID','Sample Type']]
genes['key'] = '0'
template['key'] = '0'
template_genes = template.merge(genes,how='left',on='key')
template_genes = template_genes.drop('key',axis=1)

# III. Merge template_genes with cna --> template_genes_cnv
template_genes_cnv = template_genes.merge(all_cna,how='left',on=['Patient ID','Sample ID','GENE'])

# IV. Merge template_genes_cnv with mutations --> template_genes_cnv_mutation
template_genes_cnv_mutation = template_genes_cnv.merge(mutations,how='left',on=['Patient ID','Sample ID','GENE'])

# V. Merge 'Log_ratio' and 'EFFECT' column and they are separated by ';' 
template_genes_cnv_mutation['EFFECT'] = template_genes_cnv_mutation['EFFECT'].astype(str)
for index, row in template_genes_cnv_mutation.iterrows():
    copy_number = template_genes_cnv_mutation.at[index,'Log_ratio']
    copy_num_color = translate_cnv_color(copy_number)
    effect = copy_num_color+ ";" + template_genes_cnv_mutation.at[index,'EFFECT']
    template_genes_cnv_mutation.at[index,'EFFECT'] = effect

# VI. Pivot the GENE column to the header row and rearrange the order of genes according to gene_order specified above
pivot_template_genes_cnv_mutation = template_genes_cnv_mutation.pivot(index='Sample ID', columns='GENE',values='EFFECT')
pivot_template_genes_cnv_mutation = pivot_template_genes_cnv_mutation.reset_index()
ordered = list()
ordered.append('Sample ID')
ordered.extend(gene_order)
pivot_template_genes_cnv_mutation = pivot_template_genes_cnv_mutation[ordered]

# VII. Merge the pivoted table from VI with tumor content and date --> pivot_template_genes_cnv_mutation_tc
pivot_template_genes_cnv_mutation_tc = pivot_template_genes_cnv_mutation.merge(tc_M1RP_merge, how='left',on=['Sample ID'])

"""
Get individual patient info
"""
# I. Take info from pivot_template_genes_cnv_mutation_tc
#    and divide into tc_cna table and mutations table
one_patient_all_info = pivot_template_genes_cnv_mutation_tc[pivot_template_genes_cnv_mutation_tc['Patient ID'].isin([patientID])]
one_patient_tc_cna = one_patient_all_info.copy()
one_patient_mutations = one_patient_all_info.copy()
for col in one_patient_tc_cna.columns.tolist()[1:num_genes]:
    one_patient_tc_cna[col] = one_patient_tc_cna[col].str.split(';').str[0]

for col in one_patient_mutations.columns.tolist()[1:num_genes]:
    one_patient_mutations[col] = one_patient_mutations[col].str.split(';').str[1]

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
gs.update(hspace = 0.07, wspace = 0.04, right = 0.97, top = 0.97, bottom = 0.25, left = 0.17)

axpb = fig.add_subplot(gs[0,0])
axrp = fig.add_subplot(gs[0,1])
axmln = fig.add_subplot(gs[0,2])
axmb = fig.add_subplot(gs[0,3])
axcfdna = fig.add_subplot(gs[0,4])

axpb1 = fig.add_subplot(gs[1,0])
axrp1 = fig.add_subplot(gs[1,1])
axmln1 = fig.add_subplot(gs[1,2])
axmb1 = fig.add_subplot(gs[1,3])
axcfdna1 = fig.add_subplot(gs[1,4])

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
    
    tc_ax.bar(tc_cna_table['Sample ID'], tc_cna_table['Final tNGS_TC'], zorder = 10, color = 'k')
    tc_ax.set_ylim(0,100)
    tc_ax.tick_params(bottom = False, labelbottom = False)
    tc_ax.spines['bottom'].set_zorder(100)
    tc_ax.spines['top'].set_visible(False)
    tc_ax.spines['right'].set_visible(False)
    tc_ax.grid(axis = 'y', linestyle = 'dashed', linewidth = 0.5, color = '0.7', zorder = 0)

    for index, row in tc_cna_table.iterrows():
        for i, col in enumerate(tc_cna_table.columns.tolist()[1:num_genes]):
            onco_axis.bar(row['Sample ID'], 0.8, bottom = i, color = row[col], zorder = 10)
    
    for index, row in mut_table.iterrows():
        for i, col in enumerate(mut_table.columns.tolist()[1:num_genes]):
            mut_table[col] = mut_table[col].astype(str)
            cell = mut_table.loc[index, col]
            jitter = 0
            more_than_two_mutations = 0
            if '>1' in cell:
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
                onco_axis.scatter(index+jitter, i+0.4+jitter, marker = 's', s = 100, color = cell1, zorder = 20, linewidth = 0.4)
                onco_axis.scatter(index-jitter, i+0.4-jitter, marker = 's', s = 100, color = cell2, zorder = 20, linewidth = 0.4)
            elif more_than_two_mutations != 0:
                onco_axis.scatter(index, i+0.4, marker = '^', s = 100, color = 'black', zorder = 40, linewidth = 0.4)
            else:
                onco_axis.scatter(index, i+0.4, marker = 's', s = 100, color = cell, zorder = 20, linewidth = 0.4)
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
    


# II. Ploting PB
ploting_oncoprint(pb_one_patient_tc_cna, pb_one_patient_mutations, first_plot=True, tc_ax=axpb, onco_axis=axpb1)
ploting_oncoprint(rp_one_patient_tc_cna, rp_one_patient_mutations, first_plot=False, tc_ax=axrp, onco_axis=axrp1)
ploting_oncoprint(mln_one_patient_tc_cna, mln_one_patient_mutations, first_plot=False, tc_ax=axmln, onco_axis=axmln1)
ploting_oncoprint(mb_one_patient_tc_cna, mb_one_patient_mutations, first_plot=False, tc_ax=axmb, onco_axis=axmb1)
ploting_oncoprint(cfdna_one_patient_tc_cna, cfdna_one_patient_mutations,first_plot=False, tc_ax=axcfdna, onco_axis=axcfdna1)

ax5 = fig.add_subplot(gs[1,5])
labels = ['Amplification','Gain', 'Deletion', 'Deep deletion','Non-frameshift indel','Missense', 'Truncation', 'Silent','>2 mutations']
handles = [Patch(color = '#EE2D24', linewidth = 0), Patch(color = '#F59496', linewidth = 0),
           Patch(color = '#9CC5E9', linewidth = 0), Patch(color = '#3F60AC', linewidth = 0),
           Patch(color = '#a9a9a9', linewidth = 0), Patch(color = '#79B443', linewidth = 0),
           Patch(color = '#FFC907', linewidth = 0), Patch(color = '#605857', linewidth = 0),
           mlines.Line2D([], [], color='black', markeredgecolor='k', marker='^', lw=0, markersize=13)]
ax5.legend(handles, labels, loc = 'center left', handlelength = 0.8, frameon = False,fontsize=legend_font_size)
ax5.spines['bottom'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)

filename = 'try_sep_' + cohort + '_' + patientID + '.pdf'
fig.savefig(filename)