# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 22:11:12 2020

@author: sng
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from datetime import datetime as dt
import natsort as ns


cohort = 'M1RP'
patientID = 'ID32'

all_target_genes = True    #Mark false if only want genes with mutations,not all target genes

#TODO: Filter out failed samples for M1B when more samples arrive
#TODO: Make option for font adjustment
#TODO: Could add bars over the top under ct_est plot to group samples


# =============================================================================
# Import data
# =============================================================================
#Get tumour content estimates
tc_est = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc_est.columns = tc_est.iloc[0]
tc_est = tc_est.drop(tc_est.index[0])

tc_est = tc_est[tc_est['Cohort'] == cohort]
tc_est = tc_est[tc_est['Patient ID'] == patientID]
#get rid of % sign on tc estimates and turn into float
if cohort != 'M1B':
    tc_est['Final tNGS_TC'] = list(map(lambda x: x[:-1], tc_est['Final tNGS_TC'].values))
    tc_est['Final tNGS_TC'] = [float(x) for x in tc_est['Final tNGS_TC'].values]

#Get somatic mutation data as dataframe with string values
mut = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/%s_mutations.tsv' % cohort, sep = '\t')
mut = mut[mut['Patient ID'] == patientID]
mut = mut.applymap(str)

#Get copy number data 
cn_cfdna = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_cfdna_cna.tsv' % cohort, sep = '\t')
cn_ffpe = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_ffpe_cna.tsv' % cohort, sep = '\t')

#Get times of sample collection 
timestamp = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/Clinical data tables/PB and RP dates.csv')

beta = pd.read_excel('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/%s_betastasis_all.xlsx' % cohort)
beta = beta[beta["EFFECT"] != "Intergenic"]
beta = beta.set_index('GENE', drop = False)
search = [patientID, 'GENE', 'EFFECT']
beta = beta.loc[:, beta.columns.str.contains('|'.join(search))]

# =============================================================================
# Filter copy number data
# =============================================================================
cn_cfdna = cn_cfdna[cn_cfdna['Patient ID'] == patientID]
cn_cfdna['Copy_num'] = cn_cfdna['Copy_num'].astype(str)

cn_ffpe = cn_ffpe[cn_ffpe['Patient ID'] == patientID]
cn_ffpe['Copy_num'] = cn_ffpe['Copy_num'].astype(str)

# =============================================================================
# Get list of genes to plot
# =============================================================================
# Get target genes 
targets = cn_cfdna['GENE'].unique().tolist()
cn_cfdna = cn_cfdna[cn_cfdna['Copy_num'] != '0']
cn_ffpe = cn_ffpe[cn_ffpe['Copy_num'] != '0']
# Exclude genes with no mutation or copy number changes
if all_target_genes == False:
    genes = mut['GENE'].tolist()
    genes = [x for x in genes if x != 'nan']
    genes = genes + cn_cfdna['GENE'].tolist()
    genes = genes + cn_ffpe['GENE'].tolist()
    genes = [x for x in targets if x in genes]
else: 
    genes = targets


# Filter out copy number changes from CN data that aren't in target panel
cn_cfdna = cn_cfdna[cn_cfdna['GENE'].isin(genes)]
cn_ffpe = cn_ffpe[cn_ffpe['GENE'].isin(genes)]
# =============================================================================
# Make a new dataframe to plot
# =============================================================================
mutations = tc_est['Sample ID']
mutations = pd.DataFrame(mutations)
mutations = mutations.rename(columns = {0: 'Sample ID'})
mutations = mutations.set_index('Sample ID', drop = False)

for gene in genes:
    mutations[gene] = ''

#Add mutation data to mutations dataframe
for row in mut.itertuples(index = False):
    sample = row[2]
    gene = row[7]
    mut_type = row[8]
    mut_type = mut_type.split(" ")[0]
    if ';' in gene:
        list_genes = gene.split(';')
        temp_gene = [x for x in list_genes if x in genes]
        if temp_gene:
           gene = temp_gene.pop() 
    if gene in genes:
        mutations.loc[sample, gene] = mutations.loc[sample, gene] + mut_type + ' '
    

#Add copy number data to mutations dataframe
for row in cn_cfdna.itertuples(index = False):
    sample = row[2]
    gene = row[3]
    copy_num = row[4]
    mutations.loc[sample, gene] = mutations.loc[sample, gene] + copy_num + ' '
for row in cn_ffpe.itertuples(index = False):
    sample = row[2]
    gene = row[3]
    copy_num = row[4]
    mutations.loc[sample, gene] = mutations.loc[sample, gene] + copy_num + ' '


#Reorder samples by time
timestamp = timestamp.loc[timestamp['Study ID '] == cohort + '_' + patientID]
pb_date = timestamp.iloc[0]['Date PB']
pb_date = dt.strptime(pb_date, '%d-%b-%Y')

if cohort == 'M1RP': #RP samples not here yet for M1B
    rp_date = timestamp.iloc[0]['Date RP']
    rp_date = dt.strptime(rp_date, '%d-%b-%Y')

mutations.insert(0, 'Date', '')
for i in mutations.index:
    if 'PB' in mutations.loc[i]['Sample ID']:
        mutations.at[i, 'Date'] = pb_date    
    elif 'MLN' in mutations.loc[i]['Sample ID'] or '_RP' in mutations.loc[i]['Sample ID']:
        mutations.at[i, 'Date'] = rp_date    
    elif 'cfDNA' in mutations.loc[i]['Sample ID']:
        date = mutations.loc[i]['Sample ID'][-9:]
        mutations.at[i, 'Date'] = dt.strptime(date, '%Y%b%d')
    elif 'MB' in mutations.loc[i]['Sample ID']:
        mutations.at[i, 'Date'] = dt.strptime('2018Oct', '%Y%b')
        
mutations = mutations.reset_index(drop = True)
#convert sample IDs to ordered categorical, naturally sorted 
mutations['Sample ID'] = pd.Categorical(mutations['Sample ID'], ordered=True, categories= ns.natsorted(mutations['Sample ID'].unique()))
mutations = mutations.sort_values(by = ['Date', 'Sample ID'])

# =============================================================================
# Set up oncoprint
# =============================================================================
samples = mutations['Sample ID']
mut_dot_size = 35                #Size of mutations squares - adjust if needed
fig_width = 0.3*len(samples)     #Width of plot - adjust constant if needed
fig_height = 0.3*len(genes)      #Height of plot - adjust constants if needed

fig = plt.figure(figsize = (fig_width, fig_height))
gs = gridspec.GridSpec(ncols = 2, nrows = 2, height_ratios = [1.5,20], width_ratios = [3,0.2])
gs.update(hspace = 0.02, wspace = 0, right = 0.97, top = 0.97, bottom = 0.25, left = 0.17)

if cohort == 'M1RP':
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex = ax1)
else:
    ax2 = fig.add_subplot(gs[1,0])
plt.subplots_adjust(bottom=.25, left=.25)

ax2.set_ylim(-0.1, len(genes))
ax2.set_yticks(np.arange(0.4, len(genes), 1))
ax2.set_yticklabels(genes)
ax2.tick_params(axis = 'both', length = 0)
ax2.invert_yaxis()

samples = mutations['Sample ID']
ax2.set_xlim(-0.6,len(samples))
ax2.set_xticks(np.arange(0, len(samples), 1))
ax2.set_xticklabels(samples, rotation = 90)
#ax1.xaxis.set_visible(False)


# =============================================================================
# Plot TC data
# =============================================================================
if cohort == 'M1RP':
    tc_est = tc_est.set_index('Sample ID', drop=False)
    tc_est = tc_est.reindex(samples)
    
    tc_est['Final tNGS_TC'] = tc_est['Final tNGS_TC'].fillna(0)
    ax1.bar(tc_est['Sample ID'], tc_est['Final tNGS_TC'], zorder = 100, color = '#909396')
    ax1.set_ylabel(r'$\bf{TC \%}$')
    #ax2.set_ylabel(r'$\bf{Somatic Mutations}$')
    ax1.set_ylim(0,100)
    
    #ax1.set_yscale('symlog', linthreshy=10)
    ax1.set_yticks([0, 20, 40, 60, 80, 100])
    ax1.set_yticklabels([0, 20, 40, 60, 80, 100])
    ax1.grid(axis = 'y', linestyle = 'dashed', linewidth = 0.5, color = '0.7', zorder = 0)
    
    ax1.spines['bottom'].set_zorder(100)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(bottom = False, labelbottom = False)

#horizontal line on barchart
#tc_ax.axhline(y=19.5,xmin=0,xmax=1,c="#3F60AC",linewidth=1,zorder=0)
    
    
# =============================================================================
# Plot mutation data
# =============================================================================
mutations = mutations.reset_index(drop = True)

#Unccoment line beside frameshift for splice, frameshift, and stopgain as one colour
def plot_scatter(row, string, x_offset, y_offset, zorder):
    if 'Missense' in string:
        ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#79B443', zorder = zorder, edgecolors = 'none', linewidth = 0.4)
    if 'Frameshift' in string: #or 'Stopgain' in string or 'Splice' in string:
        ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#FFC907', zorder = zorder, edgecolors = 'none', linewidth = 0.4)
    if 'Non-frameshift' in string:
        ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#8c69ff', zorder = zorder, edgecolors = 'none', linewidth = 0.4)
    if 'Stopgain' in string:
        ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#BD4398', zorder = zorder, edgecolors = 'none', linewidth = 0.4)
    if 'Splice' in string:
         ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = 'k', zorder = zorder, edgecolors = 'none', linewidth = 0.4)
    if 'Upstream' in string or 'Intronic' in string or 'Synonymous' in string or 'UTR' in string or 'Downstream' in string or 'Exonic' in string:
        ax2.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#757678', zorder = zorder, edgecolors = 'none', linewidth = 0.4)


#fill in background gray boxes
for index, row in mutations.iterrows():
    for i, col in enumerate(genes):
        ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#E6E7E8', zorder = 10)
        ax2.scatter(row.name, i+0.4, marker = 's', s = mut_dot_size, color = '#E6E7E8', zorder = 30, edgecolors = 'none', linewidth = 0)

                
                
#plotting mutations and copy number
for index, row in mutations.iterrows():
    for i, col in enumerate(genes):
        tc = tc_est.loc[row["Sample ID"]]['Final tNGS_TC']
        reads = beta.loc[beta['GENE'] == col, row["Sample ID"]].tolist()
        if reads:
            reads = [x[x.find("(")+1:x.find(")")] for x in reads]
            reads = float(min(reads))*tc/100
            if reads < 8:
                ax2.scatter(row.name, i+0.4, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 35, edgecolors = 'none', linewidth = 0)
        if mutations.at[index, col] != '':
            mut_count = sum(1 for c in mutations.at[index, col] if c.isupper())
            if 'UTR' in mutations.at[index, col]:
                mut_count = mut_count - (2 * mutations.at[index, col].count('UTR'))
            cn = mutations.at[index, col][-3:]
            m1 = mutations.at[index, col].split(' ')[0]
            m2 = mutations.at[index, col].split(' ')[1]
            if mut_count == 1:
                plot_scatter(row, m1, 0, 0.4, 100)
            elif mut_count == 2:        
                plot_scatter(row, m1, -0.09, 0.35, 100)
                plot_scatter(row, m2, 0.09, 0.49, 100)
            elif mut_count > 2:
                ax2.scatter(row['Sample ID'], i+0.37, marker = '^', s = 30, color = 'k', zorder = 100, edgecolors = 'k', linewidth = 0.8)
            if '-2' in cn:
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#3F60AC', zorder = 32)
            if '-1' in cn:
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#9CC5E9', zorder = 32)
            if '1' in cn and not ('-1' in cn):
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#F59496', zorder = 32)
            if '2' in cn and not ('-2' in cn):
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#EE2D24', zorder = 32)
        if tc == 0:     # white boxes for samples with TC too low to call CNV
            ax2.bar(row['Sample ID'], height=0.7, width=0.7, bottom = i+0.05, color = '#FFFFFF', zorder = 25)
            #For just blank, 
            #ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#FFFFFF', zorder = 25)



ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

# =============================================================================
# Legend
# =============================================================================
ax3 = fig.add_subplot(gs[1,1])
ax3.set_zorder(100)

labels = ['Missense','Frameshift indel','Stopgain','Non-frameshift indel', 'Splice', 'Silent mutation',
          'Amplification','Gain','Deletion','Deep deletion', '>2 mutations']
handles = [Patch(color = '#79B443', linewidth = 0),Patch(color = '#FFC907',linewidth = 0),
           Patch(color = '#BD4398', linewidth = 0),Patch(color = '#8c69ff',linewidth = 0),
           Patch(color = 'k',linewidth = 0), Patch(color = '#5c5c5c',linewidth = 0),
           Patch(color = '#EE2D24', linewidth = 0),Patch(color = '#F59496',linewidth = 0),
           Patch(color = '#9CC5E9', linewidth = 0),Patch(color = '#3F60AC',linewidth = 0),
           mlines.Line2D([], [], color='k', markeredgecolor='k', marker='^', lw=0, markersize=10, label='>2 mutations')]

ax3.legend(handles, labels, loc = 'center left', handlelength = 0.8, frameon = False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)

name = cohort + '_' + patientID + '_oncoprint'

#fig.savefig('C:\\Users\\Sarah\\Desktop\\testindonco.pdf', bbox_extra_artists=(ax3,), bbox_inches='tight')
fig.savefig('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/oncoprints/M1RP_marked_unable_orderedbytime/%s.pdf' %name, bbox_extra_artists=(ax3,), bbox_inches='tight')






