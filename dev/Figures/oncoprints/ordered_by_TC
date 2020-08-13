# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 22:11:12 2020
@author: sng
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.legend as legend
from datetime import datetime as dt
import natsort as ns
import re

#TODO: Filter out failed samples for M1B when more samples arrive
#TODO: Could add bars over the top under ct_est plot to group samples

# =============================================================================
# Variables and constants
# =============================================================================
cohort = 'M1RP'
patientID = 'ID21'

all_target_genes = False   #Mark false if only want genes with mutations,not all target genes
order_by_time = False       #If False, order by sample TC 

#Figure height ratios
tc_h = 3.5
clinical_h = 1.5
drivers_h = 8
passengers_h = 20

hspace = 0.05       #Space between subplots

# Detection thresholds as tc% on tc estimate plots
del_threshold = 37.5
deep_del_threshold = 19.5
gain_threshold = 46
amp_threshold = 23.5

mpl.rc('hatch', color = '#FFFFFF', linewidth = 0.2)
mpl.rcParams['font.size'] = 9
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 9
mpl.rcParams['legend.handletextpad'] = '0.7'
mpl.rcParams['legend.labelspacing'] = '0.2'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['legend.handlelength'] = 0.4
plt.rcParams['legend.handleheight'] = 0.70


driver_genes = ['TP53', 'PTEN', 'RB1', 'SPOP', 'FOXA1', 'TMPRSS2', 'BRCA2', 'NKX3-1', 'CLU', 'NCOA2', 'MYC']


# =============================================================================
# Helpers
# =============================================================================
#Helper for plotting oncoprint heatmap
#Unccoment line beside frameshift for splice, frameshift, and stopgain as one colour
def plot_scatter(fig, row, string, x_offset, y_offset, zorder):
    if 'dep' not in string:
        hatch = ''
    else:
        hatch = '////////'
    if 'Missense' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#79B443', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)
    if 'Frameshift' in string: #or 'Stopgain' in string or 'Splice' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#FFC907', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)
    if 'Non-frameshift' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#8c69ff', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)
    if 'Stopgain' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#BD4398', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)
    if 'Splice' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = 'k', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)
    if 'Upstream' in string or 'Intronic' in string or 'Synonymous' in string or 'UTR' in string or 'Downstream' in string or 'Exonic' in string:
        fig.scatter(row.name+x_offset, i+y_offset, marker = 's', s = mut_dot_size, color = '#757678', zorder = zorder, edgecolors = 'none', linewidth = 0.4, hatch = hatch)           
     
# =============================================================================
# Import data
# =============================================================================
#Get tumour content estimates
tc_est = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
#Get somatic mutation data as dataframe with string values
mut = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/%s_mutations.tsv' % cohort, sep = '\t')
#Get copy number data 
cn_cfdna = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_cfdna_cna.tsv' % cohort, sep = '\t')
cn_ffpe = pd.read_csv('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_ffpe_cna.tsv' % cohort, sep = '\t')
#Get times of sample collection 
timestamp = pd.read_csv('C:/Users/Sarah/Dropbox/Clinical data tables/PB and RP dates.csv')
#Get raw betastasis table for current cohort
beta = pd.read_excel('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/%s_betastasis_all.xlsx' % cohort)
#Get pathological data: Gleason grade, path TC, 
path_data = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

# =============================================================================
# Format dataframes - filtering down to current patient, changing data types
# =============================================================================
tc_est.columns = tc_est.iloc[0]
tc_est = tc_est.drop(tc_est.index[0])
tc_est = tc_est[tc_est['Cohort'] == cohort]             # Filter tc estimates for this patient only
tc_est = tc_est[tc_est['Patient ID'] == patientID]

#get rid of % sign on tc estimates and turn into float, taking empty values as 0
tc_est['Final tNGS_TC'] = tc_est['Final tNGS_TC'].fillna('0%')
tc_est['Final tNGS_TC'] = list(map(lambda x: x[:-1], tc_est['Final tNGS_TC'].values))
tc_est['Final tNGS_TC'] = [float(x) for x in tc_est['Final tNGS_TC'].values]

mut = mut[mut['Patient ID'] == patientID]
mut = mut.applymap(str)

beta = beta[beta["EFFECT"] != "Intergenic"]
search = [patientID + '_', 'GENE', 'EFFECT']
beta = beta.loc[:, beta.columns.str.contains('|'.join(search))]

path_data.columns = path_data.iloc[0]
path_data = path_data.drop(path_data.index[0])
path_data = path_data.set_index('Sample ID', drop=False)

path_data = path_data[path_data['Cohort'] == cohort]
path_data = path_data[path_data['Patient ID'] == patientID]

cn_cfdna = cn_cfdna[cn_cfdna['Patient ID'] == patientID]
cn_cfdna['Copy_num'] = cn_cfdna['Copy_num'].astype(str)

cn_ffpe = cn_ffpe[cn_ffpe['Patient ID'] == patientID]
cn_ffpe['Copy_num'] = cn_ffpe['Copy_num'].astype(str)

# =============================================================================
# Get list of genes to plot for y axis labels
# =============================================================================
# Get all target genes 
genes = cn_cfdna['GENE'].unique().tolist()

# Remove genes with no copy number changes from CN dataframes
cn_cfdna = cn_cfdna[cn_cfdna['Copy_num'] != '0']
cn_ffpe = cn_ffpe[cn_ffpe['Copy_num'] != '0']

#Change gene name in betastasis datafram for genes in format "gene1;gene2"
for index, row in beta.iterrows():
    gene = row['GENE']
    if ';' in gene:
        if gene.split(';')[0] in genes:
            beta.loc[index, 'GENE'] = gene.split(';')[0]
        if gene.split(';')[1] in genes:
            beta.loc[index, 'GENE'] = gene.split(';')[1]

beta = beta.set_index('GENE', drop = False)

# =============================================================================
# Make a new dataframes to organize mutations and CN data per gene and sample
# =============================================================================
#Make new dataframe with Sample Id as index
mutations = tc_est['Sample ID']
mutations = pd.DataFrame(mutations)
mutations = mutations.rename(columns = {0: 'Sample ID'})
mutations = mutations.set_index('Sample ID', drop = False)
#Add pathological info
mutations = mutations.join(path_data['GGG'])

for gene in genes:
    mutations[gene] = ''

# Make list of sample IDs and change tc_est index to Sample ID
samples = mutations['Sample ID']
tc_est = tc_est.set_index('Sample ID', drop=False)
tc_est = tc_est.reindex(samples)


#Add mutation data to mutations dataframe in format 'mut_mutation1 mut_mutation2'
for row in mut.itertuples(index = False):
    sample = row[2]
    if tc_est.loc[sample]['Final tNGS_TC'] > 0:
        gene = row[7]
        mut_type = row[8]
        mut_type = mut_type.split(" ")[0]
        if ';' in gene:
            list_genes = gene.split(';')
            temp_gene = [x for x in list_genes if x in genes]
            if temp_gene:
                gene = temp_gene.pop() 
        if gene in genes:
            mutations.loc[sample, gene] = mutations.loc[sample, gene] + 'mut_' + mut_type + '|'
    

#Add insufficient info to call (uncallable) and dependently called muts to mutations dataframe
for sample in samples:
    tc = tc_est.loc[sample]['Final tNGS_TC']
    for i, col in enumerate(genes):        
        read_counts = beta.loc[beta['GENE'] == col]
        for index2, row2 in read_counts.iterrows():
            has_muts = False
            for x in row2[2:]:
                if "*" in x:
                    has_muts = True
            if has_muts and not read_counts.empty and '*' not in row2[sample]:
                total_reads = float(re.search(r'\((.*?)\)', row2[sample]).group(1))
                vaf = float(row2[sample].split('%')[0])/100
                if vaf >= 0.01 and vaf*total_reads >= 3:
                    mutations.loc[sample, col] = mutations.loc[sample, col] + 'dep_' + row2['EFFECT'] + '|'
                elif total_reads*tc/100 < 8:    
                    mutations.loc[sample, col] = mutations.loc[sample, col] + 'uncallable '
                    
#Add copy number data to mutations dataframe in format 'mutation1|CNA_-1'
for row in cn_cfdna.itertuples(index = False):
    sample = row[2]
    if tc_est.loc[sample]['Final tNGS_TC'] > 0:
        gene = row[3]
        copy_num = row[4]
        mutations.loc[sample, gene] = mutations.loc[sample, gene] + 'cna_' + copy_num + '|'
for row in cn_ffpe.itertuples(index = False):
    sample = row[2]
    if tc_est.loc[sample]['Final tNGS_TC'] > 0:
        gene = row[3]
        copy_num = row[4]
        mutations.loc[sample, gene] = mutations.loc[sample, gene] + 'cna_' + copy_num + '|'


# =============================================================================
# Order samples by time or TC content
# =============================================================================
if order_by_time == True:
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
    samples = mutations['Sample ID']
else:
    mutations = mutations.join(tc_est['Final tNGS_TC'])
    mutations = mutations.sort_values(by = ['Final tNGS_TC'], ascending = False)
    mutations = mutations.reset_index(drop = True)
    samples = mutations['Sample ID']
    
    
#Separate driver genes into new dataframe
new_columns = driver_genes.copy()
new_columns.insert(0, 'Sample ID')
driver_mutations = mutations[new_columns].copy()

for gene in driver_genes:
    del mutations[gene]    
    
# If all_target_genes = False, delete genes without CNA or mutations
if all_target_genes == False:
    mutations = mutations.loc[:, (mutations!='').any(axis=0)]
    genes = mutations.columns.tolist()
    genes.remove('GGG')
    genes.remove('Final tNGS_TC')
    genes.remove('Sample ID')
    
# =============================================================================
# Set up oncoprint
# =============================================================================
mut_dot_size = 35                #Size of mutations squares - adjust if needed
fig_width = 0.3*len(samples)     #Width of plot - adjust constant if needed
fig_height = 0.2*len(genes) + 5      #Height of plot - adjust constants if needed

#Set figure dimensions and ratios
fig = plt.figure(figsize = (fig_width, fig_height))
gs = gridspec.GridSpec(ncols = 2, nrows = 4, height_ratios = [tc_h, clinical_h, drivers_h, passengers_h], width_ratios = [3,0.2])
gs.update(hspace = hspace, wspace = 0, right = 0.97, top = 0.97, bottom = 0.25, left = 0.17)

#Add and space out subplots where ax1 is tc plot and ax2 is oncoprint
ax1 = fig.add_subplot(gs[0,0])
clinical = fig.add_subplot(gs[1,0], sharex = ax1)
drivers = fig.add_subplot(gs[2,0], sharex = ax1)
ax2 = fig.add_subplot(gs[3,0], sharex = ax1)
plt.subplots_adjust(bottom=.25, left=.25)

#Set y and x axes for passengers oncoprint
ax2.set_ylim(-0.1, len(genes))
ax2.set_yticks(np.arange(0.4, len(genes), 1))
ax2.set_yticklabels(genes)
ax2.tick_params(axis = 'both', length = 0)
ax2.invert_yaxis()
ax2.set_xlim(-0.6,len(samples))
ax2.set_xticks(np.arange(0, len(samples), 1))
ax2.set_xticklabels(samples, rotation = 90)
#ax1.xaxis.set_visible(False)

#Set y and x axes for drivers oncoprint
drivers.set_ylim(-0.1, len(driver_genes))
drivers.set_yticks(np.arange(0.4, len(driver_genes), 1))
drivers.set_yticklabels(driver_genes)
drivers.tick_params(axis = 'both', length = 0)
drivers.tick_params(bottom = False, labelbottom = False)
drivers.invert_yaxis()

#Set y and x axes for clinical data matrix
clinical.set_ylim(-0.1, 2)
clinical.set_yticks(np.arange(0.4, 2, 1))
clinical.set_yticklabels(['Sample Type', 'GGG'])
clinical.tick_params(axis = 'both', length = 0)
clinical.tick_params(bottom = False, labelbottom = False)
clinical.invert_yaxis()
# =============================================================================
# Plot Clinical Matrix
# =============================================================================
mutations['GGG'] = mutations['GGG'].fillna('0')

for index, row in mutations.iterrows():
        #Sample type matrix row
        sample = row['Sample ID']
        if '_MLN' in sample:
            color = '#aae6bb'
        elif '_cfDNA' in sample:
            color = '#994545'
        elif '_RP' in sample:
            color = '#f7d58b'
        elif '_PB' in sample:
            color = '#667bb3'  
        elif '_MB' in sample:
            color = '#dea2e8'
        clinical.bar(row['Sample ID'], 0.8, bottom = 0, color = color, zorder = 10)
        #GGG matrix row
        gleason = row['GGG']
        if gleason == '5':
            color = '#0a004d'
        elif gleason == '4':
            color = '#422fc2'
        elif gleason == '3':
            color = '#7b68fc'
        elif gleason == '2':
            color = '#a699ff'
        elif gleason == '1':
            color = '#dcd6ff'            
        elif gleason == '0':
            color = '#ededed'
        clinical.bar(row['Sample ID'], 0.8, bottom = 1, color = color, zorder = 15)
    
# =============================================================================
# Plot TC  barchart
# =============================================================================
# Reset index of tc_est to Sample ID
tc_est = tc_est.set_index('Sample ID', drop=False)
tc_est = tc_est.reindex(samples)

#Plot bar chart
ax1.bar(tc_est['Sample ID'], tc_est['Final tNGS_TC'], zorder = 100, color = '#585b5e')

#Set up axes        
ax1.set_ylabel(r'$\bf{TC \%}$')
ax1.set_ylim(0,100)
minor_ticks = np.arange(0, 101, 10)
ax1.set_yticks([0, 20, 40, 60, 80, 100])
ax1.set_yticks(minor_ticks, minor = True)
ax1.set_yticklabels([0, 20, 40, 60, 80, 100])
ax1.grid(which = 'both', axis = 'y', linestyle = 'dashed', linewidth = 0.3, color = '0.7', alpha = 0.2, zorder = 0)
    
ax1.spines['bottom'].set_zorder(100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.tick_params(bottom = False, labelbottom = False)

#Setting detection threshold lines on tc plot
ax1.axhline(y=del_threshold, xmin=0,xmax=1, c="#9CC5E9", linewidth=0.7, zorder=20)
ax1.axhline(y=deep_del_threshold, xmin=0, xmax=1, c="#3F60AC", linewidth=0.7, zorder=20)
ax1.axhline(y=gain_threshold, xmin=0, xmax=1, c="#F59496", linewidth=0.7, zorder=20)
ax1.axhline(y=amp_threshold, xmin=0, xmax=1, c="#EE2D24", linewidth=0.7, zorder=20)
ax1.axhline(y=50, xmin=0, xmax=1, linestyle = '--', c="#383838", linewidth=0.7, zorder=20)

# =============================================================================
# Plot mutation data onto oncoprint
# =============================================================================
#TODO: Refactor out a helper function, right now code blocks below are highly redundant
            
            
mutations = mutations.reset_index(drop = True)         

#Plotting top oncoprint for driver mutations
for index, row in driver_mutations.iterrows():
    for i, col in enumerate(driver_genes):
        #Plot background gray boxes
        drivers.bar(row['Sample ID'], 0.8, bottom = i, color = '#E6E7E8', zorder = 10)
        
        #Get mutation and cn info from mutations table
        tc = tc_est.loc[row["Sample ID"]]['Final tNGS_TC']
        
        if driver_mutations.at[index, col] != '':
            cell = driver_mutations.at[index, col]
            mut_count = cell.count('mut')
            dep_count = cell.count('dep')
            uncallable_count = cell.count('uncallable')
            cn = ''
            if 'cna_' in cell:
                cn = cell[-5:]
            
            #Plot called and dependent mutations
            if mut_count + dep_count == 1 and uncallable_count == 0:
                m1 = cell.split('|')[0]
                plot_scatter(drivers, row, m1, 0, 0.4, 100)
            elif mut_count + dep_count == 2:
                m1 = cell.split('|')[0]
                m2 = cell.split('|')[1]
                plot_scatter(drivers, row, m1, -0.09, 0.35, 100)
                plot_scatter(drivers, row, m2, 0.09, 0.49, 100)
            elif mut_count + dep_count > 2:             
                drivers.scatter(row['Sample ID'], i+0.37, marker = '^', s = 30, color = 'k', zorder = 100, edgecolors = 'k', linewidth = 0.8)
            
            #Plot uncallable white boxes
            if uncallable_count == 1 and mut_count + dep_count == 0:  
                drivers.scatter(row.name, i+0.4, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, edgecolors = 'none', linewidth = 0.4)
            if uncallable_count >= 1 and mut_count + dep_count == 1:         
                m1 = cell.split('|')[0]
                plot_scatter(drivers, row, m1, 0.09, 0.49, 100)
                drivers.scatter(row.name-0.09, i+0.35, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, edgecolors = 'none', linewidth = 0)
            if uncallable_count > 1 and mut_count == 0:
                drivers.scatter(row.name-0.09, i+0.35, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, linewidth = 0)
                drivers.scatter(row.name+0.09, i+0.49, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, linewidth = 0)
            
            #Plot copy number changes                    
            if '-2' in cn:
                drivers.bar(row['Sample ID'], 0.8, bottom = i, color = '#3F60AC', zorder = 32)
            if '-1' in cn:
                drivers.bar(row['Sample ID'], 0.8, bottom = i, color = '#9CC5E9', zorder = 32)
            if '1' in cn and not ('-1' in cn):
                drivers.bar(row['Sample ID'], 0.8, bottom = i, color = '#F59496', zorder = 32)
            if '2' in cn and not ('-2' in cn):
                drivers.bar(row['Sample ID'], 0.8, bottom = i, color = '#EE2D24', zorder = 32)
       
        if tc == 0:     # white boxes for samples with TC too low to call CNV
            drivers.bar(row['Sample ID'], height=0.7, width=0.7, bottom = i+0.05, color = '#FFFFFF', zorder = 120)


                
#plotting mutations and copy number based on mutations dataframe
for index, row in mutations.iterrows():
    for i, col in enumerate(genes):
        #Plot background gray boxes
        ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#E6E7E8', zorder = 10)
        
        #Get mutation and cn info from mutations table
        tc = tc_est.loc[row["Sample ID"]]['Final tNGS_TC']
        
        if mutations.at[index, col] != '':
            cell = mutations.at[index, col]
            mut_count = cell.count('mut')
            dep_count = cell.count('dep')
            uncallable_count = cell.count('uncallable')
            cn = ''
            if 'cna_' in cell:
                cn = cell[-5:]
            
            #Plot called and dependent mutations
            if mut_count + dep_count == 1 and uncallable_count == 0:
                m1 = cell.split('|')[0]
                plot_scatter(ax2, row, m1, 0, 0.4, 100)
            elif mut_count + dep_count == 2:
                m1 = cell.split('|')[0]
                m2 = cell.split('|')[1]
                plot_scatter(ax2, row, m1, -0.09, 0.35, 100)
                plot_scatter(ax2, row, m2, 0.09, 0.49, 100)
            elif mut_count + dep_count > 2:             
                ax2.scatter(row['Sample ID'], i+0.37, marker = '^', s = 30, color = 'k', zorder = 100, edgecolors = 'k', linewidth = 0.8)
            
            #Plot uncallable white boxes
            if uncallable_count == 1 and mut_count + dep_count == 0:  
                ax2.scatter(row.name, i+0.4, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, edgecolors = 'none', linewidth = 0.4)
            if uncallable_count >= 1 and mut_count + dep_count == 1:         
                m1 = cell.split('|')[0]
                plot_scatter(ax2, row, m1, 0.09, 0.49, 100)
                ax2.scatter(row.name-0.09, i+0.35, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, edgecolors = 'none', linewidth = 0)
            if uncallable_count > 1 and mut_count == 0:
                ax2.scatter(row.name-0.09, i+0.35, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, linewidth = 0)
                ax2.scatter(row.name+0.09, i+0.49, marker = 's', s = mut_dot_size, color = '#FFFFFF', zorder = 33, linewidth = 0)
            
            #Plot copy number changes                    
            if '-2' in cn:
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#3F60AC', zorder = 32)
            if '-1' in cn:
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#9CC5E9', zorder = 32)
            if '1' in cn and not ('-1' in cn):
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#F59496', zorder = 32)
            if '2' in cn and not ('-2' in cn):
                ax2.bar(row['Sample ID'], 0.8, bottom = i, color = '#EE2D24', zorder = 32)
       
        if tc == 0:     # white boxes for samples with TC too low to call CNV
            ax2.bar(row['Sample ID'], height=0.7, width=0.7, bottom = i+0.05, color = '#FFFFFF', zorder = 120)

for i, f in enumerate([ax2, drivers, clinical]):
    f.spines['right'].set_visible(False)
    f.spines['left'].set_visible(False)
    f.spines['bottom'].set_visible(False)
    f.spines['top'].set_visible(False)    

# =============================================================================
# Legend
# =============================================================================

#Legend for clinical matrix
clinical_legend = fig.add_subplot(gs[1,1])
labels1 = ['N/A', '1', '2', '3', '4', '5']
handles1 = [Patch(color = '#ededed', linewidth = 0), 
            Patch(color = '#dcd6ff', linewidth = 0),Patch(color = '#a699ff', linewidth = 0),
            Patch(color='#7b68fc', linewidth = 0),Patch(color = '#422fc2', linewidth = 0),
            Patch(color = '#0a004d', linewidth = 0)]
legend1 = legend.Legend(clinical_legend,handles=handles1, labels=labels1,ncol=7, loc = (0,0),
                        columnspacing=None, handlelength = 0.8, frameon = False,fontsize=8)

labels2 = ['PB', 'RP', 'MLN', 'cfDNA', 'MB']
handles2 = [Patch(color = '#667bb3', linewidth = 0), 
            Patch(color = '#f7d58b', linewidth = 0),Patch(color = '#aae6bb', linewidth = 0),
            Patch(color='#994545', linewidth = 0),  Patch(color='#dea2e8', linewidth = 0)]
legend2 = legend.Legend(clinical_legend,handles=handles2, labels=labels2,ncol=6, loc = (0,0.5), 
                        columnspacing=None, handlelength = 0.8, frameon = False,fontsize=8)

clinical_legend.add_artist(legend2)
clinical_legend.add_artist(legend1)
clinical_legend.spines['bottom'].set_visible(False)
clinical_legend.spines['left'].set_visible(False)
clinical_legend.spines['right'].set_visible(False)
clinical_legend.spines['top'].set_visible(False)
clinical_legend.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)


ax3 = fig.add_subplot(gs[3,1])
ax3.set_zorder(100)

# Oncoprint legend elements and corresponding symbols
labels = ['Missense', 'Frameshift indel', 'Stopgain', 'Non-frameshift indel', 'Splice', 'Silent mutation', 'Unable to call', 'Dependent call',
          'Amplification', 'Gain', 'Deletion', 'Deep deletion', '>2 mutations']
handles = [Patch(color = '#79B443', linewidth = 0),Patch(color = '#FFC907',linewidth = 0),
           Patch(color = '#BD4398', linewidth = 0),Patch(color = '#8c69ff',linewidth = 0),
           Patch(color = 'k',linewidth = 0), Patch(color = '#5c5c5c',linewidth = 0),
           mlines.Line2D([], [], color='none', markeredgecolor='k', marker='s', lw=0, markersize=8),
           Patch(edgecolor = '#FFFFFF', facecolor='k', linewidth = 0, hatch = '////////'),
           Patch(color = '#EE2D24', linewidth = 0),Patch(color = '#F59496',linewidth = 0),
           Patch(color = '#9CC5E9', linewidth = 0),Patch(color = '#3F60AC',linewidth = 0),
           mlines.Line2D([], [], color='k', markeredgecolor='k', marker='^', lw=0, markersize=10, label='>2 mutations')]

    
# Formatting
ax3.legend(handles, labels, loc = 'upper left', handlelength = 0.8, frameon = False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)

# TC plot legend elements and corresponding symbols
ax4 = fig.add_subplot(gs[0,1])
labels = ['Gain', 'Amplification', 'Deletion', 'Deep deletion', 'Distinguish deletion from deep deletion']
handles = [mlines.Line2D([], [], color='#F59496', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#EE2D24', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#9CC5E9', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#3F60AC', marker='_', lw=0, markersize=8),
           mlines.Line2D([], [], color='#383838', linestyle='dashed', lw=0.8, markersize=8)]
ax4.legend(handles, labels, loc = 'center left', title='Detection thresholds', title_fontsize =9, handlelength = 0.9, frameon = False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False)



# =============================================================================
# Save figure 
# =============================================================================
name = cohort + '_' + patientID + '_oncoprint'

fig.savefig('C:\\Users\\Sarah\\Desktop\\%s.pdf' %name, bbox_extra_artists=(clinical_legend, ax4), bbox_inches='tight')
#fig.savefig('C:/Users/Sarah/Dropbox/Ghent M1 2019/sandbox/oncoprints/M1B_orderedbyTC/%s.pdf' %name, bbox_extra_artists=(ax4,), bbox_inches='tight')

