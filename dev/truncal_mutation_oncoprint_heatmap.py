# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:31:59 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from colour import Color

# =============================================================================
# Constants
# =============================================================================

min_shared_threshold = 0.3
min_truncal_threshold = 0.8
cohort = 'M1RP'

# ============================================================================
# Helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False))]

def get_color(per, colorslist):
    return colorslist[int(per*100-6)];

def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])

def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)

# =============================================================================
# Import data
# =============================================================================
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/%s_mutations.tsv' % cohort, sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[tc['Cohort'] == cohort]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
if cohort == 'M1RP':
    tc = tc[tc['Patient ID'] != 'ID8']

# =============================================================================
# Keep tumor content above 15%
# =============================================================================

tc = tc[tc['Final tNGS_TC'] >= (2/(1+1/0.08))]
muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]

# =============================================================================
# Keep coding mutations
# =============================================================================

muts = keepCodingMutations(muts)

# =============================================================================
# Get patient sample counts
# =============================================================================

pt_sample_counts = tc.groupby('Patient ID').count()[['Sample ID']]
pt_sample_counts = pt_sample_counts[pt_sample_counts['Sample ID'] > 1]
if cohort == 'M1RP':
    pt_sample_counts = pt_sample_counts[~pt_sample_counts.index.isin(['ID25','ID20'])]
pt_sample_counts.columns = ['Sample count']
pt_sample_counts = pt_sample_counts.reset_index()

# =============================================================================
# Get mutation counts
# =============================================================================

muts = muts[muts['Patient ID'].isin(pt_sample_counts['Patient ID'])]
mut_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
mut_counts.columns = ['Mutant samples']
mut_counts = mut_counts.reset_index()

# =============================================================================
# Merge sample counts onto mut tests
# =============================================================================

mut_counts = mut_counts.merge(pt_sample_counts, on = 'Patient ID')
mut_counts['Percent'] = (mut_counts['Mutant samples'] / mut_counts['Sample count']).round(2)
tmp = mut_counts.groupby(['Patient ID','GENE']).count()[['POSITION']].rename(columns = {'POSITION':'Gene_mut_count'}).reset_index()
mut_counts = mut_counts.merge(tmp, on = ['Patient ID','GENE'])

# =============================================================================
# Set up oncoprint
# =============================================================================

genes = muts['GENE'].unique().tolist()
pts = pt_sample_counts['Patient ID'].unique().tolist()

matrix = pd.DataFrame(index = genes, columns = pts)
matrix2 = pd.DataFrame(index = genes, columns = pts)


# =============================================================================
# Get color list
# =============================================================================

colorslist = linear_gradient('#FFC0CB', '#b20000', 95)['hex']

# =============================================================================
# Fill matrix
# =============================================================================

null_color = '#efefef'


for index, row in mut_counts.iterrows():
    if not pd.isnull(matrix.at[row['GENE'],row['Patient ID']]):
        matrix2.at[row['GENE'],row['Patient ID']] = get_color(row['Percent'], colorslist)
    matrix.at[row['GENE'],row['Patient ID']] = get_color(row['Percent'], colorslist)
    
matrix = matrix.fillna(null_color)

# =============================================================================
# Create reorder matrix
# =============================================================================

order_matrix = matrix.copy()
order_matrix = order_matrix.replace({null_color:np.nan})
gene_order = order_matrix.count(axis = 1)
matrix['order'] = gene_order
matrix = matrix.sort_values('order',ascending = False).drop('order', axis = 1)

matrix2['order'] = gene_order
matrix2 = matrix2.sort_values('order',ascending = False).drop('order', axis = 1)

order_matrix = matrix.copy()
for index, row in mut_counts.iterrows():
    order_matrix.at[row['GENE'],row['Patient ID']] = row['Percent']
    
order_matrix = order_matrix.replace({'#efefef':0})
for i, (index, row) in enumerate(order_matrix.iterrows()):
    for col in order_matrix.columns:
        order_matrix.at[index,col] = order_matrix.at[index,col]*10**(len(order_matrix) - i)
pt_order = order_matrix.sum()
pt_order = pt_order.sort_values(ascending = False)
matrix = matrix[pt_order.index.tolist()]
matrix2 = matrix2[pt_order.index.tolist()]

# =============================================================================
# Plot matrix
# =============================================================================

fig, [ax,ax2] = plt.subplots(ncols = 2, gridspec_kw={'width_ratios':[1,0.1]}, figsize = (8,5))

for y, (index, row) in enumerate(matrix.iterrows()):
    for x, col in enumerate(matrix.columns):
        if pd.isnull(matrix2.at[index, col]):
            color1 = matrix.at[index, col]
            ax.bar(x, 0.8, bottom = y, color = color1, align = 'edge', linewidth = 0)
        else:
            color1 = matrix.at[index, col]
            color2 = matrix2.at[index, col]
            ax.fill([x,x+0.8,x+0.8],[y,y,y+0.8], facecolor = color1, linewidth = 0, edgecolor = '#efefef', zorder = 100)
            ax.fill([x,x,x+0.8],[y,y+0.8,y+0.8],facecolor = color2, linewidth = 0, edgecolor = '#efefef', zorder = 100)
            ax.plot([x,x+0.8],[y,y+0.8], c = 'w', zorder = 1000, lw = 0.5)

ax.set_yticks(np.arange(0.4,len(genes),1))
ax.set_yticklabels(gene_order.sort_values(ascending = False).index)
ax.set_xticks(np.arange(0.4,len(pts),1))
ax.set_xticklabels(pt_order.index, rotation = 90)
ax.set_title(cohort)

ax.set_xlim(-0.1, len(pts))
ax.set_ylim(len(genes),-0.1)

ax2.yaxis.set_visible(False)
ax2.xaxis.set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

handles = [mpatches.Patch(color = '#efefef'),
           mpatches.Patch(color = colorslist[0]),
           mpatches.Patch(color = colorslist[4]),
           mpatches.Patch(color = colorslist[24]),
           mpatches.Patch(color = colorslist[34]),
           mpatches.Patch(color = colorslist[44]),
           mpatches.Patch(color = colorslist[54]),
           mpatches.Patch(color = colorslist[64]),
           mpatches.Patch(color = colorslist[74]),
           mpatches.Patch(color = colorslist[84]),
           mpatches.Patch(color = colorslist[94]),]
labels = ['Not present',
          '5% of samples',
          '20%',
          '30%',
          '40%',
          '50%',
          '60%',
          '70%',
          '80%',
          '90%',
          '100%',]

ax2.legend(handles, labels, fontsize = 10, handlelength = 0.8, labelspacing = 0)
plt.tight_layout()

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/%struncal_mutation_oncoprint.pdf' % cohort)
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/%struncal_mutation_oncoprint.png' % cohort, dpi = 400)