# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:47:52 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

pt = 'ID10'

# =============================================================================
# Helpers
# =============================================================================

    # =============================================================================
    # make cytoband file into df with arm starts and ends
    # inputs = path to cytoband file and segfile data frame
    # outputs = chromosome arm list, chromosome list and chr_arm defining arm boundaries
    # =============================================================================
'''turn cytoband file from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
to a dataframe with columns ['chr', 'start', 'end', 'arm', 'chr_arm'].
p starts at 0, goes to start of centromere. q starts at q centromere and goes to max position in
the seg file position'''
def make_centro(centro):
    centro.columns = ["chr", "start", "end", "band", "annot"] #centromere position file
    #get rid of stupid chr in chrom names
    centro["chr"] = centro["chr"].replace("chr", "", regex=True)
    #acen is a stain that ids centromeres.
    #This removes all the cyto bands that arent centromeres
    maxseg = centro.groupby("chr")["end"].max().reset_index()
    centro = centro.query('annot == "acen"')
    centro["arm"] = centro["band"].replace("[0-9.]", "", regex=True)
    #p arm goes 0 to centromere
    #q goes centromere to max seg position
    centro.loc[centro["arm"] == "p", "start"] = 0
    centro = pd.merge(centro, maxseg, how="left", on="chr").dropna()
    centro.loc[centro["arm"] == "q", "end_x"] = centro["end_y"]
    centro = centro.rename(columns={"end_x":"end"})
    centro = centro.drop(columns=["band", "annot", "end_y"]).reset_index(drop=True)
    # chr_arm for later groupbys
    centro["chr_arm"] = centro[["chr", "arm"]].apply(lambda x: "".join(x), axis=1)
    centro = centro.rename(columns={"start":"arm_start", "end":"arm_end"})
    # type checks, if no x and y pandas turns chr to ints
    centro["arm_start"] = centro["arm_start"].astype(dtype=int, errors="ignore")
    centro["arm_end"] = centro["arm_end"].astype(dtype=int, errors="ignore")
    centro["chr"] = centro["chr"].astype(dtype=str, errors="ignore")
    return centro;

def get_chr_arm_positions(chr_arm, chr_arms, pt_samples):
    chr_arm_start = chr_arms.loc[chr_arms['chr_arm'] == chr_arm]['arm_start'].tolist()[0]
    chr_arm_end = chr_arms.loc[chr_arms['chr_arm'] == chr_arm]['arm_end'].tolist()[0]
    return chr_arm, pd.DataFrame(columns = pt_samples, index = np.arange(chr_arm_start,chr_arm_end,1))

def split_row(intersect, chr_arms, idx):
    row = intersect.iloc[idx]
    chr_arms = chr_arms[chr_arms['chr'] == row['chromosome']]
    arm_break = chr_arms['arm_end'].tolist()[0]
    insert = pd.DataFrame(index = [idx+0.25,idx-0.25], 
                          columns = intersect.columns.tolist())
    insert['chromosome'] = row['chromosome']
    for sample in samples: 
        insert[sample] = row[sample]
    insert['start.pos'] = [row['start.pos'],arm_break]
    insert['end.pos'] = [arm_break, row['end.pos']]
    insert['chr_arm'] = [row['chromosome']+'p', row['chromosome']+'q']
    intersect = intersect.append(insert).sort_index()
    return intersect;

# =============================================================================
# Create patient column
# =============================================================================

sample_name_map = pd.read_csv("/groups/wyattgrp/users/amurtha/Ghent/pten_breakpoint_analysis/WES_serverNameMap.tsv", sep = '\t')

sample_name_map['Patient'] = sample_name_map['Sample ID correct'].str.split('_').str[1]

# =============================================================================
# Get pt samples
# =============================================================================

samples = sample_name_map[sample_name_map['Patient'] == pt]
samples = samples['Sample ID server']

# =============================================================================
# Import intersection and sample level cn data
# =============================================================================

intersect = pd.read_csv("/groups/wyattgrp/users/amurtha/Ghent/data/wxs_segments_sequenza/%s_intersect.tsv" % pt, 
                        sep = '\t', header = None, names = ['chromosome','start.pos','end.pos', 'cn'])
intersect = intersect.drop('cn', axis = 1)

bed_loc = "/groups/wyattgrp/users/amurtha/Ghent/data/wxs_segments_sequenza/"
cn_dict = dict(zip(samples,
                   [pd.read_csv(bed_loc + "%s_segments.txt.BED" % sample, 
                                sep = '\t', 
                                header = None, 
                                names = ['chromosome','start.pos','end.pos', 'cn']) 
                    for sample in samples]))

# =============================================================================
# Loop over segments; loop over samples
# =============================================================================

for index, seg in intersect.iterrows():
    chrom = seg['chromosome']
    seg_start = seg['start.pos']
    seg_end = seg['end.pos']
    for sample in samples:
        df = cn_dict.get(sample).copy()
        df = df[df['chromosome'] == chrom]
        df = df[(df['start.pos'] <= seg_start)&(df['end.pos'] >= seg_end)]
        if len(df) == 1:
            cn = df.reset_index().iloc[0]['cn']
            intersect.at[index, sample] = cn
        del df

del cn_dict
            
# =============================================================================
# drop Y chromosome
# =============================================================================

intersect = intersect[intersect['chromosome'] != 'Y']

# =============================================================================
# Merge chomosome segment onto intersect dataframe
# =============================================================================

centro = pd.read_csv("/groups/wyattgrp/reference/cytobands/cytoBandhg38.txt", sep="\t", low_memory=False, names=["chr", "start", "end", "band", "annot"], dtype={"chr":str, "start":np.int32, "end":np.int32, "band":str, "annot":str})

chr_arms = make_centro(centro)

for index, row in intersect.iterrows():
    tmp = chr_arms.copy()
    chrom = row['chromosome']
    start = row['start.pos']
    end = row['end.pos']
    tmp = tmp[tmp['chr'] == chrom]
    tmp = tmp[(tmp['arm_start'] <= start)&(tmp['arm_end'] >= end)]
    if len(tmp) != 1:
        split_row(intersect, chr_arms, index)
    else:
        intersect.at[index, 'chr_arm'] = tmp['chr_arm'].tolist()[0]
  
chr_arms = chr_arms[chr_arms['chr_arm'].isin(intersect['chr_arm'])]
chr_arms['length'] = chr_arms['arm_end'] - chr_arms['arm_start']
      
# =============================================================================
# Create plot
# =============================================================================

fig = plt.figure()

gs = gridspec.GridSpec(len(chr_arms), len(samples), height_ratios = chr_arms['length'])

