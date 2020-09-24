# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:19:01 2020

@author: amurtha
"""

import pandas as pd
import itertools
import numpy as np

# =============================================================================
# Import WES samples
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/wxs_muts_all.tsv', sep = '\t')

muts = muts[['Patient ID','Sample ID']].drop_duplicates()
patients = muts['Patient ID'].unique().tolist()
    
# =============================================================================
# Create chromosome arms dataframe
# =============================================================================
 
centro = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/wes_cn_analysis/cytobands/cytoBandhg38.txt', sep="\t", low_memory=False, names=["chr", "start", "end", "band", "annot"], dtype={"chr":str, "start":np.int32, "end":np.int32, "band":str, "annot":str})
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
	return centro

chr_arms = make_centro(centro)

chr_arms = chr_arms[~chr_arms['chr'].isin(['X','Y'])]
chroms = chr_arms[["chr", "chr_arm"]]
chroms["chrom"] = chroms["chr"].replace("X", "23").replace("Y", "24").astype(int)
chroms = chroms.sort_values(["chrom", "chr_arm"])
chr_list = [chrom for chrom in chroms["chr"].drop_duplicates().tolist()]
chr_arm_list = [chrarm for chrarm in chroms["chr_arm"].drop_duplicates().tolist()]

# =============================================================================
# Create dictionary of dataframes in each arm, 
# each dataframe will have an index of all positions in that arm
# =============================================================================

for pt in patients:
    print(pt)
    pt_samples = muts[muts['Patient ID'] == pt]['Sample ID'].unique().tolist()
    arm_dict = {}
    for chr_arm in chr_arm_list:
        chr_arm_start = chr_arms.loc[chr_arms['chr_arm'] == chr_arm]['arm_start'].tolist()[0]
        chr_arm_end = chr_arms.loc[chr_arms['chr_arm'] == chr_arm]['arm_end'].tolist()[0]
        arm_dict[chr_arm] = pd.DataFrame(columns = pt_samples, index = np.arange(chr_arm_start,chr_arm_end,1))