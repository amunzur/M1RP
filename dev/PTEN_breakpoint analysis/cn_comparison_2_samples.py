# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:12:46 2020

@author: amurtha
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sample_name_map = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Exome Analysis (Elie)/WES_serverNameMap.tsv', sep = '\t').set_index('Sample ID correct')

s1='M1RP_ID10_MB1'
s2='M1RP_ID10_PB2'

s1_v2 = sample_name_map.at[s1, 'Sample ID server']
s2_v2 = sample_name_map.at[s2, 'Sample ID server']

# =============================================================================
# Import segment files
# =============================================================================

s1_segments = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Exome Analysis (Elie)/sequenza_dbox/%s_segments.txt' % s1_v2, sep = '\t')
s2_segments = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Exome Analysis (Elie)/sequenza_dbox/%s_segments.txt' % s2_v2, sep = '\t')

# =============================================================================
# Adjust chromosome column
# =============================================================================

s1_segments['chromosome'] = s1_segments['chromosome'].replace({'X':23})
s1_segments = s1_segments[s1_segments['chromosome'] != 'Y']
s1_segments['chromosome'] = s1_segments['chromosome'].astype(np.float64)

s2_segments['chromosome'] = s2_segments['chromosome'].replace({'X':23})
s2_segments = s2_segments[s2_segments['chromosome'] != 'Y']
s2_segments['chromosome'] = s2_segments['chromosome'].astype(np.float64)

chr_list = s1_segments['chromosome'].append(s2_segments['chromosome']).unique().tolist()

# =============================================================================
# split any segments that do not share start or ends
# =============================================================================

segments = s1_segments[['chromosome','start.pos','end.pos','CNt']]
segments = segments.rename(columns = {'CNt':'CN1'})

#return true and the copy number of that segment if there is a exact match for the current segment in s2.segments. Otherwise return false and None. Return type is tuple
def check_overlap(seg, s2_segments):
    chrom = seg.chromosome
    start = seg['start.pos']
    end = seg['end.pos']
    s2_segments = s2_segments[s2_segments['chromosome'] == chrom]
    s2_segments = s2_segments[(s2_segments['start.pos'] == start)&(s2_segments['end.pos'] == end)]
    if len(s2_segments) == 1:
        return (True, s2_segments.reset_index().at[0,'CNt']); 
    else:
        return (False, None);

# Return True is there is a matching segment in s2 that shares start.pos with seg and the s2 segment is shorter than the s1 segment. Return false if not
def check_share_start_s1Longer(seg, s2_segments):
    chrom = seg.chromosome
    start = seg['start.pos']
    end1 = seg['end.pos']
    s2_segments = s2_segments[s2_segments['chromosome'] == chrom]
    s2_segments = s2_segments[(s2_segments['start.pos'] == start)&(s2_segments['end.pos'] < end1)]
    if len(s2_segments) == 0:
        return False;
    else:
        return True;
    
# Return true if there is a shared start and s2 is longer than s1. Otherwise return false
def check_share_start_s2Longer(seg, s2_segments):
    chrom = seg.chromosome
    start = seg['start.pos']
    end1 = seg['end.pos']
    s2_segments = s2_segments[s2_segments['chromosome'] == chrom]
    s2_segments = s2_segments[(s2_segments['start.pos'] == start)&(s2_segments['end.pos'] > end1)]
    if len(s2_segments) == 0:
        return False;
    else:
        return True;

# Return segments adjusted with a split in the s1 segment that shares a start point and extents past the corrosponding s2 segment. 
def split_s1_shared_start(index, seg1, segments, s2_segments):
    print(index)
    chrom = seg1.chromosome
    start = seg1['start.pos']
    end1 = seg1['end.pos']
    cn1 = seg1['CN1']
    s2_segments = s2_segments[s2_segments['chromosome'] == chrom]
    s2_segments = s2_segments[(s2_segments['start.pos'] == start)&(s2_segments['end.pos'] < end1)].reset_index()
    seg2 = s2_segments.iloc[0]
    end2 = seg2['end.pos']
    start2 = seg2['start.pos']
    cn2 = seg2['CNt']
    insert = pd.DataFrame({'chromosome':[chrom,chrom],'start.pos':[start,start2],'end.pos':[end2,end1],'CN1':[cn1,cn1],'CN2':[cn2,np.nan]}, index = [index-0.25,index+0.25])
    segments = segments.drop(index, axis = 0).append(insert).sort_index()
    # print(85)
    return segments;

# Return s2_segments adjusted with a split in the s2 segment that shares a start point and extents past the corrosponding s1 segment. 
def split_s2_shared_start(s1_seg, s2_segments):
    chrom = s1_seg['chromosome']
    start1 = s1_seg['start.pos']
    end1 = s1_seg['end.pos']
    s2_seg = s2_segments[(s2_segments['chromosome'] == chrom)&(s2_segments['start.pos'] == start1)].copy()
    index = s2_seg.index.tolist()[0]
    s2_seg = s2_segments.iloc[int(index)]
    end2 = s2_seg['end.pos']
    start2
    cn = s2_seg['CNt']
    insert = pd.DataFrame({'chromosome':[chrom,chrom],'start.pos':[start1,start2],'end.pos':[end1,end2], 'CNt':[cn,cn]}, index = [index-0.25, index+0.25])
    s2_segments = s2_segments.drop(index, axis = 0).append(insert).sort_index()
    return s2_segments

    
    
finished_indexes = []

while len(finished_indexes) < len(segments):
    for index, segment in segments.iterrows():
        if index in finished_indexes: continue;
        overlap = check_overlap(segment, s2_segments)
        if overlap[0]:
            segments.at[index, 'CN2'] = overlap[1]
            finished_indexes.append(index)
        elif check_share_start_s1Longer(segment, s2_segments):
            segments = split_s1_shared_start(index, segment, segments, s2_segments)
            break;
        elif check_share_start_s2Longer(segment, s2_segments):
            s2_segments = split_s2_shared_start(segment, s2_segments)
            break;
        
        finished_indexes.append(index)
        
print(segments.count())
    

    
    