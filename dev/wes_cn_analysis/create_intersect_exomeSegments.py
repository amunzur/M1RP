# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:48:25 2020

@author: amurtha
"""

import pandas as pd
import subprocess

samples = pd.read_csv("/groups/wyattgrp/users/amurtha/Ghent/pten_breakpoint_analysis/WES_serverNameMap.tsv", sep  = '\t')

samples['Patient'] = samples['Sample ID correct'].str.split('_').str[1]

for pt in samples['Patient'].unique():
    pt_samples = samples[samples['Patient'] == pt]
    sample_list = pt_samples['Sample ID server']
    sample_list = sample_list + '_segments.txt.BED'
    sample_list = sample_list.to_list()
    cmd = "bedtools intersect -a %s -b %s > %s_intersect.tsv" % (sample_list[0], ' '.join(sample_list[1:]), pt)
    # cmd = ['bedtools','intersect', '-a', sample_list[0]] + sample_list.tolist()[1:] + 
    subprocess.run(cmd, shell = True, cwd = '/groups/wyattgrp/users/amurtha/Ghent/data/wxs_segments_sequenza', check = True, text = True)
    