# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:45:39 2020

@author: amurtha
"""

import pandas as pd
import os

path = None
file_type = '.igv'

if file_type == '.pdf':
    path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)'
if file_type == '.igv':
    path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks'

    files = os.listdir(path)

files = [file.split('.')[0] for file in files if file_type in file]
files = [file.split('_logratio')[0] for file in files if '_logratio' in file]

original = files

files = [file.split('_WES')[0] for file in files]

files = pd.Series(files)

for index, val in files.iteritems():
    if 'Ghent_' in val:
        files.at[index] = val.split('Ghent_')[1]
        val = files.at[index]
    if '_C' in val:
        files.at[index] = '_'.join([val.split('_')[0],val.split('_')[2]])
        
        
df = pd.DataFrame({'orig':original, 'new': files})

if file_type == '.pdf':
    for index, row in df.iterrows():
        original = path+'/'+row['orig']+file_type
        output = path+'/M1RP_'+row['new']+file_type
        try:
            os.rename(original,output)
        except WindowsError:
            continue

if file_type == '.igv':
    for index, row in df.iterrows():
        original = path+'/'+row['orig']+'_logratio'+file_type
        output = path+'/M1RP_'+row['new']+'_logratio'+file_type
        try:
            os.rename(original,output)
        except WindowsError:
            continue
        original = path+'/'+row['orig']+'_hetz_snp'+file_type
        output = path+'/M1RP_'+row['new']+'_hetz_snp'+file_type
        try:
            os.rename(original,output)
        except WindowsError:
            continue
        

# =============================================================================
# To rename columns in .igv files to updated names
# =============================================================================

    
path = 'C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks'

for file in os.listdir(path):
    df = pd.read_csv('/'.join([path,file]), sep = '\t')
    col = file
    col = col.split('_logratio.igv')[0]
    col = col.split('_hetz_snp.igv')[0]
    df.columns = df.columns.tolist()[0:-1]+[col]
    df.to_csv('/'.join([path,file]), sep = '\t', index = None)