# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:29:47 2020

@author: amurtha
"""

import pandas as pd
import numpy as np


# =============================================================================
# Check samples missing from the TC estimation table
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

tf = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tf.columns = tf.iloc[0]
tf = tf.drop(tf.index[0])

#Remove WBC samples
samples = samples[samples['Sample type'] != 'PBMC']
#Remove missing sample names
samples = samples[~samples['Sample ID'].isnull()]
#Remove failed samples
samples = samples[samples['Targeted sequencing'] != 'Failed']
# Remove samples that have not been sequenced
samples = samples[~samples['Targeted sequencing'].isnull()]

tc_missing_samples = [sample for sample in samples['Sample ID'].tolist() if not sample in tf['Sample ID'].tolist()]
samples1_missing_samples = [sample for sample in tf['Sample ID'].tolist() if not sample in samples['Sample ID'].tolist()]



# =============================================================================
# Check samples missing from the sequencing metrics table
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

metrics = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=878831703')

#Remove failed samples
samples = samples[samples['Targeted sequencing'] != 'Failed']
#Remove missing sample names
samples = samples[~samples['Sample ID'].isnull()]
# Remove samples that have not been sequenced
samples = samples[~samples['Targeted sequencing'].isnull()]
# Remove pending samples
samples = samples[samples['Targeted sequencing'] != 'Pending']
# Remove "withoutPrecellys" samples
metrics = metrics[~metrics['Sample ID'].str.contains('WithoutPrecellys')]
# Remove WES for now. will be checked below
metrics = metrics[metrics['Targeted or WES'] != 'Whole exome']

metrics_missing_samples = [sample for sample in samples['Sample ID'].tolist() if not sample in metrics['Sample ID'].tolist()]
samples2_missing_samples = [sample for sample in metrics['Sample ID'].tolist() if not sample in samples['Sample ID'].tolist()]

# =============================================================================
# Check WES for sequencing metrics table
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

metrics = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=878831703')

#Remove failed samples
samples = samples[samples['Targeted sequencing'] != 'Failed']
#Remove missing sample names
samples = samples[~samples['Sample ID'].isnull()]
# Remove pending samples
samples = samples[samples['Whole-exome sequencing'] == 'Completed']
# Remove "withoutPrecellys" samples
metrics = metrics[~metrics['Sample ID'].str.contains('WithoutPrecellys')]
# Remove WES for now. will be checked below
metrics = metrics[metrics['Targeted or WES'] == 'Whole exome']

