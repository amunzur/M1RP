# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:10:26 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import sys

# =============================================================================
# Constants
# =============================================================================

if len(sys.argv) > 1:
    cohort = sys.argv[1]
    alpha = float(sys.argv[2])
else:
    cohort = 'M1RP'
    alpha = 0.001

# =============================================================================
# Helpers
# =============================================================================

def get_copy_number(row):
    direction = -1 if row['Log_ratio'] < 0 else 1
    if direction == 1: 
        if row['Log_ratio'] - row['lr_mean'] > 1.1:
            return 4;
        else:
            return 3
    elif direction == -1:
        if row['Log_ratio'] - row['lr_mean'] < -1.1:
            return 0;
        else:
            return 1
    

# =============================================================================
# Main
# =============================================================================

'''
Add cn_call column to copy number table. For each row where p_value < alpha, 
create a distribution of the expected p-values for 1 copy change. If p_value 
lies below that distribution, consider the change 2 (or more in case of gains). 
Those changes will be given values of -2 or 4. Other significant p values will
be given the copy number value 1, 3. Copy neutral will be given copy number 2. 
If the tumor content is < the min_tc required to reach the current alpha value, 
copy number is assigned -1.
'''

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_%s_FFPE_cna_abridged.tsv' % cohort, sep = '\t')

cn['Adjusted_copy_num'] = np.nan

for index, row in cn.iterrows():
    if row['p_val'] < alpha:
        cn.at[index, 'Adjusted_copy_num'] = get_copy_number(row)
    elif row['Final tNGS_TC'] < row['min_tc_1loss']:
        cn.at[index, 'Adjusted_copy_num'] = -1
    else:
        cn.at[index, 'Adjusted_copy_num'] = 2
        
cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Targeted_sequencing_error_correction/CANdy_%s_FFPE_cna_abridged.tsv' % cohort, sep = '\t', index = None)