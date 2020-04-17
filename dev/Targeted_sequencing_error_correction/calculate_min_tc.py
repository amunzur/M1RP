# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 13:27:05 2020

@author: amurtha
"""

import pandas as pd
import scipy.stats as stats

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t')

min_p = 0.001

copy_change = -1
cn['min_tc_1loss'] = (2**(-1*stats.norm.ppf(1 - min_p) * cn['lr_std'] + cn['lr_mean']) - 1)/(((cn['ploidy']+copy_change).clip(lower = 0)/cn['ploidy']) - 1)

copy_change = -2
cn['min_tc_2loss'] = (2**(-1*stats.norm.ppf(1 - min_p) * cn['lr_std'] + cn['lr_mean']) - 1)/(((cn['ploidy']+copy_change).clip(lower = 0)/cn['ploidy']) - 1)


copy_change = 1
cn['min_tc_1gain'] = (2**(1*stats.norm.ppf(1 - min_p) * cn['lr_std'] + cn['lr_mean']) - 1)/(((cn['ploidy']+copy_change).clip(lower = 0)/cn['ploidy']) - 1)

copy_change = 2
cn['min_tc_2gain'] = (2**(1*stats.norm.ppf(1 - min_p) * cn['lr_std'] + cn['lr_mean']) - 1)/(((cn['ploidy']+copy_change).clip(lower = 0)/cn['ploidy']) - 1)

cn.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/Targeted_sequencing_error_correction/cn_melted_withError.tsv', sep = '\t', index = None)