# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:03:37 2020

@author: amurtha
"""

import pandas as pd
import datetime

clin = pd.read_excel('C:/Users/amurtha/Dropbox/Ghent M1 2019/Anonymized clinical data M1RP + sample annotation.xlsx', sheet_name = 'M1-cohort')

clin = clin[clin.index < 39]

clin['Date CRPC'] = clin['Date CRPC'].replace('-', datetime.datetime(1900,1,1))
clin['Date RP'] = clin['Date RP'].replace ('St-Lucas', datetime.datetime(1900,1,1))

clin['RP_to_CRCP'] = (clin['Date CRPC'] - clin['Date RP']).dt.days / 30