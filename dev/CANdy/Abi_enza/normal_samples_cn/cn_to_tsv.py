# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:07:11 2020

@author: amurtha
"""
import pandas as pd

def meltCN(filepath):
    df = pd.read_csv(filepath, delimiter = '\t', index_col=None)
    df = pd.melt(df, id_vars=['GENE', 'CHROMOSOME', 'START', 'END'])
    df.rename(columns={'value': 'Copy_num'}, inplace=True)
    df.rename(columns={'variable': 'Sample ID'}, inplace=True)
    df['Cohort'] = 'Abienza'
    df['Patient ID'] = df['Sample ID'].str.split('-').str[:2].str.join('-')
    df['Log_ratio'] = df['Copy_num'].str.split(':').str[1]
    df['Copy_num'] = df['Copy_num'].str.split(':').str[0]
    df[['Copy_num','Log_ratio']] = df[['Copy_num','Log_ratio']].apply(pd.to_numeric)
    df = df[['Cohort','Patient ID','Sample ID', 'GENE', 'Copy_num', 'Log_ratio', 'CHROMOSOME', 'START', 'END']]
    return df;

abienza = meltCN('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Abi-enza_log_ratios.tsv')

abienza.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Abi-enza_logratios_melted.tsv', sep = '\t', index = None)



