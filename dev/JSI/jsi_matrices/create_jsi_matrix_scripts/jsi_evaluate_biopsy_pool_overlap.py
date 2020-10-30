# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 18:16:42 2020

@author: ewarner
"""

# This script will take the output of the JCI matrix script and plot a simple barplot 
# to determine the median overlap of a random biopsy or an optimal biopsy against the overall
# measure of heterogenety in the patient.

# =============================================================================
# Import Modules
# =============================================================================
from matplotlib import gridspec
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# =============================================================================
# Load dataframes
# =============================================================================
df = pd.read_csv('G:\Evan MSc\Ghent M1\Fall2020_Updated_Analysis\Clean_data\jsi_matrix_only_against_VarPool.tsv', sep='\t', index_col=0)
df2 = pd.read_csv('G:\Evan MSc\Ghent M1\Fall2020_Updated_Analysis\Clean_data\jsi_muts_per_patient.tsv', sep='\t', index_col=0)

# =============================================================================
# Melt table into useable format
# =============================================================================
df_new = df.copy()
df_new = df_new.reset_index()
df_new = pd.melt(df_new, id_vars = 'index')
df_new = df_new.dropna(how='any')
df_new= df_new.replace(to_replace=r'_variantPool', value='', regex=True)
df_new.columns = ('Patient', 'Sample', 'JCI_pool_score')

# =============================================================================
# Calculate Median/Mean/Max JCI pool scores
# =============================================================================
df_new['Median'] = df_new.groupby('Patient')['JCI_pool_score'].transform('median')
df_new['Mean'] = df_new.groupby('Patient')['JCI_pool_score'].transform('mean')
df_new['Max'] = df_new.groupby('Patient')['JCI_pool_score'].transform('max')

# =============================================================================
# Merge with #snv and ctNDA% data
# =============================================================================
df2 = df2.reset_index()
df2 = df2.replace(to_replace=r'_variantPool', value='', regex=True)
df2.columns = ['Patient', 'number_snps_used', 'ctDNA%']

df_summary = pd.merge(df_new, df2, on='Patient', how='left')

df_plot = df_summary[['Patient', 'Max', 'number_snps_used','ctDNA%']]
df_plot = df_plot.drop_duplicates()
df_plot['category'] = 'Optimal Biopsy'
df_plot.columns = ['Patient', 'Value', 'number_snps_used', 'ctDNA%', 'category']

df_average = df_summary[['Patient', 'Mean', 'number_snps_used', 'ctDNA%']]
df_average = df_average.drop_duplicates()
df_average['category'] = 'Mean biopsy'
df_average.columns = ['Patient', 'Value', 'number_snps_used', 'ctDNA%', 'category']

df_plot = pd.concat([df_plot, df_average], axis=0)

del df, df_new, df2
# =============================================================================
# Plot Data
# =============================================================================
(fig, ax) = plt.subplots(figsize=(3, 3))

sns.swarmplot(x='category',
              y='Value',
              data=df_plot,
              size=3,
              color='k')

ax.tick_params(labelrotation=45)

plt.show()





