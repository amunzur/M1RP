# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:23:57 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

# =============================================================================
# Constants
# =============================================================================

gene = 'FOXA1'
sample = 'Ghent_ID14_MLN4'
tc_sample = 'ID14_C14_MLN4'

# =============================================================================
# Import data
# =============================================================================

genes = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/refGene.txt', sep = '\t', 
                    low_memory = False, 
                    names = ['bin','name','chr', 'strand', 'txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], 
                    usecols = ['bin','name','chr','strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts','exonEnds','score', 'name2', 'cdsStartStat', 'cdsEndStat','exonFrames'], 
                    dtype = {'bin': np.int32,'name': np.str,'chr': np.str,'strand': np.str, 'txStart': np.int32, 'txEnd': np.int32, 'cdsStart': np.int32, 'cdsEnd': np.int32, 'exonCount': np.int32, 'exonStarts': np.str,'exonEnds': np.str,'score': np.int32, 'name2': np.str, 'cdsStartStat': np.str, 'cdsEndStat': np.str,'exonFrames': np.str})

snps = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks/'+sample+'_hetz_snp.igv', sep = '\t')
lrs = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Copy number analysis - All/Copy number plots (whole exome)/igv_tracks/'+sample+'_logratio.igv', sep = '\t')
tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/tumor_fraction.tsv', sep = '\t')

lrs = lrs.rename(columns = {sample: 'log_ratio'})
snps = snps.rename(columns = {sample: 'vaf'})

tc = tc[tc['Sample'] == tc_sample].reset_index()
tc = tc['tumor_fraction'][0]

# =============================================================================
# Keep only gene
# =============================================================================

genes = genes[genes['name2'] == gene]

# =============================================================================
# Get chr, start and end of gene. Add 100000 base pairs to get a larger pickture of region
# =============================================================================

chrom = genes['chr'].unique().tolist()[0]
start = genes['txStart'].min() - 1000000
end = genes['txEnd'].max() + 1000000

# =============================================================================
# Keep copy number data within start and end and on the correct chromosome
# =============================================================================

lrs = lrs[lrs['CHROM'] == chrom]
lrs = lrs[lrs['START'] >= start]
lrs = lrs[lrs['END'] <= end]

snps = snps[snps['CHROM'] == chrom]
snps = snps[snps['START'] >= start]
snps = snps[snps['END'] <= end]

lr_avg = lrs['log_ratio'].median()
copies = 1 * (1 + (2**lr_avg - 1) / tc) if chrom == 'chrX' else 2 * (1 + (2**lr_avg - 1) / tc)

snp_avg = snps['vaf'].median()

print(lr_avg, copies)
print(snp_avg)