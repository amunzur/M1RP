# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:14:14 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

gene_err = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/normal_samples_cn/normal_samples_cn.tsv', sep = '\t', index_col = 'GENE')

tp53 = gene_err.loc['TP53']

lr_dist = [np.random.normal(tp53.lr_mean, tp53.lr_std) for i in range(1000000)]

fig, ax = plt.subplots()
ax.hist(lr_dist, bins = 700)

ax.plot([tp53.lr_mean, tp53.lr_mean], [0,10000])
ax.plot([stats.norm.ppf(0.001, loc = tp53.lr_mean, scale = tp53.lr_std),stats.norm.ppf(0.001, loc = tp53.lr_mean, scale = tp53.lr_std)], [0,10000])
ax.plot([stats.norm.ppf(1-0.001, loc = tp53.lr_mean, scale = tp53.lr_std),stats.norm.ppf(1-0.001, loc = tp53.lr_mean, scale = tp53.lr_std)], [0,10000])


ax.get_yaxis().set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(-0.7,0.6)
ax.set_ylim(0,6000)
ax.set_xlabel('Log-ratio')
ax.set_title('TP53')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/normal_samples_cn/TP53_normal_distribution.png', dpi = 200)