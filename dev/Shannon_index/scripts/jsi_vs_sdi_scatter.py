# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 18:12:26 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from matplotlib.patches import Rectangle

jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/patient_mean_jsi.tsv', sep = '\t')
sdi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/Shannon_index/scripts/targeted_sdi.tsv', sep = '\t')

jsi = jsi.rename(columns = {'Unnamed: 0':'Patient ID'})
sdi = sdi.rename(columns = {'Unnamed: 0':'Patient ID'})

jsi = jsi.merge(sdi, on = 'Patient ID')


jsi_iqr = (np.percentile(np.array(jsi['JSI']),10),np.percentile(np.array(jsi['JSI']),90))
sdi_iqr = (np.percentile(np.array(jsi['sdi']),10),np.percentile(np.array(jsi['sdi']),90))


fig,ax = plt.subplots()

ax.scatter(jsi['JSI'], jsi['sdi'], c = 'k', s = 10, clip_on = False)
ax.add_patch(Rectangle((0,0), 1, sdi_iqr[0], color = 'k', alpha = 0.25))
ax.add_patch(Rectangle((0,0), jsi_iqr[0], 3, color = 'k', alpha = 0.25))
ax.add_patch(Rectangle((0,sdi_iqr[1]), 1, 2, color = 'k', alpha = 0.25))
ax.add_patch(Rectangle((jsi_iqr[1],0), 1, 3, color = 'k', alpha = 0.25))

ax.set_ylabel('Shannon diversity index')
ax.set_xlabel('Jaccard similarity index')
ax.set_xlim(right = 1)
ax.set_ylim(bottom = 0)

lin = stats.linregress(jsi['JSI'],jsi['sdi'])

ax.text(x = 1.0, y = 2.5, s = 'r = %.2f, p = %.3f' % (lin[2],lin[3]), ha = 'right')

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/Shannon_index/scripts/jsi_vs_sdi_scatter.pdf')