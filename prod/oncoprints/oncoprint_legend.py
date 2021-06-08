# -*- coding: utf-8 -*-
"""
Created on Mon May 31 13:26:36 2021

@author: amurtha
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def p(c):
    return Patch(facecolor = c, lw = 0)

def l(c,m='s'):
    return Line2D([],[], lw = 0, color = c, marker = m, markeredgewidth=0)

NE = Patch(facecolor='#E6E7E8',hatch='\\\\\\\\\\\\',edgecolor = 'w', lw = 0)

handles = [p('w'),p('#EE2D24'),p('#F59496'),p('#E6E7E8'),p('#9CC5E9'),p('#3F60AC'),NE,
           p('w'),l('#79B443'),l('#FFC907'),l('#BD4398'),l('k'),l('k',m='*'),
           p('w'),p('green'),p('yellow'),p('purple'),p('red'),]

labels = [r'$\bf{Copy\ number}$','Amplification','Gain','Copy neutral','Loss','Deep deletion','Not evaluable',
          r'$\bf{Mutations}$','Missense','Truncating','In-frame InDel','Somatic','Germline',
          r'$\bf{Sample\ type}$','Primary biopsy','Radical prostatectomy','Metastatic tissue','cfDNA',]


fig,ax = plt.subplots(figsize = (1, 2.5))
ax.set_visible(False)

fig.legend(handles, labels, fontsize = 6, handlelength = 0.8, ncol = 1, frameon = False, loc = 'upper left', bbox_to_anchor = (-0.06,1))

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Oncoprints/oncoprint_legend_vertical.pdf')