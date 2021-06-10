# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:23:24 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# 
# =============================================================================

wes = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/exome_sequencing_metrics_June2021.tsv', sep = '\t')

# =============================================================================
# Plot coverage
# =============================================================================

fig,[ax1,ax2,ax3,ax4] = plt.subplots(nrows = 4, sharex = True)

ax1.bar(wes['SAMPLE'], wes['TOTAL READS (MILLION)'])

ax2.bar(wes['SAMPLE'], wes['COVERAGE (MEDIAN)'])
ax2.scatter(wes['SAMPLE'], wes['COVERAGE (10%)'], label = '10th%-ile', zorder = 100, s = 10, lw = 0)

ax3.bar(wes['SAMPLE'], wes['Median mutant reads'])
ax3.scatter(wes['SAMPLE'], wes['10% mutant reads'], label = '10th%-ile', zorder = 100, s = 10, lw = 0)

ax4.bar(wes['SAMPLE'], wes['Additional depth needed (median 8 mutant reads)'])
ax4.scatter(wes['SAMPLE'], wes['Additional depth needed (90% 8 mutant reads)'], zorder = 100, s = 10, lw = 0)

ax3.plot([-1,len(wes)+1],[8,8], lw = 0.8, marker = None, ls = 'dashed', zorder = 1000, color = 'k')

# =============================================================================
# 
# =============================================================================

ax1.set_ylabel('Total reads\n(millions)', fontsize = 6)
ax2.set_ylabel('Median coverage', fontsize = 6)
ax3.set_ylabel('Median\nmutant reads', fontsize = 6)
ax4.set_ylabel('Additional depth\nrequired', fontsize = 6)

ax1.tick_params(bottom = False, labelbottom = False, labelsize = 6)
ax2.tick_params(bottom = False, labelbottom = False, labelsize = 6)
ax3.tick_params(bottom = False, labelbottom = False, labelsize = 6)
ax4.tick_params(bottom = False, labelbottom = False, labelsize = 6)

ax1.set_yscale('symlog', base = 10, linthresh = 10)
# ax2.set_yscale('symlog', base = 10, linthresh = 10)
ax3.set_yscale('symlog', base = 10, linthresh = 10)
ax4.set_yscale('symlog', base = 10, linthresh = 10)

ax2.set_ylim(0,100)


ax4.set_xlim(-0.5, len(wes))

fig.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Extraction and Sequencing/Figures/WES_coverage_summary.pdf')