# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:01:33 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Polygon,Patch
import scipy.stats as stats

# =============================================================================
# Baseline ctDNA analysis
# =============================================================================

clin = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=33162023')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Keep only baseline m1rp cfDNA samples
# =============================================================================

clin = clin[clin['Cohort'] == 'M1RP']
samples = samples[samples['Cohort'] == 'M1RP']

samples = samples[samples['Sample Category'] == "cfDNA"]
samples = samples[samples['Status at collection'].isin(['Radical prostatectomy','Pre ADT + abiraterone','Diagnostic biopsy'])]
samples = samples.sort_values(['Patient ID','Status at collection'], ascending = [True,False]).drop_duplicates(['Patient ID'])

tc = tc[tc['Sample ID'].isin(samples['Sample ID'])]
clin = clin[clin['Patient ID'].isin(samples['Patient ID'])]

# =============================================================================
# Merge dataframes together
# =============================================================================

tc = tc.merge(clin[['Patient ID','CHAARTED classification','No. bone metastases','Serum PSA at diagnosis','Time from ADT to CRPC (months)','CRPC',]], on = 'Patient ID', how = 'left')

# =============================================================================
# Set up figure
# =============================================================================

fig,[ax1,ax2,ax3] = plt.subplots(ncols = 3, sharey = True, gridspec_kw={'width_ratios':[1,2,2]})

# =============================================================================
# Plot all baseline ctDNA fractions
# =============================================================================

tc_all = tc.copy()
tc_all['x'] = 1
tc_all['x'] = tc_all['x'].apply(lambda x: x + random.uniform(-0.2,0.2))

ax1.scatter(tc_all['x'], tc_all['Final tNGS_TC'], s = 15, c = 'k', alpha = 0.8, lw = 0)
ax1.boxplot(tc_all['Final tNGS_TC'], showfliers = False)

ax1.set_ylabel('ctDNA fraction')
ax1.set_xticklabels(['All baseline samples\nn=%i' % len(tc_all)])

ax1.set_ylim(-0.02,1)

# =============================================================================
# Plot ctDNA fractions by CHAARTED categories
# =============================================================================

tc_dv = tc.copy()
tc_dv['CHAARTED classification'] = tc_dv['CHAARTED classification'].replace({'Low volume':1, 'High volume':2})

tc_dv['x'] = tc_dv['CHAARTED classification'].apply(lambda x: x + random.uniform(-0.2,0.2))

low_volume = tc_dv[tc_dv['CHAARTED classification'] == 1]['Final tNGS_TC']
high_volume = tc_dv[tc_dv['CHAARTED classification'] == 2]['Final tNGS_TC']

ax2.scatter(tc_dv['x'], tc_dv['Final tNGS_TC'], s = 15, c = 'k', alpha = 0.8, lw = 0)
ax2.boxplot([low_volume, high_volume], showfliers = False)

ax2.set_xlabel('CHAARTED classification')
ax2.set_xticklabels(['Low volume\nn=%i' % len(low_volume),'High volume\nn=%i' % len(high_volume)])

p = stats.mannwhitneyu(low_volume, high_volume, alternative = 'less')
ax2.text(0.75, 0.95, "Mann-Whitney U: One sided\np = %.2f" % p[1])

# =============================================================================
# ctDNA fraction by bone met
# =============================================================================

tc_bone = tc.copy()
tc_bone['No. bone metastases'] = tc_bone['No. bone metastases'].apply(lambda x: 1 if x == 0 else 2)

tc_bone['x'] = tc_bone['No. bone metastases'].apply(lambda x: x + random.uniform(-0.2,0.2))

no_bones = tc_bone[tc_bone['No. bone metastases'] == 1]['Final tNGS_TC']
bones = tc_bone[tc_bone['No. bone metastases'] == 2]['Final tNGS_TC']

ax3.scatter(tc_bone['x'], tc_bone['Final tNGS_TC'], s = 15, c = 'k', alpha = 0.8, lw = 0)
ax3.boxplot([no_bones, bones], showfliers = False)

ax3.set_xlabel('Bone metastasis')
ax3.set_xticklabels(['No\nn=%i' % len(no_bones),'Yes\nn=%i' % len(bones)])

p = stats.mannwhitneyu(no_bones, bones, alternative = 'less')
ax3.text(0.75, 0.95, "Mann-Whitney U: One sided\np = %.2f" % p[1])

fig.tight_layout()

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/cfDNAfraction_breakdown.pdf')

# =============================================================================
# cfDNA by scatter and PSA
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(tc['Serum PSA at diagnosis'], tc['Final tNGS_TC'], s = 15, c = 'k', alpha = 0.8, lw = 0)

ax.set_ylabel('ctDNA fraction')
ax.set_xlabel('Serum PSA at diagnosis')

r = stats.linregress(tc['Serum PSA at diagnosis'], tc['Final tNGS_TC'])[2]
ax.text(150,0.8,"r=%.2f" % r)

fig.tight_layout()
fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/cfDNAfractionVsPSA.pdf')

# =============================================================================
# Plot waterfall plot for time from ADT to CRPC colored by positive ctDNA or not
# =============================================================================

fig,ax = plt.subplots()

tc_wf = tc.copy()
tc_wf['Color'] = 'Grey'
tc_wf.loc[tc_wf['Final tNGS_TC'] > 0, 'Color'] = 'Red'

tc_wf = tc_wf.sort_values(['CRPC','Time from ADT to CRPC (months)'], ascending = True)
tc_wf['y'] = np.arange(0,len(tc_wf),1)

ax.barh(tc_wf['Patient ID'], tc_wf['Time from ADT to CRPC (months)'], color = tc_wf['Color'])
for index, row in tc_wf.iterrows():
    if row['CRPC'] != 1:
        time = row['Time from ADT to CRPC (months)']
        y = row['y']
        p = Polygon([[time,y-0.4],[time,y+0.4],[time+2,y]], color = row['Color'])
        ax.add_patch(p)
        
ax.set_xlabel('Time from ADT to CRPC (months)')

handles = [Patch(color = 'red'),Patch(color = 'grey')]
labels = ['Baseline ctDNA positive','Baseline ctDNA negative']

ax.legend(handles, labels, handlelength = 0.8)

fig.tight_layout()
fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/cfDNAfractionVsTimeToCRPC.pdf')