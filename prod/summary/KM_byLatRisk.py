# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:59:49 2021

@author: amurtha
"""


import pandas as pd
import numpy as np
import lifelines
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
# from lifelines.utils import datetimes_to_durations
# from lifelines import NelsonAalenFitter
# from lifelines.utils import survival_table_from_events
from lifelines import AalenAdditiveFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
from lifelines.plotting import add_at_risk_counts



clin = pd.read_excel("C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx")

# =============================================================================
# Calculate number of months between collection and last follow up
# =============================================================================

clin = clin[clin['mHSPC_adt_administered'] == 1]
clin.loc[clin['crpc_status'] == 0, 'crpc_date'] = clin['dolf']

clin['Time_to_CRPC'] = ((clin['crpc_date'] - clin['mHSPC_adt_start']) / np.timedelta64(1, 'M'))
clin['Time_to_CRPC'] = clin['Time_to_CRPC'].astype(float)

# =============================================================================
# Create Categorical data from ctDNA faction
# =============================================================================


lat_0 = clin[clin['latitude_raw'] < 2]
lat_23 = clin[clin['latitude_raw'] >= 2]

# =============================================================================
# Prepare plots
# =============================================================================

fig1, ax = plt.subplots(1, figsize=(2.5,2.5))

# =============================================================================
# Plot ctDNA below median on kmf1
# =============================================================================

kmf1 = KaplanMeierFitter()
color = 'grey'
defective_patient_number = str(len(lat_0))
defective_label = str(str('Low risk ')+r"(n="+defective_patient_number+")")
T = lat_0['Time_to_CRPC'].round(3)
C = lat_0['crpc_status'].astype(np.int32)
kmf1.fit(T, event_observed = C, label = defective_label)
kmf1.plot(ax=ax,show_censors = True, ci_show = False, color = color, lw = 1)
# lat_0_median=kmf1.median

# =============================================================================
# Plot ctDNA above median on kmf2
# =============================================================================

kmf3 = KaplanMeierFitter()
color = 'red'
defective_patient_number = str(len(lat_23))
defective_label = str(str('High risk ')+r"(n="+defective_patient_number+")")
T = lat_23['Time_to_CRPC'].round(3)
C = lat_23['crpc_status'].astype(np.int32)
kmf3.fit(T, event_observed = C, label = defective_label)
kmf3.plot(ax=ax,show_censors = True, ci_show = False, color = color, lw = 1)
# lat_23_median=kmf3.median


# =============================================================================
# Run Cox progression
# =============================================================================

df_test = clin[['Time_to_CRPC','crpc_status','latitude_raw']].copy()
df_test['latitude_raw'] = df_test['latitude_raw'].replace({3:2})
df_test = df_test.dropna()
cph = CoxPHFitter()
cph.fit(df_test,'Time_to_CRPC',event_col = 'crpc_status', show_progress = True, step_size = 0.001)
coxPHsummary = cph.summary
coxPHsummary['Variable'] = coxPHsummary.index

# =============================================================================
# Adjust plots aethetics
# =============================================================================

ax.plot([0,0],[0,0],color = 'w',alpha=0,label = 'p = '+str(round(coxPHsummary.p[0],3)))
legend = ax.legend(fontsize = 6, loc = 'upper right')
plt.setp(legend.get_title(),fontsize=7)

ax.set_xlim(0,60)
ax.set_xticks(np.arange(0,61,12))

ax.set_ylabel('Survival fraction', fontsize=8)
ax.set_xlabel('ADT to CRPC (months)', fontsize=8)
ax.tick_params(axis='y',which="both", length=5, width=0.5, color='k', direction="out",
               rotation=0, left=True, reset=False, labelleft=True, labelsize=8, labelrotation=0)
ax.tick_params(axis='x',which="both", length=5, width=0.5, color='k', direction="out",
               rotation=0, left=True, reset=False, labelleft=True, labelsize=8, labelrotation=0)
ax.grid(b=False)
ax.grid(axis='y', color='0.85', linewidth=0.5, linestyle='dotted', zorder=-1)
ax.grid(axis='x', color='0.85', linewidth=0.5, linestyle='dotted', zorder=-1, clip_on=False)

# add_at_risk_counts(kmf1, kmf2, kmf3, ax=ax, rows_to_show =['At risk'])

plt.tick_params(axis='both', which='major', labelsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

plt.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/summary/KM_byLatRisk.pdf')