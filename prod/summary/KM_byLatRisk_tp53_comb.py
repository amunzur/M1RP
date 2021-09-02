# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:07:16 2021

@author: amurtha
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:05:22 2021

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
df_TS = pd.read_csv(r'G:\Evan MSc\Ghent M1\Fall2020_Updated_Analysis\Clean_data\automated_input_TS_oncoprint_primary_mar2021.tsv', sep='\t')

# =============================================================================
# Calculate number of months between collection and last follow up
# =============================================================================

clin = clin[clin['mHSPC_adt_administered'] == 1]
clin.loc[clin['crpc_status'] == 0, 'crpc_date'] = clin['dolf']

clin['Time_to_CRPC'] = ((clin['crpc_date'] - clin['date_rp']) / np.timedelta64(1, 'M'))
clin['Time_to_CRPC'] = clin['Time_to_CRPC'].astype(float)

clin = clin.merge(df_TS[['PATIENT','TP53']], left_on = 'patient_id', right_on = 'PATIENT', how = 'left')
clin.loc[clin['TP53'] < 2, 'TP53'] = 0
clin['TP53'] = clin['TP53'].fillna(0)

# =============================================================================
# Create Categorical data from ctDNA faction
# =============================================================================

lat_0_intact = clin[(clin['latitude_raw'] < 2)&(clin['TP53'] != 2)]

lat_23_hit = clin[(clin['latitude_raw'] >= 2)|(clin['TP53'] == 2)]

# =============================================================================
# Prepare plots
# =============================================================================

fig1, ax = plt.subplots(1, figsize=(2.5,2.5))

# =============================================================================
# Plot ctDNA below median on kmf1
# =============================================================================

kmf2 = KaplanMeierFitter()
color = 'grey'
defective_patient_number = str(len(lat_0_intact))
defective_label = str(str('Low risk & TP53 intact ')+r"(n="+defective_patient_number+")")
T = lat_0_intact['Time_to_CRPC'].round(3)
C = lat_0_intact['crpc_status'].astype(np.int32)
kmf2.fit(T, event_observed = C, label = defective_label)
kmf2.plot(ax=ax,show_censors = True, ci_show = False, color = color, lw = 1)
# lat_0_median=kmf1.median


# =============================================================================
# Plot ctDNA above median on kmf2
# =============================================================================

kmf3 = KaplanMeierFitter()
color = 'purple'
defective_patient_number = str(len(lat_23_hit))
defective_label = str(str('High risk or TP53 null ')+r"(n="+defective_patient_number+")")
T = lat_23_hit['Time_to_CRPC'].round(3)
C = lat_23_hit['crpc_status'].astype(np.int32)
kmf3.fit(T, event_observed = C, label = defective_label)
kmf3.plot(ax=ax,show_censors = True, ci_show = False, color = color, lw = 1)
# lat_23_median=kmf3.median


# =============================================================================
# Run Cox progression
# =============================================================================

df_test = clin[['Time_to_CRPC','crpc_status','latitude_raw', 'TP53']].copy()
df_test['latitude_raw'] = df_test['latitude_raw'].replace({3:2})

df_test['risk_comb'] = 1
df_test.loc[(df_test['TP53'] == 2)|(df_test['latitude_raw'] >= 2), 'risk_comb'] = 2
df_test = df_test.drop(['latitude_raw','TP53'], axis = 1)
df_test = df_test.dropna()
cph = CoxPHFitter()
cph.fit(df_test,'Time_to_CRPC',event_col = 'crpc_status', show_progress = True, step_size = 0.001)
coxPHsummary = cph.summary
coxPHsummary['Variable'] = coxPHsummary.index

# =============================================================================
# Adjust plots aethetics
# =============================================================================

ax.plot([0,0],[0,0],color = 'w',alpha=0,label = 'p = '+str(round(coxPHsummary.p[0],3)))
legend = ax.legend(fontsize = 5, loc = 'lower left', )
plt.setp(legend.get_title(),fontsize=7)

ax.set_xlim(0,60)
ax.set_xticks(np.arange(0,61,12))

ax.set_ylabel('Survival fraction', fontsize=8)
ax.set_xlabel('RP to CRPC (months)', fontsize=8)
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

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/KM_byLatRisk_TP53_comb.pdf')
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/KM_byLatRisk_TP53_comb.png')