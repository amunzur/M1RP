import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
from matplotlib import gridspec

# Input data
cnv_comp = pd.read_csv("https://docs.google.com/spreadsheets/d/1zi7UiPteDA3VU4jtcp1dxET4UGCMKvtLqab-7G13_oU/export?format=csv&gid=1176347443")

# Plotting overall concordance

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot()
ax.scatter(cnv_comp['Log_ratio_Targeted'],cnv_comp['Log_ratio_WES'],s=2)
# for index, row in cnv_comp.iterrows():
#     if -1.5 < cnv_comp.at[index,'WES_take_two'] < 0.1 and  cnv_comp.at[index,'Log_ratio_Targeted'] < -1:
#         ax.annotate(cnv_comp.at[index,"Gene"], (cnv_comp.at[index,'Log_ratio_Targeted'],cnv_comp.at[index,'WES_take_two']),fontsize=5)
ax.set_xlim(-4.5,2.5)
ax.set_ylim(-4.5,2.5)
ax.set_xlabel("Targeted Logratio")
ax.set_ylabel("WES Logratio")
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, c='#A3A2A2', ls='--', alpha=0.75, zorder=0)


# Plotting per patient comparison
colors = ['r','b','#318f5d','#FFC300','#8314e3',"#FC8404",'k','#f0027f','#13bed1','#d113bb','#798a7c','#30a9ff']
patients = set(cnv_comp['Patient ID'])

for p in patients:
    p_data = cnv_comp[cnv_comp['Patient ID'].isin([p])]
    groups = p_data.groupby(['Sample ID'])
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot()
    i = 0
    for name, group in groups:
        ax.scatter(group.Log_ratio_Targeted, group.Log_ratio_WES, color = colors[i],s=8,label=name)
        i = i + 1
    ax.legend(loc='lower right')
    ax.set_title(p)
    ax.set_ylabel('WES logratio')
    ax.set_xlabel('Targeted logratio')

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, c='#A3A2A2', ls='--', alpha=0.75, zorder=0)