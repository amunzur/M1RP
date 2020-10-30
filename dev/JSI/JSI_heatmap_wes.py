# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:51:37 2020

@author: amurtha
"""

import seaborn as sns
import matplotlib.gridspec
import pandas as pd
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

class SeabornFig2Grid():

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        elif isinstance(self.sg, sns.matrix.ClusterGrid):
            self._moveclustergrid()
        # self._finalize()
        
    def _moveclustergrid(self):
        """ Move Clustergrid """
        self.subgrid = gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=self.subplot, height_ratios =[1,3])
        self._moveaxes(self.sg.ax_col_dendrogram, self.subgrid[0,0])
        self._moveaxes(self.sg.ax_heatmap, self.subgrid[1,0])

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())
        
# =============================================================================
# Constants
# =============================================================================

pt = 'ID10'
cohort = 'M1RP'

# =============================================================================
# Import matrix
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/wxs_muts_all.tsv', sep = '\t')
matrix_wes = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/jsi_matrices/jsi_matrix_wxs.tsv', sep = '\t', index_col = 'Unnamed: 0')
matrix_tngs = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/jsi_matrices/jsi_matrix_mutationOnly.tsv', sep = '\t', index_col = 'Unnamed: 0')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[tc['Cohort'] == cohort]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Final tNGS_TC'] > .1]

muts = muts[muts['Sample ID'].isin(tc['Sample ID'])]



for pt in muts['Patient ID'].unique().tolist():
# for pt in ['ID10']:
    print(pt)
    samples = muts[muts['Patient ID'] == pt]['Sample ID'].unique().tolist()
    if len(samples) <= 2:
        continue

    matrix_wes_tmp = matrix_wes[matrix_wes.index.isin(samples)][samples].fillna(1)
    matrix_tngs_tmp = matrix_tngs[matrix_tngs.index.isin(samples)][samples].fillna(1)
    
    # =============================================================================
    # Removing M1RP_IDX from columns and rows for clarity
    # =============================================================================
    for i, col in enumerate(matrix_wes_tmp.columns):
        if 'cfDNA' not in col:
            col = col.split("_")[2]
        else:
            col = col.split("_")[2] + col.split("_")[3]
        matrix_wes_tmp = matrix_wes_tmp.rename(columns={matrix_wes_tmp.columns[i]: col})
        
    matrix_wes_tmp = matrix_wes_tmp.reset_index()
    for index, row in matrix_wes_tmp.iterrows():
        name = matrix_wes_tmp.at[index,'index']
        if 'cfDNA' not in name:
            matrix_wes_tmp.at[index,'index'] = name.split("_")[2]
        else:
            matrix_wes_tmp.at[index,'index'] = name.split("_")[2] + name.split("_")[3]
    matrix_wes_tmp = matrix_wes_tmp.set_index('index')
    
    for i, col in enumerate(matrix_tngs_tmp.columns):
        if 'cfDNA' not in col:
            col = col.split("_")[2]
        else:
            col = col.split("_")[2] + col.split("_")[3]
        matrix_tngs_tmp = matrix_tngs_tmp.rename(columns={matrix_tngs_tmp.columns[i]: col})
        
    matrix_tngs_tmp = matrix_tngs_tmp.reset_index()
    for index, row in matrix_tngs_tmp.iterrows():
        name = matrix_tngs_tmp.at[index,'index']
        if 'cfDNA' not in name:
            matrix_tngs_tmp.at[index,'index'] = name.split("_")[2]
        else:
            matrix_tngs_tmp.at[index,'index'] = name.split("_")[2] + name.split("_")[3]
    matrix_tngs_tmp = matrix_tngs_tmp.set_index('index')
    
    # =============================================================================
    # Create and plot per patient matrix
    # =============================================================================
    
    # ploting heatmap and dendrogram
    w = sns.clustermap(matrix_wes_tmp,
                       figsize=(4,4),
                       cmap="BuPu",
                       cbar_pos=None,
                       square = True)
    w.ax_row_dendrogram.set_visible(False)
    
    ax = w.ax_heatmap
    ax.set_ylabel("")
    plt.setp(w.ax_heatmap.get_xticklabels(), rotation=90,fontsize=8)
    plt.setp(w.ax_heatmap.get_yticklabels(), rotation=0, fontsize = 8)
    w.ax_col_dendrogram.set_title(pt + ' wes')
    
    samples = list()
    a = w.ax_heatmap.get_yticklabels()
    for x in a:
        samples.append(x.get_text())
        
    # ploting heatmap and dendrogram
    t = sns.clustermap(matrix_tngs_tmp,
                       figsize=(4,4),
                       cmap="BuPu",
                       cbar_pos=None,
                       square = True)
    t.ax_row_dendrogram.set_visible(False)
    
    ax = t.ax_heatmap
    ax.set_ylabel("")
    plt.setp(t.ax_heatmap.get_xticklabels(), rotation=90,fontsize=8)
    plt.setp(t.ax_heatmap.get_yticklabels(), rotation=0, fontsize = 8)
    t.ax_col_dendrogram.set_title(pt + ' targeted')
    
    samples = list()
    a = t.ax_heatmap.get_yticklabels()
    for x in a:
        samples.append(x.get_text())
    
    
    fig = plt.figure(figsize=(6,4))
    gs = gridspec.GridSpec(1, 2)
    
    mg0 = SeabornFig2Grid(w, fig, gs[0,0])
    mg1 = SeabornFig2Grid(t, fig, gs[0,1])
    fig.subplots_adjust(wspace = 1, right = 1.2, hspace = 0, bottom = 0.3)
    
    fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/wes_matrices/%s_%s_wes_jsiMatrix.pdf' % (cohort, pt), bbox_to_inches = 'tight')
    fig.set_size_inches(12,8)
    fig.subplots_adjust(wspace = 0.5, right = 2.8, hspace = 0, bottom = 0.55, top = 1.9)
    plt.setp(fig.axes[1].get_xticklabels(), rotation=90,fontsize=16)
    plt.setp(fig.axes[1].get_yticklabels(), rotation=0,fontsize=16)
    plt.setp(fig.axes[3].get_xticklabels(), rotation=90, fontsize=16)
    plt.setp(fig.axes[3].get_yticklabels(), rotation=0, fontsize=16)
    plt.setp(fig.axes[0].set_title(pt+ ' wes'), fontsize = 16)
    plt.setp(fig.axes[2].set_title(pt+ ' targeted'), fontsize = 16)
    fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/wes_matrices/%s_%s_wes_jsiMatrix.png' % (cohort, pt))