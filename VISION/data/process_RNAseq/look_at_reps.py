#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import itertools


abundances_df = pd.read_csv("txi.kallisto.abundances.csv", index_col=0)
new_abundances_df = pd.DataFrame(index=abundances_df.index)
full_column_names = np.array(abundances_df.columns)
column_names = np.array([ct for exp, ct in (name.split('_') for name in abundances_df.columns)])
'''compare the replicates'''
ct_rep_OI = ['CMP', 'ERY', 'GMP', 'G1E', 'ER4', 'CFUE', 'IMK', 'CFUMK', 'LSK', 'MEP', 'CLP', 'MON', 'B', 'CD4', 'CD8', 'NK', 'NEU'] #these are ones with definite IDEAS data
straight_scatter = []
combo_scatters = []
k=0
for ct in ct_rep_OI:
    where_col = np.where(column_names == ct)[0]
    if len(where_col) > 2:
        '''There are 4 of these (each with 3 reps). let's average them and store them in the new df. let's also scatter plot among the combinations'''
        where_col_names = [full_column_names[where_col_val] for where_col_val in where_col]
        combo_scatters.append([])
        to_set = list(itertools.combinations(where_col_names, 2))
        combo_scatters[k] = to_set
        k+=1
        new_abundances_df[ct] = abundances_df[full_column_names[where_col]].mean(axis=1)

    elif len(where_col) == 2:
        '''there are 8 of these. let's average them and store them in the new df. but let's also scatter plot them'''
        new_abundances_df[ct] = abundances_df[full_column_names[where_col]].mean(axis=1)
        straight_scatter.append((full_column_names[where_col][0], full_column_names[where_col][1]))
    else: #they only show up one time
        '''there are 5 of these. Let's just add them to the new df'''
        new_abundances_df[ct] = abundances_df[full_column_names[where_col]]

print(new_abundances_df)

def add_trendline(x, y, ax):
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    ax.plot(x,p(x),color='dodgerblue', linestyle='--', linewidth=0.5)
    ax.text(-6, 14, "y=%.3fx+(%.3f)"%(z[0], z[1]), color='dodgerblue')

def add_one_one_line(x_s, y_s, ax):
    ax.plot(x_s, y_s, color='darkgoldenrod', linestyle='--', linewidth=0.5)
    ax.text(-6, 9.5, "one-to-one", color='darkgoldenrod')

def add_correlation(x, y, ax):
    pearsonr, pval = stats.pearsonr(x, y)
    ax.text(-6, 11.75, "Pearsonr: %.3f"%pearsonr, color='saddlebrown')

fig, axes = plt.subplots(3, 4, figsize=(18,18))
fig.suptitle('Gene level expression of biological replicates\nPlotting log2(TPM+0.001)', fontsize=16)
for k, combo_scatter in enumerate(combo_scatters): #This will loop four times
    for l in range(len(combo_scatter)):
        col_0, col_1 = combo_scatter[l][0], combo_scatter[l][1]
        x, y = np.log2(abundances_df[col_0]+0.001), np.log2(abundances_df[col_1]+0.001)
        min_ax_val = np.minimum(np.amin(x), np.amin(y))
        max_ax_val = np.maximum(np.amax(x), np.amax(y))
        if l == 0:
            axes[l, k].set_title(col_0.split('_')[1], fontsize=14, fontweight='bold')
        axes[l, k].scatter(x,y, alpha=0.2, s=1, color='black')
        add_trendline(x, y, axes[l, k])
        add_one_one_line(np.linspace(min_ax_val, max_ax_val, 200), np.linspace(min_ax_val, max_ax_val, 200), axes[l, k])
        add_correlation(x, y, axes[l, k])
        axes[l, k].set_xlabel(col_0, fontsize=12)
        axes[l, k].set_ylabel(col_1, fontsize=12)
        axes[l, k].set_xlim(min_ax_val, max_ax_val)
        axes[l, k].set_ylim(min_ax_val, max_ax_val)
        axes[l, k].set_box_aspect(1)
    plt.tight_layout()
    fig.savefig('scatter_rna_rep_all3ct.png')
    plt.close(fig)

fig, axes = plt.subplots(4, 2, figsize=(12,12))
fig.suptitle('Gene level expression of biological replicates\nPlotting log2(TPM+0.001)')
for i in range(len(straight_scatter)):
    col_0, col_1 = straight_scatter[i][0], straight_scatter[i][1]
    x, y = np.log2(abundances_df[col_0]+0.001), np.log2(abundances_df[col_1]+0.001)
    min_ax_val = np.minimum(np.amin(x), np.amin(y))
    max_ax_val = np.maximum(np.amax(x), np.amax(y))
    m, n = int(np.floor(i/2)), i%2
    axes[m, n].set_title(col_0.split('_')[1], fontsize=12)
    axes[m, n].scatter(x, y, alpha=0.2, s=1, color='black')
    add_trendline(x, y, axes[m, n])
    add_one_one_line(np.linspace(min_ax_val, max_ax_val, 200), np.linspace(min_ax_val, max_ax_val, 200), axes[m, n])
    add_correlation(x, y, axes[m, n])
    axes[m, n].set_xlabel(col_0, fontsize=12)
    axes[m, n].set_ylabel(col_1, fontsize=12)
    axes[m, n].set_xlim(min_ax_val, max_ax_val)
    axes[m, n].set_ylim(min_ax_val, max_ax_val)
    axes[m, n].set_box_aspect(0.75)
plt.tight_layout()
fig.savefig('scatter_rna_rep_all2ct.png')
plt.close(fig)
#
# '''look at some informed comparisons'''
# '''--> ERY, ERYA, & ERYB'''
# '''--> LSK, HSC, LTHSC'''
