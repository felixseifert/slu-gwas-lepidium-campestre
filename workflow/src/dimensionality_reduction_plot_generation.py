# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Generation of plots from dimensionality reduction analysis (MDS, PCA)

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import animation
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

# parse command line arguments
parser = argparse.ArgumentParser(description='Generates plot and 3d animation for dimensionality reduction analysis (MDS, PCA).')
parser.add_argument('--coord', type=str, help='coordinates (tsv file)')
parser.add_argument('--cluster', type=str, help='cluster assignment for colouring (optional)')
parser.add_argument('--eigenval', type=str, help='eigenvalues (tsv file)')
parser.add_argument('--output', type=str, help='output file prefix (without ending)')
parser.add_argument('--samples', type=str, help='sample names, one-per-line (txt file)')
args = parser.parse_args()

if args.coord is None or args.output is None:
    parser.print_help()
    sys.exit()

# data preprocessing
coord_df = pd.read_csv(args.coord, index_col=0, header=0, delimiter='\t')

sample_names_df = None
if args.samples is not None:
    sample_names_df = pd.read_csv(args.samples, index_col=None, header=None, delimiter='\t')

if args.cluster is None:
    coord_df['cluster'] = 1
else:
    cluster_df = pd.read_csv(args.cluster, index_col=None, header=None, delimiter='\t')
    
    if sample_names_df is not None:
        cluster_df.index = sample_names_df[0]
    
    cluster_df.columns = ['cluster']
    coord_df = pd.merge(coord_df, cluster_df, left_index=True, right_index=True)

# import explained variance information
eigenvalue_df = None
if args.eigenval is not None:
    eigenvalue_df = pd.read_csv(args.eigenval, index_col=0, header=0, delimiter='\t')
    explained_variance_df = eigenvalue_df.abs() * 100 / eigenvalue_df.abs().sum()

# 2d plot
x_label = 'PC1'
y_label = 'PC2'

if args.eigenval is not None:
    x_label = 'PC1 (' + str('%.2f' % round(explained_variance_df.iloc[0], 2)) + ' %)'
    y_label = 'PC2 (' + str('%.2f' % round(explained_variance_df.iloc[1], 2)) + ' %)'

colors = ["#2121d9", "#9999ff", "#df0101", "#04b404", "#fffb23", "#ff9326", "#a945ff", "#0089b2"]
color_map = sns.color_palette(colors)

sns_plot = sns.relplot(x=coord_df['PC1'], y=coord_df['PC2'], alpha=0.5, edgecolor='none', hue=coord_df['cluster'], legend=False, palette=color_map, s=10)
sns_plot.set_axis_labels(x_label, y_label)
sns_plot.savefig(args.output + '_2d_plot.png', dpi=300)
sns_plot.savefig(args.output + '_2d_plot.svg', dpi=300)

for index, line in enumerate(coord_df.index):
    plt.annotate(line, (coord_df['PC1'][index], coord_df['PC2'][index]), color = 'gray', ha = 'left', va = 'center', size = 4, textcoords = 'offset points', xytext = (2, 2))
    
sns_plot.savefig(args.output + '_2d_plot_with_label.png', dpi=300)
sns_plot.savefig(args.output + '_2d_plot_with_label.svg', dpi=300)
plt.clf()

# 3d plot
sns.set_style("whitegrid", {'axes.grid' : False})
fig = plt.figure()

ax = Axes3D(fig)
ax.view_init(elev = 30, azim = 45)
ax.scatter(coord_df['PC1'], coord_df['PC2'], coord_df['PC3'], alpha = 0.75, marker = '.', s = 10)

x_label = 'PC1'
y_label = 'PC2'

if args.eigenval is not None:
    x_label = 'PC1 (' + str('%.2f' % round(explained_variance_df.iloc[0], 2)) + ' %)'
    y_label = 'PC2 (' + str('%.2f' % round(explained_variance_df.iloc[1], 2)) + ' %)'
    z_label = 'PC3 (' + str('%.2f' % round(explained_variance_df.iloc[2], 2)) + ' %)'

ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.set_zlabel(z_label)

plt.savefig(args.output + '_3d_plot.png', dpi=300)
plt.savefig(args.output + '_3d_plot.svg', dpi=300)
plt.clf()

# 3d animation
fig = plt.figure()
ax = Axes3D(fig)

def init():
    # uses labels from 3d plot
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    ax.scatter(coord_df['PC1'].to_numpy(), coord_df['PC2'].to_numpy(), coord_df['PC3'].to_numpy(), marker='.', s=5, c='mediumblue', alpha=0.75)

    return fig,

def animate(i):
    ax.view_init(elev=10., azim=(i - 90))

    return fig,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
anim.save(args.output + '_3d_animation.mp4', dpi=300, bitrate=2500, fps=30, extra_args=['-vcodec', 'libx264'])
plt.close()

# matrix plot
sns_plot = sns.pairplot(coord_df[coord_df.columns[0:6]], diag_kind="kde", plot_kws=dict(alpha=0.5, edgecolor="none"))
sns_plot.map_upper(sns.kdeplot, alpha=0.25, cmap=sns.color_palette("coolwarm", as_cmap=True), levels=[0.1, 0.25, 0.5, 0.75, 1], shade="True", thresh = 0.1)

sns_plot.savefig(args.output + '_matrix_plot.png', dpi=300)
sns_plot.savefig(args.output + '_matrix_plot.svg', dpi=300)
plt.close()

# scree plot
plt.xlabel('number of coordinates')
plt.ylabel('eigenvalue')

plt.plot(range(1,31), eigenvalue_df[1:31], 'r-o')
sns.despine()

plt.savefig(args.output + '_scree_plot.png', dpi=300)
plt.savefig(args.output + '_scree_plot.svg', dpi=300)
plt.clf()

# variance explained plot
cumulative_explained_variance = np.cumsum(explained_variance_df)

plt.xlabel('number of coordinates')
plt.ylabel('cumulative explained variance [%]')

plt.step(range(1, cumulative_explained_variance.size + 1), cumulative_explained_variance.values)
sns.despine()
plt.xlim([0, cumulative_explained_variance.size])
plt.ylim([0, 100])

plt.savefig(args.output + '_cumulative_explained_variance.png', dpi=300)
plt.savefig(args.output + '_cumulative_explained_variance.svg', dpi=300)
plt.clf()
