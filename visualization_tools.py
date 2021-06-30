import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.patches import Rectangle, Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def plot_enrichment_table(enrichment_table, direction, interesting_pathways, save_dir=None, experiment_names=None, title=None):
    fig, ax = plt.subplots()

    enriched_clusters = np.nonzero(np.sum(enrichment_table, axis=0) != 0)[0]
    found_pathways = np.nonzero(np.sum(enrichment_table, axis=1) != 0)[0]
    enrichment_table = enrichment_table[:, enriched_clusters][found_pathways, :]
    interesting_pathways_filtered = {x: xx for x, xx in enumerate(interesting_pathways) if x in found_pathways}
    annotation_map = (np.round(enrichment_table, 3)).astype(str)
    annotation_map[annotation_map == '0.0'] = ''
    y_ticks = [x[:60] for x in interesting_pathways_filtered.values()]
    important_indexes = np.where(enrichment_table > 1.3)

    # set low propagation scores to be negative in order to color them blue
    enrichment_table[np.logical_not(direction)] = -enrichment_table[np.logical_not(direction)]

    # set color bar size and location
    cax = inset_axes(ax,
                     width="100%",  # width: 40% of parent_bbox width
                     height="70%",  # height: 10% of parent_bbox height
                     loc='lower left',
                     bbox_to_anchor=(1.1, 0.2, 0.2, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )

    heatmap = sns.heatmap(enrichment_table, fmt=".4s", yticklabels=y_ticks, xticklabels=experiment_names,
                                 cbar_ax=cax, annot=annotation_map, cmap="coolwarm",
                                 linewidths=.1, linecolor='gray',
                                    cbar_kws={'label': '-log10(adj_p)'}, square=True, ax=ax)

    # circle significant scores (<0.05)
    for i in range(len(important_indexes[0])):
        heatmap.add_patch(Rectangle((important_indexes[1][i], important_indexes[0][i]), 1, 1, fill=False, edgecolor='blue', lw=3))
    # set font size of colobar ticks
    heatmap.figure.axes[-1].yaxis.label.set_size(20)

    # set colorbar tick values to be positive in both direction
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(np.round(np.abs(cbar.get_ticks()),1))

    # rotate test names
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')

    # add square color legend
    legend_handles = [Patch(fill =False, edgecolor='blue', label='adjusted p_value<0.05'),
                      Patch(fill =True, color='red', label='high propagation score'),
                      Patch(fill =True, color='blue', label='low propagation score')]
    # locate legend
    ax.legend(handles=legend_handles, bbox_to_anchor=[1, 0, 1, 1], ncol=1, loc='lower left', fontsize=14)

    fig.set_size_inches(18.5, 10.5, forward=True)
    plt.savefig(save_dir, bbox_inches='tight')
