import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.patches import Rectangle, Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def plot_enrichment_table(enrichment_table, adj_p_mat, direction, interesting_pathways, save_dir=None, experiment_names=None,
                          title=None, res_type=None, minimum_pathway_log_p_value_for_display=None):
    fig, ax = plt.subplots()

    enriched_clusters = np.nonzero(np.sum(enrichment_table, axis=0) != 0)[0]
    found_pathways = np.nonzero(np.sum(enrichment_table, axis=1) != 0)[0]
    enrichment_table = enrichment_table[:, enriched_clusters][found_pathways, :]
    interesting_pathways_filtered = {x: xx for x, xx in enumerate(interesting_pathways) if x in found_pathways}
    
    # set values lower than a given threshold to zero:
    if minimum_pathway_log_p_value_for_display:
        enrichment_table[abs(enrichment_table)<minimum_pathway_log_p_value_for_display] = 0.0

    annotation_map = (np.round(enrichment_table, 3)).astype(str)
    annotation_map[annotation_map == '0.0'] = ''

    y_ticks = [x[:60] for x in interesting_pathways_filtered.values()]
    important_indexes = np.where(adj_p_mat < 0.05)

    if direction is not None and res_type != 'z_score':
        # set low propagation scores to be negative in order to color them blue
        enrichment_table[np.logical_not(direction)] = -enrichment_table[np.logical_not(direction)]
        # enrichment_table[np.logical_not(direction[:, 0]), 0] = -enrichment_table[np.logical_not(direction[:, 0]), 0]
        # enrichment_table[np.logical_not(direction[:, 2]), 2] = -enrichment_table[np.logical_not(direction[:, 2]), 2]
        # enrichment_table[np.logical_not(direction[:, 3]), 3] = -enrichment_table[np.logical_not(direction[:, 3]), 3]
        # enrichment_table[np.logical_not(direction[:, 5]), 5] = -enrichment_table[np.logical_not(direction[:, 5]), 5]

    # set color bar size and location
    cax = inset_axes(ax,
                     width="100%",  # width: 40% of parent_bbox width
                     height="70%",  # height: 10% of parent_bbox height
                     loc='lower left',
                     bbox_to_anchor=(1.1, 0.2, 0.2, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )
    colorbar_edge = np.maximum(np.abs(np.min(enrichment_table)), np.abs(np.max(enrichment_table)))
    heatmap = sns.heatmap(enrichment_table, fmt=".4s", yticklabels=y_ticks,
                                 cbar_ax=cax, annot=annotation_map, cmap="coolwarm",
                                 linewidths=.1, linecolor='gray',
                                    cbar_kws={'label': res_type}, ax=ax, vmin=np.min([0, -colorbar_edge]), vmax=colorbar_edge)
    ax.set_xticklabels(experiment_names, fontsize=24)

    # circle significant scores (<0.05)
    for i in range(len(important_indexes[0])):
        heatmap.add_patch(Rectangle((important_indexes[1][i], important_indexes[0][i]), 1, 1, fill=False, edgecolor='black', lw=3))
    # set font size of colobar ticks
    heatmap.figure.axes[-1].yaxis.label.set_size(20)

    # set colorbar tick values to be positive in both direction
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(np.round(np.abs(cbar.get_ticks()),1))

    # rotate test names
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')

    # add square color legend
    legend_handles = [Patch(fill =False, edgecolor='black', label='adj_p_value<0.05')]

    if direction is not None:
        legend_handles.append(Patch(fill =True, color='red', label='Upwards'))
        legend_handles.append(Patch(fill =True, color='blue', label='Downwards'))

    # locate legend
    ax.legend(handles=legend_handles, bbox_to_anchor=[1, 0, 1, 1], ncol=1, loc='lower left', fontsize=14)
    ax.set_title(title, fontsize=28)
    sns.set(font_scale=0.9)
    fig.set_size_inches(8.5, 18.5, forward=True)
    plt.savefig(save_dir, bbox_inches='tight')
    sns.set(font_scale=1/0.9)