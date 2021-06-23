import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.patches import Rectangle, Patch

def plot_enrichment_table(enrichment_table, direction, interesting_pathways, save_dir=None, experiment_names=None, title=None):

    plt.figure(figsize=(18, 12))

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

    heatmap = sns.heatmap(enrichment_table, fmt=".4s", yticklabels=y_ticks, xticklabels=experiment_names,
                                 annot=annotation_map, cmap="coolwarm",
                                 linewidths=.1, linecolor='gray',
                                    cbar_kws={'label': '-log10(adj_p)'}, square=True)

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
    # sns.set(font_scale=1)
    # add square color legend
    legend_handles = [Patch(fill =False, edgecolor='blue', label='adjusted p_value<0.05'),
                      Patch(fill =True, color='red', label='high propagation score'),
                      Patch(fill =True, color='blue', label='low propagation score')]
    plt.legend(handles=legend_handles, bbox_to_anchor=[4, 0], ncol=1, loc='lower right', fontsize=14,)

    plt.title(title, fontsize=20)
    plt.savefig(save_dir, bbox_inches='tight')
