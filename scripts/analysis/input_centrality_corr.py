import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.dirname(path.realpath(__file__)))))
import utils as utils
from utils import get_propagation_input, save_propagation_score, moving_average
from propagation_routines import propagate_network, generate_similarity_matrix, propagate
from args import MockArgs, DeltaArgs, CovPhosJanArgs
import numpy as np
load_scores = False
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns
font_size_1 = 20
font_size_2 = 16
def plot_correlation(args, network_graph):
    args.propagation_input_type = 'ones'
    args.get_derived_parameters(is_create_output_folder=False)
    matrix, genes = generate_similarity_matrix(network_graph, None)
    num_genes = len(genes)
    gene_id_to_idx = dict([(gene, index) for (index, gene) in enumerate(genes)])
    genes_idx_to_id = {xx: x for x, xx in gene_id_to_idx.items()}
    prior_set, prior_data, all_data = args.experiment_reader(args)
    prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
    all_genes = list(prior_data.Gene_Name)
    propagation_input = get_propagation_input(prior_gene_dict, prior_data, SCORE_TYPE,
                                          network=network_graph)
    ones_input = get_propagation_input(prior_gene_dict, prior_data, 'ones', network=network_graph)
    one_scores_10h = propagate([id for id in propagation_input.keys()], ones_input, matrix, gene_id_to_idx, num_genes, args)
    degree_dict = dict(network_graph.degree(weight=2))
    centrality_10h = np.array([one_scores_10h[gene_id_to_idx[id]] for id in propagation_input.keys()])
    degree_10h = np.array([degree_dict[id] for id in propagation_input.keys()])
    scores_arr_10h = np.array([propagation_input[id] for id in propagation_input.keys()])

    fig, ax = plt.subplots()

    y_value = np.abs(scores_arr_10h) if Graph_Y_AXIS=='abs' else scores_arr_10h
    spear_cor_10h = scipy.stats.spearmanr(centrality_10h, y_value)

    # m, b = np.polyfit(np.log(centrality_10h), np.log(np.abs(scores_arr_10h)), deg=1)
    # line_val = lambda x: np.exp((np.log(x) * m) + b)
    # y = line_val(x)

    sorted = np.argsort(centrality_10h)
    move_mean = moving_average(y_value[sorted], 50)
    plt.plot(centrality_10h[sorted], move_mean, color= 'orange', label='Moving average')

    y_label = '|log2FC|' if Graph_Y_AXIS=='abs' else 'log2FC'
    sns.scatterplot(centrality_10h, y_value, alpha=0.34, ax=ax)

    # where some data has already been plotted to ax
    handles, labels = ax.get_legend_handles_labels()

    patch_1 = mpatches.Patch(label='Spearman p_value={0:.2e}'.format(spear_cor_10h.pvalue))
    patch_2 = mpatches.Patch(label='Spearman correlation={0:.2f}'.format(spear_cor_10h.correlation))
    handles.append(patch_1)
    handles.append(patch_2)

    plt.legend(handles=handles, loc='upper right')

    x = np.linspace(np.min(centrality_10h), np.max(centrality_10h), 100)

    plt.title(args.sheet_name )
    plt.xlabel('One propagation score', fontsize=font_size_2)
    plt.ylabel(y_label, fontsize=font_size_2)
    plt.yscale('symlog')
    plt.xscale('log')
    plt.grid(alpha=1, color='black')
    # plt.legend()
    plt.tight_layout()
    plt.savefig(path.join(args.output_folder, args.sheet_name))
    plt.close()

if __name__ == '__main__':
    SCORE_TYPE = 'Score'
    Graph_Y_AXIS = 'signed' # abs, signed
    sheet_names = ['Mock_10h-Mock_24h', 'DeltaE_10h-VIC_10h', 'DeltaE_24h-VIC_24h',
                   'Alpha_10h-VIC_10h', 'Alpha_24h-VIC_24h', 'Beta_10h-VIC_10h', 'Beta_24h-VIC_24h',
                   'Gamma_10h-VIC_10h', 'Gamma_24h-VIC_24h',]
    test_name = path.basename(__file__).split('.')[0]
    args = CovPhosJanArgs(test_name, is_create_output_folder=True)
    output_folder = args.output_folder
    network_graph = utils.read_network(args.network_file_path)
    test_name = path.basename(__file__).split('.')[0]
    for sheet in sheet_names:
        args.sheet_name = sheet
        args.get_derived_parameters(is_create_output_folder=False)
        plot_correlation(args, network_graph)
