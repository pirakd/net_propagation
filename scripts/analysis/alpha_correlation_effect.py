import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.dirname(path.realpath(__file__)))))
from statistic_methods import get_sample_p_values
import utils as utils
from utils import get_propagation_input, save_propagation_score
from propagation_routines import generate_similarity_matrix, propagate
from args import MockArgs, DeltaArgs, CovJanArgs, CovPhosJanArgs
import numpy as np
from tqdm import tqdm
from utils import generate_bins
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

test_name = path.basename(__file__).split('.')[0]

args = CovPhosJanArgs(test_name, is_create_output_folder=True)
args.propagation_input_type = 'Score'
args.sheet_name = 'DeltaE_10h-VIC_10h'
args.get_derived_parameters()
title = 'Phospho, {}, {}'.format(args.propagation_input_type, args.sheet_name)

n_randomizations = 1000
precision = 6
min_bin_size = 100
alpha_list = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

network_graph = utils.read_network(args.network_file_path)
prior_set, prior_data, _ = args.experiment_reader(args)
prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
prior_set_ids = set.intersection(set(prior_gene_dict.values()), set(network_graph.nodes))

all_genes_ids_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
all_genes_ids = set.intersection(set(all_genes_ids_dict.values()), set(network_graph.nodes))
propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type,
                                          network=network_graph)
matrix, genes = generate_similarity_matrix(network_graph, None)
num_genes = len(genes)
genes_id_to_idx = dict([(gene, index) for (index, gene) in enumerate(genes)])

ones_input = {id:1 for id in propagation_input.keys()}

ipn_centrality_cor, bin_centrality_cor, ipn_score_cor, bin_score_cor = [], [], [], []
for alpha in alpha_list:
    args.alpha = alpha
    ones_centrality = propagate([id for id in propagation_input.keys()], ones_input, matrix, genes_id_to_idx, num_genes, args)
    ones_centrality_experiment = np.array([ones_centrality[genes_id_to_idx[id]] for id in propagation_input.keys()])

    genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}
    gene_scores = list()
    self_prop_scores = np.zeros(num_genes)

    true_score = propagate([id for id in propagation_input.keys()], propagation_input, matrix, genes_id_to_idx, num_genes, args)
    values = list(propagation_input.values())
    random_scores = list()
    for i in tqdm(range(n_randomizations), desc='Propagating random scores', total=n_randomizations):
        shuffled_values = np.random.choice(values, len(values))
        # values = np.random.choice(values, len(propagation_input))
        shuffled_input = dict(zip(propagation_input.keys(), shuffled_values))
        random_scores.append(propagate([id for id in propagation_input.keys()], shuffled_input, matrix, genes_id_to_idx, num_genes
                                       , args))

    random_scores = np.array(random_scores)
    p_values, _, ranks_ipn = get_sample_p_values(true_score, random_scores, two_tailed=True)


    ones_centrality_quant = np.ceil(ones_centrality_experiment*(10**precision)).astype(int)
    genes_idx_to_id_experiment = {idx: id for idx, id in enumerate(propagation_input.keys())}
    weight_dict = {genes_idx_to_id_experiment[x]: xx for x, xx in enumerate(ones_centrality_experiment)}
    bins, id_to_bin = generate_bins(weight_dict, precision, min_bin_size)

    random_scores = list()
    for i in tqdm(range(n_randomizations), desc='Propagating random scores with bins', total=n_randomizations):
        shuffled_input = {}
        for id in propagation_input.keys():
                shuffled_input[id] = propagation_input[np.random.choice(bins[id_to_bin[id]], 1)[0]]

        random_scores.append(propagate([id for id in propagation_input.keys()], shuffled_input, matrix, genes_id_to_idx, num_genes
                                       , args))
    random_scores = np.array(random_scores)
    p_values, _, ranks_ipn_bin = get_sample_p_values(true_score, random_scores, two_tailed=True)


    logfc_scores = np.array(list(propagation_input.values()))

    ranks_ipn_experiment = np.array([ranks_ipn[genes_id_to_idx[id]] for id in propagation_input])
    ranks_ipn_bin_experiment = np.array([ranks_ipn_bin[genes_id_to_idx[id]] for id in propagation_input])
    ipn_centrality_cor.append(spearmanr(ones_centrality_experiment, ranks_ipn_experiment).correlation)
    bin_centrality_cor.append(spearmanr(ones_centrality_experiment, ranks_ipn_bin_experiment).correlation)
    ipn_score_cor.append(spearmanr(logfc_scores, ranks_ipn_experiment).correlation)
    bin_score_cor.append(spearmanr(logfc_scores, ranks_ipn_bin_experiment).correlation)



font_size = 24
font_size_2 = 18
x_positions = np.arange(1,len(alpha_list)+1)


fig, axs = plt.subplots(2,1)
plt.suptitle(title)
axs[0].set_title('Centrality Correlation', fontsize=font_size)
axs[0].scatter(x_positions, ipn_centrality_cor, label='Permutation correlation')
axs[0].scatter(x_positions, bin_centrality_cor, label='Permutation with binning correlation')
axs[0].set_xlabel(r'Propagation $\alpha$', fontsize=font_size_2)
axs[0].set_ylabel('Spearmean correlation', fontsize=font_size_2)
axs[0].set_xticks(x_positions)
axs[0].set_xticklabels(labels=alpha_list)
axs[0].grid()
axs[0].legend(loc='upper right')
# labels[1] = 'a'
# labels = [a for a in alpha_list]

axs[1].set_title('Raw Scores Correlation', fontsize=font_size)
axs[1].scatter(x_positions, ipn_score_cor, label='Permutation correlation')
axs[1].scatter(x_positions, bin_score_cor, label='Permutation with binning correlation')
axs[1].set_xlabel(r'Propagation $\alpha$', fontsize=font_size_2)
axs[1].set_ylabel('Spearmean correlation', fontsize=font_size_2)
axs[1].set_xticks(x_positions)
axs[1].set_xticklabels(labels=alpha_list)
axs[1].grid()

fig.tight_layout()
fig.set_size_inches(14, 10, forward=True)
plt.savefig(path.join(args.output_folder, test_name))