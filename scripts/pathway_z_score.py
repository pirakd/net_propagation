import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, bh_correction
from os import path
from propagation_routines import propagate_network, get_genes_p_values
from visualization_tools import plot_enrichment_table
import numpy as np
from args import Args
import pickle as pl

prior_set_conditions = ['kent_vic_24h', 'kent_mock_no_vic_mock_24h',  'kent_vic_10h']
# prior_set_conditions = ['kent_mock_no_vic_mock_24h'] * 3
test_name = 'classical_enrichment'
args = Args(test_name)
title = 'Protein Abundance'
n_draws = 5000
p_vals = list()
adj_pvals = list()
direction = list()
# Read the h_sapiens network
network_graph = utils.read_network(args.network_file)

for c, condition in enumerate(prior_set_conditions):
    args.condition_function_name = prior_set_conditions[c]
    args.set_condition_function()

    # loading prior set
    prior_set, prior_data = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)
    prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
    prior_set_ids = list(prior_gene_dict.values())
    propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type)

    # Using the graph, either run the propagation or load previously acquired propagation results
    _, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, propagation_input,
                                                           prior_set=list(prior_gene_dict.values()))
    genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}


    # load pathways and their related genes
    interesting_pathways = load_interesting_pathways(args.interesting_pathway_file_dir)
    genes_of_interesting_pathways = load_pathways_genes(args.pathway_file_dir, interesting_pathways)

    # check scores of genes related to interesting pathways
    genes_idx_by_pathways = [[genes_id_to_idx[y] for y in x if y in genes_id_to_idx] for x in list(genes_of_interesting_pathways.values())]
    genes_scores_by_patways = [gene_scores[x] for x in genes_idx_by_pathways]

    total_mean = np.mean(gene_scores)
    total_sd = np.std(gene_scores, ddof=1)

    # get p_value empirically
    mean_scores = [np.mean(x) for x in genes_scores_by_patways]
    mean_stds = [np.std(x, ddof=1) for x in genes_scores_by_patways]
    random_scores = []
    for p in range(len(genes_idx_by_pathways)):
        random_scores.append([])
        random_scores[-1] = np.array([np.mean(gene_scores[np.random.randint(0,len(gene_scores), len(genes_idx_by_pathways[p]))]) for x in range(n_draws)])
    random_scores = np.array(random_scores).transpose()
    p_temp, direction_temp = get_genes_p_values(mean_scores, random_scores, two_tailed=True)
    p_vals.append(p_temp)
    direction.append(direction_temp)
    adj_pvals.append(bh_correction(p_vals[-1]))

direction = np.array(direction).transpose()
adj_pvals = np.array(adj_pvals).transpose()
with open(path.join(args.output_folder, 'enrichemnt_results.pl'), 'wb') as f:
    pl.dump({'direction': direction, 'adj_p': adj_pvals}, f)

fig_out_dir = path.join(args.output_folder, 'abundance_classical_enrichment')
plot_enrichment_table(-np.log10(adj_pvals), direction, interesting_pathways, fig_out_dir,
                      experiment_names=prior_set_conditions, title=title)