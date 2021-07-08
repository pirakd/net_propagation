""""
This script is used to generate and save a relatively small number of propagation scores
"""
import utils as utils
from utils import read_prior_set, get_propagation_input
from os import path, makedirs
from propagation_routines import propagate_network, get_genes_p_values
from args import Args
import pickle as pl

prior_set_conditions = ['kent_vic_24h', 'kent_vic_10h'] * 4
propagation_results_dir = path.join('output', 'propagation_results')
args = Args(None, is_create_output_folder=False)
sheet_names = ['Protein_Abundance','RNA'] * 2
alpha = [0.9, 0.9, 1, 1]


network_graph = utils.read_network(args.network_file)
fc_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list())
prop_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list())
for c, condition in enumerate(prior_set_conditions):

    # loading arguments
    args.sheet_name = sheet_names[c]
    args.alpha = alpha[c]
    args.condition_function_name = prior_set_conditions[c]
    args.set_condition_function()

    # loading prior set
    prior_set, prior_data = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)
    prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
    prior_set_ids = list(prior_gene_dict.values())
    propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type)

    # Propagate network
    _, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, propagation_input, args,
                                                           prior_set=list(prior_gene_dict.values()))
    genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

    # save propagation score
    makedirs(propagation_results_dir, exist_ok=True)
    propagation_results_path = path.join(propagation_results_dir,args.sheet_name + '_' + condition + '_' +
                                         str(args.alpha))

    with open(propagation_results_path, 'wb') as f:
        pl.dump({'gene_idx_to_id': genes_idx_to_id, 'gene_prop_scores' : gene_scores}, f)
