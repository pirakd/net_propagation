""""
This script is used to generate and save propagation score for a variety 
of parameters
"""
import sys
from os import path, makedirs
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
from utils import read_prior_set, get_propagation_input, save_propagation_score, get_time
from propagation_routines import propagate_network
from args import CovArgs, HDArgs, MockArgs

propagation_results_dir = path.join('output', 'propagation_results')
args = MockArgs(None, is_create_output_folder=False)
# args = HDArgs(None, False)
n_tests = 1
alpha = [0.9] * n_tests
sheet_names = ['Kent_10h-VIC_10h',
                'Kent_24h-VIC_24h']
sheet_names = ['Table_A'] * len(alpha)
#
args.remove_self_propagation = True
prior_set_conditions = [args.condition_function_name] * n_tests # ['huntington_DDA_significant'] * len(alpha)
propagation_input_type_list = ['ones'] * n_tests


network_graph = utils.read_network(args.network_file_path)
fc_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list())
prop_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list())

# Run:
for c in range(n_tests):
    # loading arguments
    args.sheet_name = sheet_names[c]
    args.alpha = alpha[c]
    args.condition_function_name = prior_set_conditions[c]
    args.set_condition_function()
    args.propagation_input_type = propagation_input_type_list[c]

    # loading prior set
    prior_set, prior_data, _ = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)
    prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
    prior_set_ids = set.intersection(set(prior_gene_dict.values()), set(network_graph.nodes))
    propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type,
                                              network=network_graph)
    # Propagate network
    _, _, genes_id_to_idx, gene_scores, self_propagation = propagate_network(network_graph, propagation_input, args)
    genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

    prop_score_file_name = '{}_{}_{}_{}'.format(args.propagation_input_type, args.sheet_name, prior_set_conditions[c], str(args.alpha))

    # save propagation score
    save_propagation_score(file_name=prop_score_file_name, propagation_scores=gene_scores, prior_set=prior_set,
                           propagation_input=propagation_input, self_propagation=self_propagation, genes_idx_to_id=genes_idx_to_id, args=args,
                           date=get_time())