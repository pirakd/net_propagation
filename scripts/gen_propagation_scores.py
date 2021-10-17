"""
Example of generating a single propagation process
"""
import sys
from os import path, makedirs
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
from utils import read_prior_set, get_propagation_input, save_propagation_score, get_time
from propagation_routines import propagate_network
from args import MockArgs

propagation_results_dir = path.join('output', 'propagation_results')
args = MockArgs(None, is_create_output_folder=False)

network_graph = utils.read_network(args.network_file_path)

# loading prior set
prior_set, prior_data, _ = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)
prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
prior_set_ids = set.intersection(set(prior_gene_dict.values()), set(network_graph.nodes))
propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type,
                                          network=network_graph)
# Propagate network
_, _, genes_id_to_idx, gene_scores, self_propagation = propagate_network(network_graph, propagation_input, args)
genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}


# save propagation score
save_propagation_score(propagation_scores=gene_scores, prior_set=prior_set,
                       propagation_input=propagation_input, self_propagation=self_propagation,
                       genes_idx_to_id=genes_idx_to_id, args=args,
                       date=get_time(), file_name=None)
