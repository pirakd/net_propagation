import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
from utils import read_prior_set, get_propagation_input, save_propagation_score, get_time
from statistic_methods import get_sample_p_values, bh_correction
from propagation_routines import propagate_network, propagate_networks, propagate_networks_parallel
from args import Args
import pickle as pl
test_name = 'main_RDPN'
args = Args(test_name)
n_processes = 1

#Read the network
network_graph = utils.read_network(args.network_file_path)


#Load prior set
prior_set, prior_data, _ = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)

prior_gene_dict = utils.convert_symbols_to_ids(prior_set)
prior_set_ids = list(prior_gene_dict.values())
propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type)

# since edges of random graphs are one we set all original edges to one
for u, v, d in network_graph.edges(data=True):
    d[2] = 1

#Use the graph, either run the propagation or load previously acquired propagation results
_, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, propagation_input, args=args, prior_set=prior_set_ids)

genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

# Propagate using randomized networks
if n_processes == 1:
    _, random_networks_scores = propagate_networks(network_graph, args, list(genes_id_to_idx.keys()), prior_set_ids,
                                                   propagation_input, args.random_networks_dir, n_networks=args.n_networks)
else:
    random_networks_scores = propagate_networks_parallel(network_graph, args, list(genes_id_to_idx.keys()), prior_set_ids,
                                                   propagation_input, args.random_networks_dir, n_networks=args.n_networks, n_processes=n_processes)

prop_score_file_name = args.sheet_name + '_' + args.condition_function_name + '_' + str(args.alpha)
# save propagation score
save_propagation_score(file_name=prop_score_file_name, propagation_scores=gene_scores, prior_set=prior_set,
                       propagation_input=propagation_input, genes_idx_to_id=genes_idx_to_id, args=args,
                       date=get_time(), save_dir=args.output_folder)

# Rank the genes in the original network compared to the random networks
p_values = get_sample_p_values(gene_scores, random_networks_scores)
adj_p_val = bh_correction(p_values)

title = ['gene_id\tp_value\tadj_p\n']
lines = ['{}\t{}\t{}\n'.format(genes_idx_to_id[i], p_values[i], adj_p_val[i]) for i in range(len(p_values))]
lines = title + lines

with open(path.join(args.output_folder, 'p_values.txt'), 'w') as f:
    f.writelines(lines)