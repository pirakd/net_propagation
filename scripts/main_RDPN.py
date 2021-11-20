import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
from utils import get_propagation_input, save_propagation_score, get_time
from statistic_methods import get_sample_p_values, bh_correction
from propagation_routines import propagate_network, propagate_networks, propagate_networks_parallel
from args import Args, MockArgs

test_name = path.basename(__file__).split('.')[0]
args = MockArgs(test_name)
n_processes = 2

#Read the network
network_graph = utils.read_network(args.network_file_path)

#Load prior set
prior_set, prior_data, _ = args.experiment_reader(args)

prior_gene_dict = utils.convert_symbols_to_ids(prior_set)
prior_set_ids = list(prior_gene_dict.values())
propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type, network_graph)

# since edges of random graphs are one we set all original edges to one
for u, v, d in network_graph.edges(data=True):
    d[2] = 1

#Use the graph, either run the propagation or load previously acquired propagation results
_, _, genes_id_to_idx, gene_scores, self_propagation = propagate_network(network_graph, propagation_input, args=args)
genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

# Propagate using randomized networks
if n_processes == 1:
    _, random_networks_scores, self_propagations = propagate_networks(network_graph, args, list(genes_id_to_idx.keys()), prior_set_ids,
                                                   propagation_input, args.random_networks_dir, n_networks=args.n_networks)
else:
    random_networks_scores, self_propagations = propagate_networks_parallel(network_graph, args, list(genes_id_to_idx.keys()), prior_set_ids,
                                                   propagation_input, args.random_networks_dir, n_networks=args.n_networks, n_processes=n_processes)

n_random_networks = random_networks_scores.shape[0]
# Rank the genes in the original network compared to the random networks
p_values, _, ranks = get_sample_p_values(gene_scores, random_networks_scores, two_tailed=True)
adj_p_val = bh_correction(p_values)
file_name = '{}_{}_{}_{}_{}_RDPN'.format(args.propagation_input_type, args.experiment_name, args.sheet_name,
                                    args.condition_function_name, str(args.alpha))

save_propagation_score(propagation_scores=gene_scores, prior_set=prior_set, propagation_input=propagation_input,
                       genes_idx_to_id=genes_idx_to_id, args=args, self_propagation=self_propagation,
                       randomization_ranks=ranks, n_randomizations=n_random_networks , scores_p_values=p_values, file_name=file_name,)
