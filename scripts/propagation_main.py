import utils as utils
from utils import read_prior_set, get_propagation_input
from os import path
from statistics import get_sample_p_values
from propagation_routines import propagate_network, propagate_networks
from args import Args

test_name = 'propagation_main'
args = Args(test_name)

#Read the network
network_graph = utils.read_network(args.network_file)


#Load prior set
prior_set, prior_data, _ = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)

prior_gene_dict = utils.convert_symbols_to_ids(prior_set)
prior_set_ids = list(prior_gene_dict.values())
propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type)


#Use the graph, either run the propagation or load previously acquired propagation results
_, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, propagation_input, args=args, prior_set=prior_set_ids)

genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

# Propagate using randomized networks
_, random_networks_scores = propagate_networks(network_graph, args, list(genes_id_to_idx.keys()), prior_set_ids,
                                               propagation_input, args.random_networks_dir, n_networks=args.n_networks)

# Rank the genes in the original network compared to the random networks
p_values = get_sample_p_values(gene_scores, random_networks_scores)

title = ['gene_id\tp_value\n']
lines = ['{}\t{}\n'.format(genes_idx_to_id[i], p_values[i]) for i in range(len(p_values))]
lines = title + lines

with open(path.join(args.output_folder, 'p_values.txt'), 'w') as f:
    f.writelines(lines)

with open(path.join(args.output_folder, 'significant_p_values.txt'), 'w') as f:
    significant_p_values = dict()
    significant_ids = dict()
    for i in range(len(p_values)):
        if (p_values[i] <= 0.05):
            significant_p_values[i] = p_values[i]
            significant_ids[i] = genes_idx_to_id[i]
    significant = ['{}\t{}\n'.format(significant_ids[k], significant_p_values[k]) for k,v in significant_p_values.items()]
    f.writelines(significant)

lines = ['{}\n'.format(p) for p in prior_set_ids]
with open(path.join(args.output_folder, 'prior_set.txt'), 'w') as f:
    f.writelines(lines)
