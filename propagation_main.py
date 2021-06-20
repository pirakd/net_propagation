import utils as utils
from utils import read_prior_set, create_output_folder
from os import path
from propagation_routines import propagate_network, propagate_networks, get_genes_p_values, propagate_networks_parallel
from prior_conditions import get_condition_function


# ~~~ general parameters ~~~
data_dir = 'data'
network_file = 'H_sapiens.net'
experiment_file = 'Table_S1_V1.xlsx'
sheet_name = 'Protein_Abundance'
random_networks_dir = 'random_networks'
condition_function_name = 'kent_mock_no_vic_mock_24h'
n_networks = 2
n_processes = 4

# ~~~ derived parameters ~~~
test_name = '{}_{}'.format(condition_function_name, n_networks)
output_folder = create_output_folder(test_name)
condition_function = get_condition_function(condition_function_name)
# get conditions on experiment
network_file = path.join(data_dir, network_file)
experiment_file_path = path.join(data_dir, experiment_file)

# Read the h_sapiens network
network_graph = utils.read_network(network_file)

# loading prior set
prior_set = read_prior_set(condition_function, experiment_file_path, sheet_name)
prior_gene_dict = utils.convert_symbols_to_ids(prior_set)
prior_set_ids = list(prior_gene_dict.keys())

# Using the graph, either run the propagation or load previously acquired propagation results
_, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, prior_set_ids)
genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

# Propagate using randomized networks
# _, random_networks_scores = propagate_networks(network_graph, prior_set_ids, random_networks_dir, n_networks=n_networks)
random_networks_scores = propagate_networks_parallel(network_graph, prior_set_ids, random_networks_dir, n_networks=n_networks, n_processes=n_processes)

# Rank the genes in the original network compared to the random networks
p_values = get_genes_p_values(gene_scores, random_networks_scores)

title = ['gene_id\tp_value\n']
lines = ['{}\t{}\n'.format(genes_idx_to_id[i], p_values[i]) for i in range(len(p_values))]
lines = title + lines

with open(path.join(output_folder, 'p_values'), 'w') as f:
    f.writelines(lines)

lines = ['{}\n'.format(p) for p in prior_set_ids]
with open(path.join(output_folder, 'prior_set'), 'w') as f:
    f.writelines(lines)
