import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
from statistic_methods import get_sample_p_values
import utils as utils
from utils import get_propagation_input, save_propagation_score
from propagation_routines import propagate_network, generate_similarity_matrix, propagate
from args import MockArgs, DeltaArgs
import numpy as np


test_name = path.basename(__file__).split('.')[0]
args = DeltaArgs(test_name, is_create_output_folder=False)
n_randomizations = 1000

network_graph = utils.read_network(args.network_file_path)
fc_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list())

prior_set, prior_data, _ = args.experiment_reader(args)
prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
prior_set_ids = set.intersection(set(prior_gene_dict.values()), set(network_graph.nodes))

all_genes = list(prior_data.Gene_Name)
all_genes_ids_dict = utils.convert_symbols_to_ids(all_genes, args.genes_names_file_path)
all_genes_ids = set.intersection(set(all_genes_ids_dict.values()), set(network_graph.nodes))

propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type,
                                          network=network_graph)
ones_input = get_propagation_input(prior_gene_dict, prior_data, 'ones', network=network_graph)

matrix, genes = generate_similarity_matrix(network_graph, None)
num_genes = len(genes)

gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])
genes_idx_to_id = {xx: x for x, xx in gene_indexes.items()}
gene_scores = list()
self_prop_scores = np.zeros(num_genes)
for idx, id in enumerate(list(ones_input.keys())):
    # if id in significant_genes_ids:
    gene_scores.append(propagate([id], ones_input, matrix, gene_indexes, num_genes, args))
    self_prop_scores[gene_indexes[id]] = gene_scores[-1][gene_indexes[id]]

one_scores = np.array(gene_scores)
inputs = np.array([val for val in propagation_input.values()])[:, np.newaxis]
true_score = np.sum(one_scores * inputs, axis=0)

random_scores = []
for i in range(n_randomizations):
    if np.mod(i, 100) == 0:
        print('finished {} randomization'.format(i))
    random_permutation = np.random.permutation(inputs)
    random_scores.append(np.sum(one_scores * random_permutation, axis=0))

random_scores = np.array(random_scores)
p_values, _, ranks = get_sample_p_values(true_score, random_scores, two_tailed=True)

n_experiments = random_scores.shape[0]
sorted_scores = np.sort(random_scores, axis=0)

file_name = '{}_{}_{}_{}_{}_IPN'.format(args.propagation_input_type, args.experiment_name, args.sheet_name,
                                        args.experiment_reader_name, str(args.alpha))

save_propagation_score(propagation_scores=gene_scores, prior_set=prior_set, propagation_input=propagation_input,
                       genes_idx_to_id=genes_idx_to_id, args=args, self_propagation=self_prop_scores,
                       randomization_ranks=ranks, n_randomizations=n_randomizations, scores_p_values=p_values,
                       file_name=file_name)