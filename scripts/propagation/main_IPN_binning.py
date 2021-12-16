import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
from statistic_methods import get_sample_p_values
import utils as utils
from utils import get_propagation_input, save_propagation_score, generate_bins
from propagation_routines import generate_similarity_matrix, propagate
from args import DeltaArgs
import numpy as np
from tqdm import tqdm

test_name = path.basename(__file__).split('.')[0]
args = DeltaArgs(test_name, is_create_output_folder=False)
n_randomizations = 1000
precision = 6
min_bin_size = 30
args.alpha = 0.9

args.get_derived_parameters()
network_graph = utils.read_network(args.network_file_path)
# network_graph = remove_outlier_vertices(network_graph, 800)
prior_set, prior_data, _ = args.experiment_reader(args)
prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
prior_set_ids = set.intersection(set(prior_gene_dict.values()), set(network_graph.nodes))

all_genes_ids_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
all_genes_ids = set.intersection(set(all_genes_ids_dict.values()), set(network_graph.nodes))

propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type,
                                          network=network_graph)
ones_input = get_propagation_input(prior_gene_dict, prior_data, 'ones',
                                          network=network_graph)

matrix, genes = generate_similarity_matrix(network_graph, None)
num_genes = len(genes)

genes_id_to_idx = dict([(gene, index) for (index, gene) in enumerate(genes)])
genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}
gene_scores = list()
self_prop_scores = np.zeros(num_genes)

true_score = propagate([id for id in propagation_input.keys()], propagation_input, matrix, genes_id_to_idx, num_genes, args)

ones_centrality = propagate([id for id in propagation_input.keys()], ones_input, matrix, genes_id_to_idx, num_genes, args)
ones_centrality_experiment = np.array([ones_centrality[genes_id_to_idx[id]] for id in propagation_input])
ones_centrality_quant = np.ceil(ones_centrality_experiment*(10**(precision))).astype(int)
genes_idx_to_id_experiment = {idx:id for idx, id in enumerate(propagation_input.keys())}
weight_dict = {genes_idx_to_id_experiment[x]: xx for x, xx in enumerate(ones_centrality_experiment)}
bins, id_to_bin = generate_bins(weight_dict, precision, min_bin_size)

random_scores = list()
for i in tqdm(range(n_randomizations), desc='Propagating random scores', total=n_randomizations):
    shuffled_input = {}
    for id in propagation_input.keys():
            shuffled_input[id] = propagation_input[np.random.choice(bins[id_to_bin[id]], 1)[0]]

    random_scores.append(propagate([id for id in propagation_input.keys()], shuffled_input, matrix, genes_id_to_idx, num_genes
                                   , args))


random_scores = np.array(random_scores)
p_values, _, ranks = get_sample_p_values(true_score, random_scores, two_tailed=True)

ranks_experiment = [ranks[genes_id_to_idx[id]] for id in propagation_input.keys()]

file_name = '{}_{}_{}_{}_{}_IPN_bin_{}'.format(args.propagation_input_type, args.experiment_name, args.sheet_name,
                                        args.experiment_reader_name, str(args.alpha), min_bin_size)

save_propagation_score(propagation_scores=true_score, prior_set=prior_set, propagation_input=propagation_input,
                       genes_idx_to_id=genes_idx_to_id, args=args, self_propagation=self_prop_scores,
                       randomization_ranks=ranks, n_randomizations=n_randomizations, scores_p_values=p_values,
                       file_name=file_name)
