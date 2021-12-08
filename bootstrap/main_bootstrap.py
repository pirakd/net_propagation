from bootstrap.methods import fetch_experiment_data, change_expression, shuffle_data, get_statistic_score
import numpy as np
import utils
from args import Args
from statistic_methods import get_sample_p_values
from propagation_routines import generate_similarity_matrix, propagate
from tqdm import tqdm
from random import shuffle
from statistic_methods import wilcoxon_rank_sums_test
import gseapy as gs
import pandas as pd
# arguments
args = Args(is_create_output_folder=False)
network = utils.read_network(args.network_file_path)
network_nodes = network.nodes
pathway_to_change = 'KEGG_TGF_BETA_SIGNALING_PATHWAY'
# pathway_to_change = 'REACTOME_COMPLEMENT_CASCADE'
change_prop = 0.2
change_amount = 0
n_bootstraps = 100
n_randomizations = 1000
raw_data_file = 'raw_data.txt'
true_labels_file = 'labels.txt'

# load data
raw_data, labels, genes_names, symbol_to_id, id_to_experiment_idx = fetch_experiment_data(raw_data_file, true_labels_file, args)
pathway_genes = utils.load_pathways_genes(args.pathway_file_dir,[pathway_to_change])
gene_candidates_to_change = pathway_genes[pathway_to_change]
gene_candidates_idx_to_change = [id_to_experiment_idx[x] for x in gene_candidates_to_change if x in id_to_experiment_idx]
num_of_genes_changed = np.ceil(change_prop * len(gene_candidates_to_change))
genes_idx_to_change = np.random.choice(gene_candidates_idx_to_change, int(num_of_genes_changed), replace=False)

network_graph = utils.read_network(args.network_file_path)
matrix, genes = generate_similarity_matrix(network_graph, None)

# inject expression in part of pathway genes
changed_data = change_expression(raw_data, labels, genes_idx_to_change, change_amount)

num_genes = len(genes)
id_to_network_idx = dict([(gene, index) for (index, gene) in enumerate(genes)])

gene_score = get_statistic_score(changed_data, labels, 's2n')
propagation_input = {id: gene_score[id_to_experiment_idx[id]] for id in id_to_experiment_idx.keys()}
true_scores = propagate([id for id in propagation_input.keys()], propagation_input, matrix, id_to_network_idx, num_genes, args)


tests = ['fc_wilcox', 'prop_wilcox', 'fc_gsea', 'prop_prerank_gsea', 'fc_prerank_gsea']
aggregated_results = np.zeros((n_bootstraps, 5))


for i in tqdm(range(n_bootstraps), desc='bootstrap loop', total=n_bootstraps, position=0, leave=True):
    # bootstrap
    shuffled_data, shuffled_labels = shuffle_data(changed_data, labels)
    bootstrap_gene_score = get_statistic_score(shuffled_data, shuffled_labels, 's2n')
    # propagate
    bootstrap_propagation_input = {id: bootstrap_gene_score[id_to_experiment_idx[id]] for id in id_to_experiment_idx.keys()}
    bootstraped_true_score = propagate([id for id in bootstrap_propagation_input.keys()], bootstrap_propagation_input, matrix, id_to_network_idx, num_genes, args)

    values = list(bootstrap_propagation_input.values())
    random_scores = list()
    for j in tqdm(range(n_randomizations), desc='Propagating randomization for bootstrap {}'.format(i),
                  total=n_randomizations, position=0, leave=True):
        shuffle(values)
        shuffled_input = dict(zip(bootstrap_propagation_input.keys(), values))
        random_scores.append(propagate([id for id in shuffled_input.keys()], shuffled_input, matrix, id_to_network_idx,
                                       num_genes, args))

    random_scores = np.array(random_scores)
    _, _, ranks = get_sample_p_values(bootstraped_true_score, random_scores, two_tailed=True)

    pathway_indexes, background_indexes = [], []
    fc_pathway_indexes, fc_background_indexes = [], []
    for id in bootstrap_propagation_input.keys():
        if id in gene_candidates_to_change:
            pathway_indexes.append(id_to_network_idx[id])
            fc_pathway_indexes.append(id_to_experiment_idx[id])
        else:
            background_indexes.append(id_to_network_idx[id])
            fc_background_indexes.append(id_to_experiment_idx[id])

    pathway_scores, background_scores= ranks[pathway_indexes], ranks[background_indexes]
    result_wilcoxon_propagation = wilcoxon_rank_sums_test(pathway_scores, background_scores).p_value

    fc_pathway_scores, fc_background_scores = bootstrap_gene_score[fc_pathway_indexes], bootstrap_gene_score[fc_background_indexes]
    result_wilcoxon_fc = wilcoxon_rank_sums_test(fc_pathway_scores, fc_background_scores).p_value

    idx_to_experiment_id = {xx:x for x,xx in id_to_experiment_idx.items()}
    data = pd.DataFrame(shuffled_data)
    data.insert(loc=0, column='gene_id', value=[idx_to_experiment_id[i] for i in range(len(idx_to_experiment_id))])
    results_gsea_fc = gs.gsea(data=data, gene_sets={pathway_to_change: gene_candidates_to_change }, cls=shuffled_labels, method='s2n', outdir='test',
                permutation_num=1000,
                no_plot=False, verbose=True, weighted_score_type=1, permutation_type='phenotype').results[pathway_to_change]['pval']

    experiment_genes_scores = []
    for id in bootstrap_propagation_input.keys():
        experiment_genes_scores.append(bootstrap_gene_score[id_to_experiment_idx[id]])

    data = pd.DataFrame(experiment_genes_scores)
    data.insert(loc=0, column='gene_id', value=[id for id in bootstrap_propagation_input.keys()])
    results_fc_prerank_gsea = gs.prerank(rnk=data, gene_sets={pathway_to_change: gene_candidates_to_change }, outdir='test',
                no_plot=False, verbose=True, permutation_num=1000).results[pathway_to_change]['pval']

    experiment_genes_ranks = []
    for id in bootstrap_propagation_input.keys():
        experiment_genes_ranks.append(ranks[id_to_network_idx[id]])

    data = pd.DataFrame(experiment_genes_ranks)
    data.insert(loc=0, column='gene_id', value=[id for id in bootstrap_propagation_input.keys()])
    results_propagation_prerank_gsea = gs.prerank(rnk=data, gene_sets={pathway_to_change: gene_candidates_to_change },
                                                  outdir='test',
                no_plot=False, verbose=True, permutation_num=1000).results[pathway_to_change]['pval']

    aggregated_results[i, :] = [result_wilcoxon_fc, result_wilcoxon_propagation, results_gsea_fc,
                                results_fc_prerank_gsea,  results_propagation_prerank_gsea]

header = '\t'.join(tests)
np.savetxt('change_prop_{}_change_amount_{}'.format(change_prop, change_amount), aggregated_results, fmt='%.2e' ,
           header=header, delimiter='\t')