import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, load_file,\
    load_propagation_scores
from os import path
from visualization_tools import plot_enrichment_table
import numpy as np
from args import Args
from statistic_methods import bh_correction, get_stat_test_func, proportion_test
test_name = path.basename(__file__).split('.')[0]
from prior_conditions import huntington_DDA_significant, huntington_DDA, colorectal_cancer_significant, colorectal_cancer
from collections import namedtuple
Results = namedtuple('Results', 'p_value direction adj_p_value')
args = Args(test_name)
colls_names = ['Significant genes proportion test', 'Log2FC signed', 'Propagation signed',
              'Propagation unsigned']

alpha_values_list = [0.8, 0.9]
min_genes_in_pathway = 10
min_genes_for_display = 3

# colorectal cancer
# all_genes_condition = colorectal_cancer
# significant_genes_condition = colorectal_cancer_significant

# huntingtons
all_genes_condition = huntington_DDA
significant_genes_condition = huntington_DDA_significant

n_tests = len(alpha_values_list)
n_colls = len(colls_names)

# Huntington
interesting_pathways_keywords = ['KEGG']
sheet_names_list = ['Table_A'] * n_tests
statistic_test_list = ['wilcoxon_rank_sums_test'] * n_tests

# load network
network_graph = utils.read_network(args.network_file)
network_genes_ids = set(list(network_graph.nodes()))

for i in range(n_tests):
    args.alpha = alpha_values_list[i]
    args.sheet_name = sheet_names_list[i]
    statistic_test = get_stat_test_func(statistic_test_list[i])
    title = ('Pathway Enrichment')
    fig_name = 'enrichment_{:2d}.pdf'.format(np.int(args.alpha*100))

    # read significant genes according to experiment
    # significant_genes, _, _ = read_prior_set(huntington_DDA_significant, args.experiment_file_path, args.sheet_name)
    significant_genes, _, _ = read_prior_set(significant_genes_condition, args.experiment_file_path, args.sheet_name)
    # get all genes list
    # all_genes, all_data, _ = read_prior_set(huntington_DDA, args.experiment_file_path, args.sheet_name)
    all_genes, all_data, _ = read_prior_set(all_genes_condition, args.experiment_file_path, args.sheet_name)
    all_genes = list(all_data.Gene_Name)
    all_genes_dict = utils.convert_symbols_to_ids(all_genes, args.genes_names_file_path)

    # get excel scores for all genes which appear in the network
    reference_scores = get_propagation_input(all_genes_dict, all_data, 'Score')
    reference_scores = {x: xx for x, xx in reference_scores.items() if x in network_genes_ids}

    # load propagation scores
    args.propagation_input_type = 'abs_score_all'
    unsigned_gene_scores, _, _ = load_propagation_scores(args, normalize_score=True)
    args.propagation_input_type = 'Score'
    signed_gene_scores, genes_idx_to_id, genes_id_to_idx = load_propagation_scores(args, normalize_score=True)

    # load all MSigDB pathways and their associated genes
    genes_by_pathway_all = load_pathways_genes(args.pathway_file_dir)
    names_to_ids = utils.convert_symbols_to_ids(genes_names_file_path=args.genes_names_file_path)
    ids_to_names = {xx: x for x, xx in names_to_ids.items()}
    all_pathways_genes_by_name = {pathway: [ids_to_names[id] for id in genes if id in ids_to_names]
                                 for pathway, genes in genes_by_pathway_all.items()}

    # filter pathway according to a defined condition
    genes_by_pathway_interesting = {x: xx for x, xx in genes_by_pathway_all.items() if any(key in x for key in interesting_pathways_keywords)}
    interesting_pathways_names = list(genes_by_pathway_interesting.keys())

    filtered_genes_by_pathway = {pathway: [id for id in genes_by_pathway_interesting[pathway] if id in genes_id_to_idx] for pathway
                                      in interesting_pathways_names}
    reference_filtered_genes_by_pathway = {pathway: [id for id in filtered_genes_by_pathway[pathway] if id in reference_scores] for pathway
                                           in interesting_pathways_names}
    pathways_with_many_genes = [pathway_name for pathway_name in filtered_genes_by_pathway.keys() if
                                (len(filtered_genes_by_pathway[pathway_name]) >= min_genes_in_pathway)]


    # get scores of genes related to interesting pathway
    prop_filtered_signed_genes_scores_by_filtered_pathways = {pathway:[signed_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                       in pathways_with_many_genes}
    prop_filtered_unsigned_genes_scores_by_filtered_pathways = {pathway:[unsigned_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                       in pathways_with_many_genes}

    reference_filtered_genes_scores_by_filtered_pathways = {pathway: [reference_scores[id] for id in reference_filtered_genes_by_pathway[pathway]] for pathway
                                                            in pathways_with_many_genes}
    prop_unsigned_scores_dict = Results([], [], [])
    prop_signed_scores_dict = Results([], [], [])
    fc_unsigned_scores_dict = Results([], [], [])
    fc_signed_scores_dict = Results([], [], [])
    proportion_scores_dict = Results([],[],[])

    # propotion test for significant experiment genes
    proportion_res = proportion_test(significant_genes, all_pathways_genes_by_name, pathways_with_many_genes)
    proportion_scores_dict.p_value.extend([proportion_res[pathway] for pathway in pathways_with_many_genes])
    proportion_scores_dict.direction.extend([True] * len(pathways_with_many_genes))
    proportion_scores_dict.adj_p_value.extend(bh_correction(proportion_scores_dict.p_value))

    for p,pathway in enumerate(pathways_with_many_genes):

        reference_genes = [id for id in reference_scores.keys() if id not in reference_filtered_genes_by_pathway[pathway]]
        reference_background_scores = [reference_scores[id] for id in reference_genes]
        path_scores = reference_filtered_genes_scores_by_filtered_pathways[pathway]

        if len(path_scores) >= min_genes_for_display:
        # signed log2FC Wilcoxon rank sum test
            signed_fc_res = statistic_test(path_scores, reference_background_scores, alternative='two-sided')
            fc_signed_scores_dict.p_value.append(signed_fc_res.p_value)
            fc_signed_scores_dict.direction.append(signed_fc_res.directionality)

            dir = np.mean(path_scores) > np.mean(reference_background_scores)

        else:
            fc_signed_scores_dict.p_value.append(1.001)
            fc_signed_scores_dict.direction.append(True)
            fc_unsigned_scores_dict.p_value.append(1.001)
            fc_unsigned_scores_dict.direction.append(True)
            dir = True


        prop_random_genes = [id for id in genes_id_to_idx.keys() if id not in filtered_genes_by_pathway[pathway]]
        # signed propagation Wilcoxon rank sum test
        prop_signed_res = statistic_test(prop_filtered_signed_genes_scores_by_filtered_pathways[pathway],
                                  [signed_gene_scores[genes_id_to_idx[id]] for id in prop_random_genes],
                                  alternative='two-sided')
        prop_signed_scores_dict.p_value.append(prop_signed_res.p_value)
        prop_signed_scores_dict.direction.append(prop_signed_res.directionality)

        # unsigned propagation Wilcoxon rank sum test
        prop_unsigned_res = statistic_test(prop_filtered_unsigned_genes_scores_by_filtered_pathways[pathway],
                                  [unsigned_gene_scores[genes_id_to_idx[id]] for id in prop_random_genes],
                                  alternative='greater')
        prop_unsigned_scores_dict.p_value.append(prop_unsigned_res.p_value)
        prop_unsigned_scores_dict.direction.append(dir)

    p_vals_mat = np.ones((len(pathways_with_many_genes), n_colls))
    adj_p_vals_mat = np.ones_like(p_vals_mat)
    directions_mat = np.zeros_like(p_vals_mat)

    results = [proportion_scores_dict, fc_signed_scores_dict, prop_signed_scores_dict,
               prop_unsigned_scores_dict]

    for i in range(n_colls):
        p_vals_mat[:, i] = results[i].p_value
        adj_p_vals_mat[:, i] = bh_correction(results[i].p_value)
        directions_mat[:, i] = results[i].direction

    keep_rows = np.nonzero(np.any(adj_p_vals_mat<0.01 , axis=1 ))[0]

    pathways_with_many_genes = [pathways_with_many_genes[x] for x in keep_rows]
    p_vals_mat = p_vals_mat[keep_rows, :]
    adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
    directions_mat = directions_mat[keep_rows, :]
    res = -np.log10(p_vals_mat)
    fig_out_dir = path.join(args.output_folder, fig_name)
    plot_enrichment_table(res ,adj_p_vals_mat, directions_mat, pathways_with_many_genes, fig_out_dir,
                          experiment_names=colls_names, title=title, res_type='-log10(p_val)')