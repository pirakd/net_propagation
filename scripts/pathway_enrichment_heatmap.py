import sys
from os import path, makedirs
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, load_file,\
    load_propagation_scores
from visualization_tools import plot_enrichment_table
import numpy as np
from args import Args, CovArgs
from statistic_methods import bh_correction, get_stat_test_func, proportion_test
test_name = path.basename(__file__).split('.')[0]
from prior_conditions import huntington_DDA_significant, huntington_DDA, colorectal_cancer_significant,\
    colorectal_cancer, cov_data, cov_data_significant
from collections import namedtuple
Results = namedtuple('Results', 'p_value direction adj_p_value')

test_name = 'ones normalization'
args = CovArgs(test_name)
# colls_names = ['Significant genes proportion test', 'Log2FC signed', 'Propagation signed',
#               'Propagation unsigned']
# colls_names = ['Log2FC signed', 'Propagation signed', 'Propagation unsigned']
colls_names = ['Log2FC signed', 'Propagation signed']

propagation_weight = 0.9
alpha_values_list = [0.9]
min_genes_in_pathway = 5
min_genes_for_display = 3
minimum_pathway_adjusted_p_value_for_display = 0.05
normalization_method = 'ones'
add_self_prop_to_norm_factor = True

interesting_database = ['_']
# interesting_pathways_keywords = ['IMMUNE','cytokine']
interesting_pathways_keywords = ['_']

add_self_propagation = True
include_prior_scores = False
only_excel_genes = True
only_significant = True

all_genes_condition = cov_data
significant_genes_condition = cov_data_significant
signed_input = 'Score'
unsigned_input = 'abs_Score'
n_tests = len(alpha_values_list)
n_colls = len(colls_names)
statistic_test_list = ['wilcoxon_rank_sums_test'] * n_tests
experiment_list = ['Kent_10h-VIC_10h', 'Kent_24h-VIC_24h' ]
# load network
network_graph = utils.read_network(args.network_file_path)
network_genes_ids = set(list(network_graph.nodes()))

pathways_to_display = set()
for i in range(n_tests):
    prop_unsigned_scores_dict = dict()
    prop_signed_scores_dict = dict()
    fc_unsigned_scores_dict = dict()
    fc_signed_scores_dict = dict()
    proportion_scores_dict = dict()
    for e, experiment in enumerate(experiment_list):
        args.sheet_name = experiment
        args.alpha = alpha_values_list[i]
        args.normalization_method = normalization_method
        args.add_self_prop_to_norm_factor = add_self_prop_to_norm_factor

        statistic_test = get_stat_test_func(statistic_test_list[i])
        title = ('Pathway Enrichment, {}'.format(args.sheet_name))
        fig_name = 'enrichment_{:2d}.pdf'.format(np.int(args.alpha*100))

        # read significant genes according to experiment
        # significant_genes, _, _ = read_prior_set(huntington_DDA_significant, args.experiment_file_path, args.sheet_name)
        significant_genes, _, _ = read_prior_set(significant_genes_condition, args.experiment_file_path, args.sheet_name)
        # get all genes list
        # all_genes, all_data, _ = read_prior_set(huntington_DDA, args.experiment_file_path, args.sheet_name)
        all_genes, all_data, _ = read_prior_set(all_genes_condition, args.experiment_file_path, args.sheet_name)
        all_genes_dict = utils.convert_symbols_to_ids(all_genes, args.genes_names_file_path)

        # get excel scores for all genes which appear in the network
        reference_scores = get_propagation_input(all_genes_dict, all_data, 'Score', network=network_graph)
        reference_scores = {x: xx for x, xx in reference_scores.items() if x in network_genes_ids}

        # load propagation scores
        args.propagation_input_type = unsigned_input
        unsigned_gene_scores, _, _, unsigned_propagation_input = load_propagation_scores(args, add_self_propagation=add_self_propagation, normalize_score=True)
        args.propagation_input_type = signed_input
        signed_gene_scores, genes_idx_to_id, genes_id_to_idx, signed_propagation_input =\
            load_propagation_scores(args, add_self_propagation=add_self_propagation, normalize_score=True)

        if only_excel_genes:
            genes_id_to_idx = {id : idx for id, idx in genes_id_to_idx.items() if id in reference_scores}
            genes_idx_to_id = {idx: id for id, idx in genes_id_to_idx.items()}
            if include_prior_scores:
                scores = np.array([signed_gene_scores[idx] for idx in genes_idx_to_id.keys()])
                scores = scores / np.linalg.norm(scores)
                prior_scores = np.array([signed_propagation_input[id] for id in genes_id_to_idx.keys()])
                prior_scores = prior_scores / np.linalg.norm(prior_scores)
                signed_gene_scores[[x for x in genes_id_to_idx.values()]] = propagation_weight * scores + (1 - propagation_weight) * prior_scores

                scores = np.array([unsigned_gene_scores[idx] for idx in genes_idx_to_id.keys()])
                scores = scores / np.linalg.norm(scores)
                prior_scores = np.array([unsigned_propagation_input[id] for id in genes_id_to_idx.keys()])
                prior_scores = prior_scores / np.linalg.norm(prior_scores)
                unsigned_gene_scores[[x for x in genes_id_to_idx.values()]] = propagation_weight * scores + (1 - propagation_weight) * prior_scores


        # load all MSigDB pathways and their associated genes
        intersting_pathways = load_interesting_pathways(args.interesting_pathway_file_dir)
        genes_by_pathway_all = load_pathways_genes(args.pathway_file_dir, intersting_pathways)
        names_to_ids = utils.convert_symbols_to_ids(genes_names_file_path=args.genes_names_file_path)
        ids_to_names = {xx: x for x, xx in names_to_ids.items()}

        # filter pathway according to a defined condition
        genes_by_pathway_interesting = {x: xx for x, xx in genes_by_pathway_all.items() if any(key in x for key in interesting_database)}
        genes_by_pathway_interesting = {x: xx for x, xx in genes_by_pathway_interesting.items() if any(str.upper(key) in x for key in interesting_pathways_keywords)}
        interesting_pathways_names = list(genes_by_pathway_interesting.keys())

        filtered_genes_by_pathway = {pathway: [id for id in genes_by_pathway_interesting[pathway] if id in genes_id_to_idx] for pathway
                                          in interesting_pathways_names}
        reference_filtered_genes_by_pathway = {pathway: [id for id in filtered_genes_by_pathway[pathway] if id in reference_scores] for pathway
                                               in interesting_pathways_names}
        pathways_with_many_genes = [pathway_name for pathway_name in filtered_genes_by_pathway.keys() if
                                    (len(filtered_genes_by_pathway[pathway_name]) >= min_genes_in_pathway)]

        pathways_to_display.update(pathways_with_many_genes)
        # get scores of genes related to interesting pathway
        prop_filtered_signed_genes_scores_by_filtered_pathways = {pathway:[signed_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                           in pathways_with_many_genes}
        prop_filtered_unsigned_genes_scores_by_filtered_pathways = {pathway:[unsigned_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                           in pathways_with_many_genes}

        reference_filtered_genes_scores_by_filtered_pathways = {pathway: [reference_scores[id] for id in reference_filtered_genes_by_pathway[pathway]] for pathway
                                                                in pathways_with_many_genes}
        prop_unsigned_scores_dict[experiment] = dict()
        prop_signed_scores_dict[experiment] = dict()
        fc_signed_scores_dict[experiment] = dict()
        proportion_scores_dict[experiment] = dict()

        # propotion test for significant experiment genes
        network_genes_names = [ids_to_names[id] for id in genes_id_to_idx.keys() if id in ids_to_names]
        all_pathways_genes_by_name = {pathway: [ids_to_names[id] for id in genes if (id in ids_to_names)]
                                     for pathway, genes in genes_by_pathway_all.items()}
        significant_in_net = [gene_name for gene_name in significant_genes if gene_name in network_genes_names]
        proportion_res = proportion_test(significant_in_net, all_pathways_genes_by_name, pathways_with_many_genes,
                                         constrain_to_netowrk_genes=True, network_genes=network_genes_names)

        for p,pathway in enumerate(pathways_with_many_genes):

            proportion_scores_dict[experiment][pathway] = dict()
            proportion_scores_dict[experiment][pathway]['p_value'] = [proportion_res[pathway] for pathway in pathways_with_many_genes]
            proportion_scores_dict[experiment][pathway]['direction'] = [True] * len(pathways_with_many_genes)
            proportion_scores_dict[experiment][pathway]['adj_p_value'] = bh_correction(proportion_scores_dict[experiment][pathway]['p_value'])

            reference_genes = [id for id in reference_scores.keys() if id not in reference_filtered_genes_by_pathway[pathway]]
            reference_background_scores = [reference_scores[id] for id in reference_genes]
            path_scores = reference_filtered_genes_scores_by_filtered_pathways[pathway]
            fc_signed_scores_dict[experiment][pathway] = dict()
            prop_signed_scores_dict[experiment][pathway] = dict()
            prop_unsigned_scores_dict[experiment][pathway] = dict()

            if len(path_scores) >= min_genes_for_display:
            # signed log2FC Wilcoxon rank sum test
                if pathway == 'REACTOME_INFECTIOUS_DISEASE':
                    a=1
                signed_fc_res = statistic_test(path_scores, reference_background_scores, alternative='two-sided')
                fc_signed_scores_dict[experiment][pathway]['p_value'] = signed_fc_res.p_value
                fc_signed_scores_dict[experiment][pathway]['direction'] = signed_fc_res.directionality
                dir = np.mean(path_scores) > np.mean(reference_background_scores)

            else:
                fc_signed_scores_dict[experiment][pathway]['p_value'] = 1.001
                fc_signed_scores_dict[experiment][pathway]['direction'] = True

                dir = True


            prop_random_genes = [id for id in genes_id_to_idx.keys() if id not in filtered_genes_by_pathway[pathway]]
            # signed propagation Wilcoxon rank sum test
            prop_signed_res = statistic_test(prop_filtered_signed_genes_scores_by_filtered_pathways[pathway],
                                      [signed_gene_scores[genes_id_to_idx[id]] for id in prop_random_genes],
                                      alternative='two-sided')
            prop_signed_scores_dict[experiment][pathway]['p_value'] =prop_signed_res.p_value
            # prop_signed_scores_dict.direction.append(prop_signed_res.directionality)
            prop_signed_scores_dict[experiment][pathway]['direction'] = dir

            # unsigned propagation Wilcoxon rank sum test
            prop_unsigned_res = statistic_test(prop_filtered_unsigned_genes_scores_by_filtered_pathways[pathway],
                                      [unsigned_gene_scores[genes_id_to_idx[id]] for id in prop_random_genes],
                                      alternative='greater')
            prop_unsigned_scores_dict[experiment][pathway]['p_value'] = prop_unsigned_res.p_value
            prop_unsigned_scores_dict[experiment][pathway]['direction'] = dir


    # ~~~~~~~
    # ~~~~~~~
    # ~~~~~~
    pathways_to_display = np.sort(list(pathways_to_display))
    p_vals_mat = np.ones((len(pathways_to_display), n_colls * len(experiment_list)))
    adj_p_vals_mat = np.ones_like(p_vals_mat)
    directions_mat = np.zeros_like(p_vals_mat)


    coll_names_in_heatmap = []
    for e, experiment in enumerate(experiment_list):
        # results = [fc_signed_scores_dict[experiment], prop_signed_scores_dict[experiment],
        # prop_unsigned_scores_dict[experiment]]
        results = [fc_signed_scores_dict[experiment], prop_signed_scores_dict[experiment]]
        for i in range(n_colls):
            coll_names_in_heatmap.append('{}, {}'.format(experiment, colls_names[i]))
            index_shift = e * n_colls
            for p, pathway in enumerate(list(pathways_to_display)):
                if pathway in results[i]:
                    p_vals_mat[p, i + index_shift] = results[i][pathway]['p_value']
                    directions_mat[p, i + index_shift] = results[i][pathway]['direction']
                else:
                    a = 1
                    # p_vals_mat[p, i + index_shift] = -1
                    # directions_mat[p, i + index_shift] = 0
            indexes = np.nonzero(p_vals_mat[:, i + index_shift] != -1)[0]
            adj_p_vals_mat[indexes, i + index_shift] = bh_correction(p_vals_mat[indexes, i + index_shift])

    if only_significant:
        keep_rows = np.nonzero(np.any(adj_p_vals_mat <= minimum_pathway_adjusted_p_value_for_display, axis=1))[0]
        n_pathways_before = len(pathways_to_display)
        pathways_to_display = [pathways_to_display[x] for x in keep_rows]

        p_vals_mat = p_vals_mat[keep_rows, :]
        adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
        directions_mat = directions_mat[keep_rows, :]
        # directions_mat[:, 0] = directions_mat[:, 1]

    title =  'Pathway Enrichment:'
    if only_excel_genes:
        # title = '{} only experiment genes'.format(title)
        pass
    if include_prior_scores:
        title = title + ' prior weighting, c={}'.format(propagation_weight)
    else:
        if add_self_propagation:
            title = title + ' with self prop'
        else:
            pass
            # title = title + ' without self prop'

    if only_significant:
        pass
        # title = ' {}  Only significant pathways({}/{}'.format(title, len(keep_rows), n_pathways_before)

    # title = title + '\n propagation from significant genes'

    row_names = ['({}) {}'.format(len(filtered_genes_by_pathway[pathway]), pathway) for pathway in list(pathways_to_display)]
    res = -np.log10(p_vals_mat)
    fig_out_dir = path.join(args.output_folder, fig_name)
    plot_enrichment_table(res, adj_p_vals_mat, directions_mat, row_names, fig_out_dir,
                          experiment_names=coll_names_in_heatmap, title=title, res_type='-log10(p_val)')