from args import Args, MockArgs
from os import path
from statistic_methods import wilcoxon_rank_sums_test, bh_correction
from utils import read_network, load_interesting_pathways, load_pathways_genes, get_propagation_input,\
    convert_symbols_to_ids, load_propagation_scores, create_output_folder, get_root_path
import numpy as np
from visualization_tools import plot_enrichment_table
from enrichment_utils import GeneralArgs, RawScoreTask, PropagationTask, PathwayResults, filter_pathways

def run(task_list, general_args):
    # load network
    network_graph = read_network(general_args.network_file_path)
    network_genes_ids = set(list(network_graph.nodes()))

    # load interesting pathways and pathway members
    if general_args.interesting_pathway_path is not None:
        interesting_pathways = load_interesting_pathways(general_args.interesting_pathway_path)
    else:
        interesting_pathways = None

    genes_by_pathway = load_pathways_genes(general_args.pathway_members_path, interesting_pathways)
    if interesting_pathways is None:
        interesting_pathways = list(genes_by_pathway.keys())

    interesting_pathways = {x for x in interesting_pathways if
                                    any(key in x for key in general_args.pathway_databases)}
    interesting_pathways = {x for x in interesting_pathways if
                                    any(str.upper(key) in x for key in general_args.pathway_keywords)}
    n_genes_per_task_per_pathway = list()
    pathways_to_display = set()
    for task in task_list:

        # load scores
        if isinstance(task, RawScoreTask):
            all_genes, all_data, _ = task.experiment_reader(task)
            all_genes_dict = convert_symbols_to_ids(all_genes, general_args.genes_names_file_path)
            scores = get_propagation_input(all_genes_dict, all_data, task.propgation_input_type, network=network_graph)

        elif isinstance(task, PropagationTask):
            res_dict = load_propagation_scores(task, add_self_propagation=task.add_self_propagation,
                                               normalization_file_name=task.normalization_file,
                                               normalize_score=task.normalize_scores,
                                               propagation_file_name=task.propagation_file)
            gene_id_to_idx = {xx: x for x, xx in res_dict['gene_idx_to_id'].items()}
            scores = res_dict[task.target_field]
            if task.constrain_to_experiment_genes:
                gene_id_to_idx = {id: idx for id, idx in gene_id_to_idx.items() if id in res_dict['propagation_input']}
            scores = {id:scores[idx] for id, idx in gene_id_to_idx.items()}
        else:
            assert 0, 'invalid task'
        # filter genes for each pathway
        genes_by_pathway_filtered = {pathway: [id for id in genes_by_pathway[pathway] if id in scores] for pathway
                                          in interesting_pathways}
        # keep only pathway with certain amount of genes
        pathways_with_many_genes = [pathway_name for pathway_name in genes_by_pathway_filtered.keys() if
                                    (len(genes_by_pathway_filtered[pathway_name]) >= general_args.minimum_genes_per_pathway)]

        # statistic test
        for p, pathway in enumerate(pathways_with_many_genes):
            pathways_to_display.add(pathway)
            pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
            background_scores = [scores[id] for id in scores.keys() if id not in genes_by_pathway_filtered[pathway]]
            result = task.statistic_test(pathway_scores, background_scores)
            task.results[pathway] = PathwayResults()
            task.results[pathway].p_value = result.p_value
            task.results[pathway].direction = result.directionality
        n_genes_per_task_per_pathway.append({pathway: len(genes_by_pathway_filtered[pathway]) for pathway in genes_by_pathway_filtered})

    pathways_to_display = np.sort(list(pathways_to_display))
    p_vals_mat = np.ones((len(pathways_to_display), len(task_list)))
    adj_p_vals_mat = np.ones_like(p_vals_mat)
    directions_mat = np.zeros_like(p_vals_mat)

    coll_names_in_heatmap = []
    for t, task in enumerate(task_list):
        indexes = []
        for p, pathway in enumerate(pathways_to_display):
            if pathway in task.results:
                indexes.append(p)
                p_vals_mat[p, t] = task.results[pathway].p_value
                directions_mat[p, t] = task.results[pathway].direction
        adj_p_vals_mat[indexes, t] = bh_correction(p_vals_mat[indexes, t])
        coll_names_in_heatmap.append(task.name)

    if general_args.display_only_significant_pathways:
        keep_rows = np.nonzero(np.any(adj_p_vals_mat <= general_args.significant_pathway_threshold, axis=1))[0]
        pathways_to_display = [pathways_to_display[x] for x in keep_rows]

        p_vals_mat = p_vals_mat[keep_rows, :]
        adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
        directions_mat = directions_mat[keep_rows, :]

    if p_vals_mat.shape[0] > general_args.maximum_number_of_pathways and general_args.merge_similar_pathways:
        keep_rows = \
            np.sort(filter_pathways(p_vals_mat, adj_p_vals_mat, pathways_to_display, n_genes_per_task_per_pathway, general_args))


    if p_vals_mat.shape[0] > general_args.maximum_number_of_pathways:
        candidates = np.min(p_vals_mat, axis=1)
        ind = np.sort(np.argpartition(candidates, general_args.maximum_number_of_pathways)[:general_args.maximum_number_of_pathways])
        p_vals_mat, adj_p_vals_mat, directions_mat = p_vals_mat[ind, :], adj_p_vals_mat[ind, :], \
                                                     directions_mat[ind, :]
        pathways_to_display = [pathways_to_display[x] for x in ind]



    row_names = ['({}) {}'.format(len(genes_by_pathway_filtered[pathway]), pathway) for pathway in pathways_to_display]

    res = -np.log10(p_vals_mat)
    fig_out_dir = path.join(general_args.output_path, general_args.figure_name)
    n_pathways_before = len(pathways_with_many_genes)
    plot_enrichment_table(res, adj_p_vals_mat, directions_mat, row_names, fig_out_dir,
                          experiment_names=coll_names_in_heatmap,
                          title=general_args.figure_title + '{}/{}'.format(len(row_names), n_pathways_before),
                          res_type='-log10(p_val)')


if __name__ == '__main__':
    args = MockArgs(is_create_output_folder=False)
    propagation_scores_file = 'Score_mock_scores_Table_A_cov_data_0.9'
    normalization_score_file = 'Score_mock_scores_Table_A_cov_data_0.9'

    task_1 = RawScoreTask(name='Raw Enrichment', experiment_file_path=args.experiment_file_path,
                          sheet_name='Table_A', statistic_test=wilcoxon_rank_sums_test,
                          experiment_reader=args.experiment_reader, propagation_input_type='Score')

    task_2 = PropagationTask(name='Propagation Enrichment', propagation_file=propagation_scores_file,
                             propagation_folder='propagation_scores',
                             normalization_file=normalization_score_file, statistic_test=wilcoxon_rank_sums_test,
                             target_field='gene_prop_score', normalization_method='EC', constrain_to_experiment_genes=True)

    general_args = GeneralArgs(args.network_file_path, genes_names_path=args.genes_names_file_path,
                               interesting_pathway_path=args.interesting_pathway_file_dir,
                               pathway_members_path=args.pathway_file_dir)

    task_list = [task_1, task_2]
    run(task_list, general_args)