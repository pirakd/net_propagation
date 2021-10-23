from args import Args, MockArgs
from os import path
from statistic_methods import wilcoxon_rank_sums_test, bh_correction
from utils import read_network, load_interesting_pathways, load_pathways_genes, get_propagation_input, read_prior_set,\
    convert_symbols_to_ids, load_propagation_scores, create_output_folder
import numpy as np
from visualization_tools import plot_enrichment_table


class PropagationTask:
    def __init__(self, name, propagation_file, normalization_file, statistic_test, target_field,
                 normalization_method, constrain_to_experiment_genes, add_self_propagation=False,
                 normalize_scores=True):
        self.name = name
        self.propagation_file = propagation_file
        self.normalization_file = normalization_file
        self.statistic_test = statistic_test
        self.target_field = target_field
        self.normalize_scores = normalize_scores
        self.normalization_method = normalization_method
        self.constrain_to_experiment_genes = constrain_to_experiment_genes
        self.add_self_propagation = add_self_propagation
        self.add_self_propagation_to_norm_factor = add_self_propagation
        self.results = dict()


class RawScoreTask:
    def __init__(self, name, score_file_path, statistic_test, condition_function, propagation_input_type,
                 constrain_to_network_genes ):
        self.name = name
        self.score_file_path = score_file_path
        self.statistic_test = statistic_test
        self.condition_function = condition_function
        self.propgation_input_type = propagation_input_type
        self.constrain_to_experiment = constrain_to_network_genes
        self.results = dict()


class GeneralArgs:
    def __init__(self, network_path, genes_names_path, interesting_pathway_path, pathway_members_path,
                 output_folder_name=None, figure_name=None, figure_title='Pathway Enrichment'):

        self.minimum_gene_per_pathway = 1
        self.display_only_significant_pathways = False
        self.network_file_path = network_path
        self.genes_names_file_path = genes_names_path
        self.pathway_databases = ['_']
        self.pathway_keywords = ['_']
        self.significant_pathway_threshold = 5e-2
        if output_folder_name is None:
            output_folder_name = path.basename(__file__).split('.')[0]
        self.output_path, _ = create_output_folder(output_folder_name)
        self.figure_name = figure_name if figure_name is not None else 'figure'
        self.interesting_pathway_path = interesting_pathway_path
        self.pathway_members_path = pathway_members_path
        self.figure_title = figure_title


class PathwayResults:
    def __init__(self):
        self.p_value = None
        self.direction = None
        self.adj_p_value = None


def run(task_list, general_args):

    # load network
    network_graph = read_network(general_args.network_file_path)
    network_genes_ids = set(list(network_graph.nodes()))

    # load interesting pathways and pathway members
    intersting_pathways = load_interesting_pathways(general_args.interesting_pathway_path)
    intersting_pathways = {x for x in intersting_pathways if
                                    any(key in x for key in general_args.pathway_databases)}
    intersting_pathways = {x for x in intersting_pathways if
                                    any(str.upper(key) in x for key in general_args.pathway_keywords)}
    genes_by_pathway = load_pathways_genes(general_args.pathway_members_path, intersting_pathways)
    names_to_ids =  convert_symbols_to_ids(genes_names_file_path=general_args.genes_names_file_path)
    ids_to_names = {xx: x for x, xx in names_to_ids.items()}
    pathways_to_display = set()

    for task in task_list:

        # load scores
        if isinstance(task, RawScoreTask):
            all_genes, all_data, _ = read_prior_set(task.condition_function, args.experiment_file_path, args.sheet_name)
            all_genes_dict = convert_symbols_to_ids(all_genes, general_args.genes_names_file_path)
            scores = get_propagation_input(all_genes_dict, all_data, task.propgation_input_type, network=network_graph)

        elif isinstance(task, PropagationTask):
            args = Args(is_create_output_folder=False)
            args.normalization_method = task.normalization_method
            args.add_self_prop_to_norm_factor = task.add_self_propagation_to_norm_factor
            scores, genes_idx_to_id, genes_id_to_idx, propagation_input = load_propagation_scores(args,
                                                                      add_self_propagation=task.add_self_propagation,
                                                                      normalization_file_name=task.normalization_file,
                                                                      normalize_score=task.normalize_scores,
                                                                      propagation_file_name=task.propagation_file)
            if task.constrain_to_experiment_genes:
                genes_id_to_idx = {id: idx for id, idx in genes_id_to_idx.items() if id in propagation_input}
            scores = {id:scores[idx] for id, idx in genes_id_to_idx.items()}
        else:
            assert 0, 'invalid task'
        # filter genes for each pathway
        genes_by_pathway_filtered = {pathway: [id for id in genes_by_pathway[pathway] if id in scores] for pathway
                                          in intersting_pathways}
        # keep only pathway with certain amount of genes
        pathways_with_many_genes = [pathway_name for pathway_name in genes_by_pathway_filtered.keys() if
                                    (len(genes_by_pathway_filtered[pathway_name]) >= general_args.minimum_gene_per_pathway)]

        # statistic test
        for p, pathway in enumerate(pathways_with_many_genes):
            pathways_to_display.add(pathway)
            pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
            background_scores = [scores[id] for id in scores.keys() if id not in genes_by_pathway_filtered[pathway]]
            result = task.statistic_test(pathway_scores, background_scores)
            task.results[pathway] = PathwayResults()
            task.results[pathway].p_value = result.p_value
            task.results[pathway].direction = result.directionality

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
        n_pathways_before = len(pathways_to_display)
        pathways_to_display = [pathways_to_display[x] for x in keep_rows]

        p_vals_mat = p_vals_mat[keep_rows, :]
        adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
        directions_mat = directions_mat[keep_rows, :]

    row_names = ['({}) {}'.format(len(genes_by_pathway_filtered[pathway]), pathway) for pathway in pathways_to_display]
    res = -np.log10(p_vals_mat)
    fig_out_dir = path.join(general_args.output_path, general_args.figure_name)
    plot_enrichment_table(res, adj_p_vals_mat, directions_mat, row_names, fig_out_dir,
                          experiment_names=coll_names_in_heatmap, title=general_args.figure_title, res_type='-log10(p_val)')


if __name__ == '__main__':
    args = MockArgs(is_create_output_folder=False)
    propagation_scores_file = 'ones_mock_scores_Table_A_cov_data_0.9'
    normalization_score_file = 'ones_mock_scores_Table_A_cov_data_0.9'
    task_1 = PropagationTask(name = 'first', propagation_file=propagation_scores_file,
                             normalization_file=normalization_score_file, statistic_test=wilcoxon_rank_sums_test,
                             target_field='gene_score', normalization_method='EC', constrain_to_experiment_genes=True)

    general_args = GeneralArgs(args.network_file_path, genes_names_path=args.genes_names_file_path,
                               interesting_pathway_path=args.interesting_pathway_file_dir,
                               pathway_members_path=args.pathway_file_dir)

    task_list = [task_1]
    run(task_list, general_args)