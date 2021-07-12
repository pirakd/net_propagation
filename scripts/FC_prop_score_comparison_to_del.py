import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, load_file
from os import path
from propagation_routines import propagate_network
from visualization_tools import plot_enrichment_table
import numpy as np
from args import Args
import pickle as pl
from statistics import bh_correction, get_stat_test_func

test_name = path.basename(__file__).split('.')[0]
n_tests = 2
normalize_by_eig_vec_cent_list = [False] * 2
normalize_by_degree_list = [False] * 2
res_type_list = ['-log10(p)'] * 2
log_of_prop_scores = True
sheet_names_list = ['Protein_Abundance', 'Protein_Abundance', 'RNA', 'RNA']
sheet_names_list = ['Suppl. Table 4A'] * 2
statistic_test_list = ['man_whit_U_test'] * 2
prior_set_conditions = ['huntington_DIA', 'huntington_DDA']
reference_score_type_list = ['AVG Log2 Ratio','Log2FC (HD/C116)']
args = Args(test_name)
load_prop_results = True
keep_only_excel_genes = False


for i in range(n_tests):
    args.sheet_name = sheet_names_list[i]
    res_type = res_type_list[i]
    statistic_test = get_stat_test_func(statistic_test_list[i])
    normalize_by_eig_vec_cent = normalize_by_eig_vec_cent_list[i]
    normalize_by_degree = normalize_by_degree_list[i]
    title = '{} - Enrichment of Pathway Genes'.format(args.sheet_name)
    if normalize_by_eig_vec_cent:
        title = title + ' Eigvec normalized'
    fig_name = 'enrichement_{}'.format(i)
    fc_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list(), z_score=[])
    prop_scores_dict = dict(p_vals=list(), adj_p_vals=list(), direction=list(), z_score=[])
    pathways_per_condition = []

    for c, condition in enumerate(prior_set_conditions):
        reference_score_type = reference_score_type_list[c]
        args.condition_function_name = prior_set_conditions[c]
        args.set_condition_function()

        # loading prior set
        prior_set, prior_data, all_data = read_prior_set(args.condition_function, args.experiment_file_path, args.sheet_name)
        prior_gene_dict = utils.convert_symbols_to_ids(prior_set, args.genes_names_file_path)
        prior_set_ids = list(prior_gene_dict.values())
        propagation_input = get_propagation_input(prior_gene_dict, prior_data, args.propagation_input_type)


        all_genes = list(all_data.Gene)
        all_genes_dict = utils.convert_symbols_to_ids(all_genes, args.genes_names_file_path)
        all_reference_scores = get_propagation_input(all_genes_dict, all_data, reference_score_type)

        if load_prop_results:
            propagation_file_name = '{}_{}_{}'.format(args.sheet_name, condition, str(args.alpha))
            propagation_results_path = path.join(args.propagation_scores_path, propagation_file_name)
            propagation_res_dict = load_file(propagation_results_path, decompress=True)

            gene_scores = np.array(propagation_res_dict['gene_prop_scores'])
            genes_idx_to_id = propagation_res_dict['gene_idx_to_id']
            genes_id_to_idx = {xx: x for x, xx in genes_idx_to_id.items()}
        else:
            # Read the h_sapiens network
            network_graph = utils.read_network(args.network_file)
            # Using the graph, either run the propagation or load previously acquired propagation results
            _, _, genes_id_to_idx, gene_scores = propagate_network(network_graph, propagation_input, args,
                                                                   prior_set=list(prior_gene_dict.values()))
            genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}

        if normalize_by_eig_vec_cent:
            propagation_norm_file_name = '{}_{}_1'.format(args.sheet_name, condition)
            propagation_norm_res_path = path.join(args.propagation_scores_path, propagation_norm_file_name)
            norm_propagation_res_dict = load_file(propagation_norm_res_path, decompress=True)

            norm_genes_idx_to_id = norm_propagation_res_dict['gene_idx_to_id']
            norm_scores = np.array(norm_propagation_res_dict['gene_prop_scores'])
            zero_normalization_genes = np.nonzero(norm_scores == 0)[0]
            zero_prop_genes = np.nonzero(gene_scores == 0)[0]
            genes_to_delete = list(set(zero_normalization_genes).difference(zero_prop_genes))
            norm_scores[genes_to_delete] = 1
            gene_scores[gene_scores != 0] = np.array(gene_scores[gene_scores!=0] / norm_scores[gene_scores!=0])

            for gene_idx in genes_to_delete:
                gene_id = genes_idx_to_id[gene_idx]
                genes_idx_to_id.pop(gene_idx)
                genes_id_to_idx.pop(gene_id)

        if log_of_prop_scores:
            gene_scores[gene_scores != 0] = np.log(gene_scores[gene_scores != 0])
        if keep_only_excel_genes:
            genes_idx_to_id = {idx: id for idx, id in genes_idx_to_id.items() if
                                            id in propagation_input}
            genes_id_to_idx = {xx: x for x, xx in genes_idx_to_id.items()}

        # load pathways and their related genes
        interesting_pathways = load_interesting_pathways(args.interesting_pathway_file_dir)
        interesting_pathways = np.sort(interesting_pathways)
        genes_of_interesting_pathways = load_pathways_genes(args.pathway_file_dir, interesting_pathways)
        # genes_of_interesting_pathways = load_pathways_genes(args.pathway_file_dir)
        interesting_pathways = [x for x,xx in genes_of_interesting_pathways.items()]
        genes_by_pathway_filtered = [[id for id in genes_of_interesting_pathways[pathway] if id in genes_id_to_idx] for pathway
                                   in interesting_pathways]
        excel_genes_by_pathway_filtered = [[id for id in genes_of_interesting_pathways[pathway] if id in all_reference_scores] for pathway
                                   in interesting_pathways]

        pathways_with_many_genes = [x for x, xx in enumerate(genes_by_pathway_filtered) if
                                    (len(genes_by_pathway_filtered[x]) >= 5 and len(excel_genes_by_pathway_filtered[x]) >= 5)]
        genes_by_pathway_filtered = [genes_by_pathway_filtered[x] for x in pathways_with_many_genes]

        # check scores of genes related to interesting pathway
        genes_scores_by_pathway = [[gene_scores[genes_id_to_idx[id]] for id in gene_set] for gene_set
                                   in genes_by_pathway_filtered]
        all_reference_scores_by_pathway = [[all_reference_scores[id] for id in gene_set if id in all_reference_scores] for gene_set
                                    in genes_by_pathway_filtered]

        fc_p_vals, fc_dir, fc_z_score, prop_p_vals, prop_dir, prop_z_score = [], [], [], [], [], []
        for g, gene_set in enumerate(genes_by_pathway_filtered):
            random_genes = [id for id in all_reference_scores if id not in gene_set]

            fc_res = statistic_test(all_reference_scores_by_pathway[g], [all_reference_scores[id] for id in random_genes])
            fc_p_vals.append(fc_res.p_value)
            fc_z_score.append(fc_res.z_score)
            fc_dir.append(fc_res.directionality)

            random_genes = [id for id in genes_id_to_idx.keys() if id not in gene_set]

            prop_res = statistic_test(genes_scores_by_pathway[g], [gene_scores[genes_id_to_idx[id]]
                                                                     for id in random_genes])
            prop_p_vals.append(prop_res.p_value)
            prop_z_score.append(prop_res.z_score)
            prop_dir.append(prop_res.directionality)


        prop_scores_dict['direction'].append(np.array(prop_dir))
        prop_scores_dict['p_vals'].append(np.array(prop_p_vals))
        prop_scores_dict['z_score'].append(np.array(prop_z_score))
        prop_scores_dict['adj_p_vals'].append(np.array(bh_correction(prop_p_vals)))


        fc_scores_dict['direction'].append(np.array(fc_dir))
        fc_scores_dict['p_vals'].append(np.array(fc_p_vals))
        fc_scores_dict['z_score'].append(np.array(fc_z_score))
        fc_scores_dict['adj_p_vals'].append(np.array(bh_correction(fc_p_vals)))

        pathways_per_condition.append(pathways_with_many_genes)


    ## construct the enrichment table
    all_pathways = np.sort(list(set().union(*pathways_per_condition)))
    all_pathways_names = [interesting_pathways[x] for x in all_pathways]
    pathway_to_idx = {xx:x for x,xx in enumerate(all_pathways)}
    p_vals_mat = np.ones((len(all_pathways_names), len(prior_set_conditions)*2))
    adj_p_vals_mat = np.ones((len(all_pathways_names), len(prior_set_conditions)*2))
    z_score_mat = np.zeros((len(all_pathways_names), len(prior_set_conditions)*2))
    directions_mat = np.zeros((len(all_pathways_names), len(prior_set_conditions)*2))
    rows_names = []
    for i in range(len(prior_set_conditions)*2):
        index = int(i//2)
        if np.mod(i, 2):
            p_vals_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = prop_scores_dict['p_vals'][index]
            adj_p_vals_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = prop_scores_dict['adj_p_vals'][index]
            if prop_scores_dict['direction'][index] is not None:
                directions_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = prop_scores_dict['direction'][index]
            if prop_scores_dict['z_score'][index] is not None:
                z_score_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = prop_scores_dict['z_score'][index]
            rows_names.append(prior_set_conditions[index] + '_prop')
        else:
            p_vals_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = fc_scores_dict['p_vals'][index]
            adj_p_vals_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = fc_scores_dict['adj_p_vals'][index]
            if prop_scores_dict['direction'][index] is not None:
                directions_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] =fc_scores_dict['direction'][index]
            if fc_scores_dict['z_score'][index] is not None:
                z_score_mat[[pathway_to_idx[x] for x in pathways_per_condition[index]], i] = fc_scores_dict['z_score'][index]
            rows_names.append(prior_set_conditions[index] + '_logfc')

    if prop_scores_dict['direction'][index][0] is None:
        directions_mat = None

    if res_type == '-log10(p)':
        res_mat = -np.log10(p_vals_mat)
    elif res_type == '-log10(adj_p)':
        res_mat = -np.log10(adj_p_vals_mat)
    elif res_type == 'z_score':
        res_mat = z_score_mat

    fig_out_dir = path.join(args.output_folder, fig_name)
    plot_enrichment_table(res_mat ,adj_p_vals_mat, directions_mat, all_pathways_names, fig_out_dir,
                          experiment_names=rows_names, title=title, res_type=res_type)