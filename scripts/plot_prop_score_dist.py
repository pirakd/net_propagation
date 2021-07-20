import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, load_file
from os import path, makedirs
from propagation_routines import propagate_network
from os import path
import numpy as np
from args import Args
import pickle as pl
from statistic_methods import empirical_mean_diff, man_whit_U_test
test_name = path.basename(__file__).split('.')[0]

n_distributions = 5
n_tests = 2
normalize_by_eig_vec_cent_list = [False, True]
normalize_by_degree_list = [False, False]
from statistic_methods import get_stat_test_func
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches




res_type_list = ['-log10(p)'] * 2
log_of_prop_scores = True
# sheet_names_list = ['Protein_Abundance', 'Protein_Abundance', 'RNA', 'RNA']
sheet_names_list = ['Suppl. Table 4A'] * 2
statistic_test_list = ['man_whit_U_test'] * 2
prior_set_conditions = ['huntington_DIA', 'huntington_DDA']
reference_score_type_list = ['AVG Log2 Ratio','Log2FC (HD/C116)']
args = Args(test_name)
load_prop_results = True
keep_only_excel_genes = False
network_graph = utils.read_network(args.network_file)

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

        genes_by_pathway_filtered =  [[id for id in  genes_of_interesting_pathways[pathway] if id in genes_id_to_idx] for pathway
                                   in interesting_pathways]
        pathways_with_many_genes = [x for x, xx in enumerate(genes_by_pathway_filtered) if len(xx) >= 15]
        genes_by_pathway_filtered = [genes_by_pathway_filtered[x] for x in pathways_with_many_genes]

        # check scores of genes related to interesting pathway
        genes_scores_by_pathway = [[gene_scores[genes_id_to_idx[id]] for id in gene_set] for gene_set
                                   in genes_by_pathway_filtered]
        degrees = network_graph.degree(weight=2)
        degrees = {x:xx for x,xx in degrees}
        genes_degree_by_pathway = [[degrees[id] for id in gene_set] for gene_set
                                   in genes_by_pathway_filtered]

        log2fc_scores_by_pathway = [[all_reference_scores[id] for id in gene_set if id in all_reference_scores] for gene_set
                                    in genes_by_pathway_filtered]

        fc_p_vals, fc_dir, prop_p_vals, prop_dir = [], [], [], []
        for g, gene_set in enumerate(genes_by_pathway_filtered):
            random_genes = [id for id in all_reference_scores if id not in gene_set]

            fc_res = statistic_test(log2fc_scores_by_pathway[g], [all_reference_scores[id] for id in random_genes])
            fc_p_vals.append(fc_res.p_value)
            fc_dir.append(fc_res.directionality)

            random_genes = [id for id in genes_id_to_idx.keys() if id not in gene_set]
            prop_res = statistic_test(genes_scores_by_pathway[g], [gene_scores[genes_id_to_idx[id]]
                                                                     for id in random_genes])
            # prop_res = statistic_test(genes_degree_by_pathway[g], [degrees[id] for id in random_genes])


            prop_p_vals.append(prop_res.p_value)
            prop_dir.append(prop_res.directionality)

        sorted_pathways_by_p_value_idx = np.argsort(prop_p_vals)[:n_distributions]
        gene_scores_log, gene_degree_log, prop_reference_degree, prop_reference_score = [], [], [], []

        for g in sorted_pathways_by_p_value_idx:
            gene_set = genes_by_pathway_filtered[g]
            random_genes = [id for id in genes_id_to_idx.keys() if id not in gene_set]

            gene_degree_log.append(np.log10(genes_degree_by_pathway[g]))
            gene_scores_log.append(genes_scores_by_pathway[g])

            prop_reference_score.append([gene_scores[genes_id_to_idx[id]] for id in random_genes])
            prop_reference_degree.append(np.log10([degrees[id] for id in random_genes if id in degrees]))



        cols = ['Score distribution', 'Degree and Propagation Score']
        rows = [interesting_pathways[pathways_with_many_genes[x]] for x in sorted_pathways_by_p_value_idx]
        # create 3x1 subfigs

        fig, axes = plt.subplots(nrows=n_distributions, ncols=2, figsize=(12, 8))
        # plt.suptitle('Score Distribution, {}'.format(title))
        for ax, col in zip(axes[0, :], cols):
            ax.set_title(col)

        for ax, row in zip(axes[:,0], rows):
            ax.set_ylabel(row, rotation=0, size='small', ha='right')

        for d in range(axes.shape[0]):
            # create 1x3 subplots per subfig

            bins = np.maximum(int(len(gene_degree_log[d]) // 2), 10)
            axes[d, 0].hist(gene_scores_log[d], density=True, bins=bins, alpha=0.5, label ='pathway prop scores')
            axes[d, 0].hist(prop_reference_score[d], density=True, bins=bins, alpha=0.5, label ='all genes prop scores')
            axes[d, 0].set_xticks([])
            axes[d, 0].set_yticks([])

            handles = [mpatches.Patch(color='none', label='p value: {:.2e}'.format(prop_p_vals[sorted_pathways_by_p_value_idx[d]]))]
            legend2 = axes[d,0].legend(handles=handles, loc=1, )
            if d == 0:
                legend1 = axes[d,0].legend(loc='upper left')
                axes[d,0].add_artist(legend2)


            axes[d, 1].scatter(gene_degree_log[d],gene_scores_log[d])
            axes[d, 1].set_xlabel('gene degree')
            axes[d, 1].set_ylabel('prop score')
            axes[d, 1].set_xticks([])
            axes[d, 1].set_yticks([])


        plt.tight_layout()
        fig_out_dir = path.join(args.output_folder, title.replace('.','')+ '{}_{}'.format(c, i))
        plt.savefig(fig_out_dir, bbox_inches='tight')
