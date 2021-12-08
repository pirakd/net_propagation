from methods import fetch_experiment_data, change_expression
import pandas as pd
import numpy as np
import utils
from args import Args
from external_code.microarray_simulator.methods import get_statistic_score

args = Args(is_create_output_folder=False)
network = utils.read_network(args.network_file_path)
network_nodes = network.nodes

pathway_to_change = 'KEGG_TGF_BETA_SIGNALING_PATHWAY'
# pathway_to_change = 'REACTOME_COMPLEMENT_CASCADE'
change_prop = 0.25
change_amount = 0.2
n_bootstraps = 10

raw_data_file = 'raw_data.txt'
true_labels_file = 'labels.txt'

raw_data, labels, genes_names, symbol_to_id, id_to_idx = fetch_experiment_data(raw_data_file, true_labels_file, args)
n_samples = len(labels)

pathway_genes = utils.load_pathways_genes(args.pathway_file_dir,[pathway_to_change])
gene_candidates_to_change = pathway_genes[pathway_to_change]
gene_candidates_idx_to_change = [id_to_idx[x] for x in gene_candidates_to_change if x in id_to_idx]

num_of_genes_changed = np.ceil(change_prop * len(gene_candidates_to_change))
genes_idx_to_change = np.random.choice(gene_candidates_idx_to_change, int(num_of_genes_changed), replace=False)

# shuffled_labels = np.random.permutation(labels)
changed_data = change_expression(raw_data, labels, genes_idx_to_change, change_amount)

from statistic_methods import wilcoxon_rank_sums_test
pathway_indexes = [id_to_idx[x] for x in gene_candidates_to_change if x in id_to_idx]

score = get_statistic_score(changed_data, labels, type='s2n')
b = wilcoxon_rank_sums_test(score[pathway_indexes], score)

df_dict = {}
df_dict['Gene_Name'] = [x for x in genes_names if x in symbol_to_id]
for i in range(n_bootstraps):
    suffeled_indexes = np.random.choice(n_samples, replace=True, size=n_samples)
    shuffled_labels = labels[suffeled_indexes]
    shuffled_data = changed_data[:, suffeled_indexes]
    df_dict['bootstrap_{}'.format(i)] = get_statistic_score(shuffled_data, shuffled_labels, type='s2n')

df = pd.DataFrame(df_dict)
df.to_excel('{}_{}_{}.xlsx'.format(pathway_to_change, int(change_prop*100), int(change_amount*100)), index=False, sheet_name='Table_A')

# signal_to_noise = calc_signal_to_noise(raw_data, classes)
# t_test = calc_t_test(raw_data, classes)
# data = np.hstack([signal_to_noise[:, np.newaxis], t_test[:, np.newaxis]])
# # df = pd.DataFrame(data, index=genes_names,  columns=['genes_names', 'signal_to_noise', 't_test'])
#
# df_dict = {'Gene_Name': genes_names, 'signal_to_noise': signal_to_noise}
# df = pd.DataFrame(df_dict)
#
# with pd.ExcelWriter('test.xlsx') as writer:
#     for i in range(1):
#         df.to_excel(writer, sheet_name='Table_A'.format(i), index=False)