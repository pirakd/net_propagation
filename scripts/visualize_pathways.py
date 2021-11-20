import utils as utils
from utils import load_pathways_genes, load_interesting_pathways, get_propagation_input, load_file,\
    load_propagation_scores
from os import path
import numpy as np
from args import Args
test_name = path.basename(__file__).split('.')[0]
from experiment_readers import huntington_DDA_significant, huntington_DDA, colorectal_cancer_significant, colorectal_cancer, cov_data, cov_data_significant
from collections import namedtuple
from visualization_tools import visualise_pathway

Results = namedtuple('Results', 'p_value direction adj_p_value')
args = Args(test_name)
colls_names = ['Significant genes proportion test', 'Log2FC signed', 'Propagation signed',
              'Propagation unsigned']
# Huntington
interesting_database = ['REACTOME']
interesting_pathways_keywords = ['immune']

mark_second_neighbors = True
min_genes_in_pathway = 10
min_genes_for_display = 3
minimum_pathway_adjusted_p_value_for_display = 1

# huntingtons
all_genes_condition = cov_data
significant_genes_condition = cov_data_significant
signed_input = 'Score'
unsigned_input = 'abs_Score'

# load network
network_graph = utils.read_network(args.network_file_path)
network_genes_ids = set(list(network_graph.nodes()))
args.alpha = 0.95
args.sheet_name = 'India2_10h-IC19_10h'
title = 'Pathway Enrichment'
fig_name = 'enrichment_{:2d}.pdf'.format(np.int(args.alpha*100))

# read significant genes according to experiment
significant_genes, _, _ = significant_genes_condition(args)
args.experiment_reader = all_genes_condition
# get all genes list
all_genes, all_data, _ = all_genes_condition(args)
all_genes_dict = utils.convert_symbols_to_ids(all_genes, args.genes_names_file_path)

# get excel scores for all genes which appear in the network
reference_scores = get_propagation_input(all_genes_dict, all_data, 'Score', network=network_graph)
reference_scores = {x: xx for x, xx in reference_scores.items() if x in network_genes_ids}


# load propagation scores
args.propagation_input_type = unsigned_input
unsigned_gene_scores, _, _ = load_propagation_scores(args, normalize_score=True)
args.propagation_input_type = signed_input
signed_gene_scores, genes_idx_to_id, genes_id_to_idx = load_propagation_scores(args, normalize_score=True)

# load all MSigDB pathways and their associated genes
genes_by_pathway_all = load_pathways_genes(args.pathway_file_dir)
names_to_ids = utils.convert_symbols_to_ids(genes_names_file_path=args.genes_names_file_path)
ids_to_names = {xx: x for x, xx in names_to_ids.items()}

# filter pathway by DB
genes_by_pathway_interesting = {x: xx for x, xx in genes_by_pathway_all.items() if any(key in x for key in interesting_database)}

# filter pathway by keywords
genes_by_pathway_interesting = {x: xx for x, xx in genes_by_pathway_interesting.items() if any(str.upper(key) in x for key in interesting_pathways_keywords)}

interesting_pathways_names = list(genes_by_pathway_interesting.keys())
filtered_genes_by_pathway = {pathway: [id for id in genes_by_pathway_interesting[pathway] if id in genes_id_to_idx] for pathway
                                  in interesting_pathways_names}
reference_filtered_genes_by_pathway = {pathway: [id for id in filtered_genes_by_pathway[pathway] if id in reference_scores] for pathway
                                       in interesting_pathways_names}
pathways_with_many_genes = [pathway_name for pathway_name in filtered_genes_by_pathway.keys() if
                            (len(filtered_genes_by_pathway[pathway_name]) >= min_genes_in_pathway) and (len(filtered_genes_by_pathway[pathway_name]) <=50)]

# get scores of genes related to interesting pathway
prop_filtered_signed_genes_scores_by_filtered_pathways = {pathway: [signed_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                   in pathways_with_many_genes}
prop_filtered_unsigned_genes_scores_by_filtered_pathways = {pathway: [unsigned_gene_scores[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]] for pathway
                                                   in pathways_with_many_genes}
reference_filtered_genes_scores_by_filtered_pathways = {pathway: [reference_scores[id] for id in reference_filtered_genes_by_pathway[pathway]] for pathway
                                                        in pathways_with_many_genes}

# propotion test for significant experiment genes
network_genes_names = [ids_to_names[id] for id in genes_id_to_idx.keys() if id in ids_to_names]
all_pathways_genes_by_name = {pathway: [ids_to_names[id] for id in genes if (id in ids_to_names)]
                             for pathway, genes in genes_by_pathway_all.items()}
significant_in_net = [gene_name for gene_name in significant_genes if gene_name in network_genes_names]

sorted = np.argsort(np.array(list(reference_scores.values())))
ids = list(reference_scores.keys())
reference_scores_by_rank = dict()
for i in range(len(sorted)):
    reference_scores_by_rank[ids[sorted[i]]] = i

norm_func = lambda val, min_val, max_val, : (2 * ((val - min_val) / (max_val - min_val))) - 1
signed_gene_scores_sorted = np.argsort(signed_gene_scores)
signed_gene_scores_by_rank = np.zeros_like(signed_gene_scores_sorted)
for i in range(len(signed_gene_scores_sorted)):
    signed_gene_scores_by_rank[signed_gene_scores_sorted[i]] = i

reference_scores_normalized = {x: norm_func(xx, 0, len(reference_scores_by_rank)) for x, xx in
                               reference_scores_by_rank.items()}
signed_gene_scores_normalized = norm_func(signed_gene_scores_by_rank, 0, len(signed_gene_scores))

significant_in_net = [gene_name for gene_name in significant_genes if gene_name in network_genes_names]

for pathway in pathways_with_many_genes:
    genes_scores = {id: signed_gene_scores_normalized[genes_id_to_idx[id]] for id in filtered_genes_by_pathway[pathway]}
    save_dir = path.join(args.output_folder, '{}.pdf'.format(pathway))
    visualise_pathway(network_graph, filtered_genes_by_pathway[pathway], reference_scores_normalized,
                      genes_scores, pathway,
                      ids_to_names, mark_second_neighbors=True, significant_genes=significant_in_net, save_dir=save_dir )