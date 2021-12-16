from os import path
from utils import get_root_path, create_output_folder
import subprocess
import pandas as pd

class PropagationTask:
    def __init__(self, name, propagation_file, propagation_folder, normalization_file, statistic_test, target_field,
                 normalization_method, constrain_to_experiment_genes, add_self_propagation=False,
                 normalize_scores=True):
        self.name = name
        self.propagation_file = propagation_file
        self.propagation_scores_path = path.join(get_root_path(), propagation_folder)
        self.normalization_file = normalization_file
        self.statistic_test = statistic_test
        self.target_field = target_field
        self.normalize_scores = normalize_scores
        self.normalization_method = normalization_method
        self.constrain_to_experiment_genes = constrain_to_experiment_genes
        self.add_self_propagation = add_self_propagation
        self.add_self_propagation_to_norm_factor = add_self_propagation
        self.add_self_propagation_to_norm_factor = False
        self.results = dict()


class RawScoreTask:
    def __init__(self, name, experiment_file_path, sheet_name, statistic_test, experiment_reader, propagation_input_type,
                  constrain_to_network_genes=True):
        self.name = name
        self.experiment_file_path = experiment_file_path
        self.sheet_name = sheet_name
        self.statistic_test = statistic_test
        self.experiment_reader = experiment_reader
        self.propgation_input_type = propagation_input_type
        self.constrain_to_experiment = constrain_to_network_genes
        self.results = dict()


class GeneralArgs:
    def __init__(self, network_path, genes_names_path, interesting_pathway_path, pathway_members_path,
                 output_folder_name=None, figure_name=None, figure_title='Pathway Enrichment',
                 merge_similar_pathways=True):
        root_path = get_root_path()
        self.minimum_genes_per_pathway = 1
        self.display_only_significant_pathways = True
        self.network_file_path = network_path
        self.genes_names_file_path = genes_names_path
        self.pathway_databases = ['_']
        self.pathway_keywords = ['_']
        self.significant_pathway_threshold = 1
        if output_folder_name is None:
            output_folder_name = path.basename(__file__).split('.')[0]
        self.output_path, _ = create_output_folder(output_folder_name)
        self.figure_name = figure_name if figure_name is not None else 'figure'
        self.interesting_pathway_path = interesting_pathway_path
        self.pathway_members_path = pathway_members_path
        self.figure_title = figure_title
        self.simplify_dir = path.join(root_path, 'external_code')
        self.simplify_script_path = path.join(self.simplify_dir, 'simplify.R')
        self.merge_similar_pathways = merge_similar_pathways
        self.maximum_number_of_pathways = 2


class PathwayResults:
    def __init__(self):
        self.p_value = None
        self.direction = None
        self.adj_p_value = None


def filter_pathways(p_values, adj_p_values, pathway_names, n_genes, general_args:GeneralArgs):
    data_list = []
    n_groups = len(n_genes)

    for p, pathway in enumerate(pathway_names):
        for e in range(n_groups):
            if n_genes[e][pathway] > general_args.minimum_gene_per_pathway:
                id = pathway
                group = e
                p_val = p_values[p, e]
                adj_p_val = adj_p_values[p, e]
                n_genes_in_entry = n_genes[e][pathway]

                data_list.append([id, group, p_val, adj_p_val, n_genes_in_entry])

    df = pd.DataFrame(data_list, columns=['ID', 'group', 'pvalue', 'p.adjust', 'Count'])

    enrichment_path = path.join(general_args.output_path, 'enrichment.txt')
    pathway_path = general_args.pathway_members_path
    out_path = path.join(general_args.output_path, 'enrichment_simplified')
    cutoff_height = 0.8

    df.to_csv(path.join(general_args.output_path, 'enrichment.txt'), sep='\t', index_label=False)
    subprocess.run(['Rscript', general_args.simplify_script_path, pathway_path, enrichment_path, out_path, str(cutoff_height)],
                   stdout=subprocess.DEVNULL)
    simplified_pathways = pd.read_csv(out_path, sep='\t')

    b = simplified_pathways.sort_values(by=['group', 'cluster','Count'])
    c = b.groupby(by=['group','cluster'], as_index=False).last().sort_values(by=['group', 'p.adjust'])
    new_pathways_to_display = pd.unique(c['ID'])
    indexes_to_keep = [x for x,xx in enumerate(pathway_names) if xx in new_pathways_to_display]
    return indexes_to_keep
