from os import path
from utils import create_output_folder
from prior_conditions import get_condition_function

class Args:
    def __init__(self, test_name=None):
        # ~~~ general parameters ~~~
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'data'
        self.network_file = 'H_sapiens.net'
        self.experiment_file = 'Table_S1_V1.xlsx'
        self.sheet_name = 'Protein_Abundance'
        self.condition_function_name = 'kent_mock_no_vic_mock_24h'
        self.propagation_input_type = 'abs_log2FC'  # ones, logfc
        self.pathway_file = 'canonical_pathways.txt'
        self.interesting_pathway_file = 'interesting_pathways.txt'
        self.random_network_file = 'random_networks'
        self.n_networks = 1000
        # ~~~ derived parameters ~~~
        if test_name is None:
            self.test_name = 'classical_enrichment_{}'.format(self.condition_function_name)
        else:
            self.test_name = test_name

        self.data_dir = path.join(self.root_folder, 'data')
        self.output_folder = create_output_folder(self.test_name)
        self.condition_function = get_condition_function(self.condition_function_name)
        # get conditions on experiment
        self.network_file = path.join(self.data_dir, self.network_file)
        self.experiment_file_path = path.join(self.data_dir, self.experiment_file)
        self.pathway_file_dir = path.join(self.data_dir, self.pathway_file)
        self.interesting_pathway_file_dir = path.join(self.data_dir, self.interesting_pathway_file)
        self.random_networks_dir = path.join(self.root_folder, self.random_network_file)
    def set_condition_function(self):
        self.condition_function = get_condition_function(self.condition_function_name)
