from os import path
from utils import create_output_folder
from prior_conditions import get_condition_function
from datetime import datetime


class Args:
    def __init__(self, test_name=None, is_create_output_folder=True):
        # ~~~ general parameters ~~~
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'data'

        # # corona virus
        # self.network_file = 'H_sapiens.net'
        # self.experiment_file = 'Table_S1_V1.xlsx'
        # self.propagation_input_type = 'abs_log2FC' #ones, logfc
        # self.condition_function_name = 'kent_mock_no_vic_mock_24h'
        # self.interesting_pathway_file = 'interesting_pathways.txt'

        # # Huntington
        # self.network_file = 'HD.net'
        # self.experiment_file = 'HD_scores.xlsx'
        # self.condition_function_name = 'huntington_DDA_significant'
        # self.propagation_input_type = 'Absolute Log2FC (HD/C116)'
        # self.sheet_name = 'Suppl. Table 4A'
        # self.interesting_pathway_file = 'interesting_pathways_HD.txt'

        # colorectal_data
        self.data_file = 'colorectal_data'
        self.network_file = 'H_sapiens.net'
        self.experiment_file = 'fc_scores.xlsx'
        self.condition_function_name = 'colorectal_cancer'
        self.propagation_input_type = 'abs_log2FC'
        self.sheet_name = 'EV_negative'
        self.interesting_pathway_file = None

        self.pathway_file = 'canonical_pathways.txt'
        self.random_network_file = 'random_networks'
        self.genes_names_file = 'genes_names_to_ids'
        self.n_networks = 1000
        self.propagation_folder = 'propagation_scores'
        self.date = None

        # propagation parameters
        self.alpha = 0.9
        self.n_max_iterations = 500
        self.convergence_th = 1e-8


        # ~~~ derived parameters ~~~
        if test_name is None:
            self.test_name = 'classical_enrichment_{}'.format(self.condition_function_name)
        else:
            self.test_name = test_name

        self.data_dir = path.join(self.root_folder, self.data_file)
        if is_create_output_folder:
            self.output_folder, self.date = create_output_folder(self.test_name)
        else:
            self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
            self.output_folder = None

        self.condition_function = get_condition_function(self.condition_function_name)
        self.network_file = path.join(self.data_dir, self.network_file)
        self.experiment_file_path = path.join(self.data_dir, self.experiment_file)
        self.pathway_file_dir = path.join(self.data_dir, self.pathway_file)

        if self.interesting_pathway_file:
            self.interesting_pathway_file_dir = path.join(self.data_dir, self.interesting_pathway_file)
        self.random_networks_dir = path.join(self.root_folder, self.random_network_file)
        self.genes_names_file_path = path.join(self.data_dir, self.genes_names_file)
        self.propagation_scores_path = path.join(self.root_folder, self.propagation_folder)


    def set_condition_function(self):
        self.condition_function = get_condition_function(self.condition_function_name)
