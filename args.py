from os import path
from experiment_readers import get_experiment_reader
from datetime import datetime
from utils import create_output_folder


class Args:
    def __init__(self, test_name='test', is_create_output_folder=True):
        # here we define default parameters that will optionally be overwritten by child classes

        # ~~~ general parameters ~~~
        self.propagation_folder = 'propagation_scores'
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'inputs'
        self.network_file = 'H_sapiens.net'
        self.experiment_file = 'scores'
        self.propagation_input_type = 'ones'
        self.sheet_name = 'Table_A'
        self.experiment_reader_name = 'cov_data'
        self.interesting_pathway_file = 'interesting_pathways.txt'
        self.genes_names_file = 'H_sapiens_symbol'
        self.experiment_file = 'scores.xlsx'
        self.pathway_file = 'canonical_pathways.txt'
        self.random_network_file = 'random_networks'
        self.n_networks = 1000
        self.propagation_folder = 'propagation_scores'
        self.date = None

        # propagation parameters
        self.alpha = 0.95
        self.n_max_iterations = 200
        self.convergence_th = 1e-5
        self.test_name = test_name
        self.remove_self_propagation = True

        # post processing arguments
        self.normalization_method = 'ones'  #ones\EC
        self.add_self_prop_to_norm_factor = False

        # init derived params
        self.data_dir = None
        self.output_folder = None
        self.experiment_reader = None
        self.network_file_path = None
        self.experiment_file_path = None
        self.pathway_file_dir = None
        self.interesting_pathway_file_dir = None
        self.genes_names_file_path = None
        self.propagation_scores_path = None
        self.random_networks_dir = None
        self.get_derived_parameters(is_create_output_folder)

        # ~~~ derived parameters ~~~
    def get_derived_parameters(self, is_create_output_folder=True):
        self.set_experiment_reader()
        self.experiment_name = self.experiment_file.split('.')[0]

        self.data_dir = path.join(self.root_folder, self.data_file)
        if is_create_output_folder:
            self.output_folder, self.date = create_output_folder(self.test_name)
        else:
            self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
            self.output_folder = None

        self.experiment_reader = get_experiment_reader(self.experiment_reader_name)
        self.network_file_path = path.join(self.data_dir, 'networks', self.network_file)
        self.experiment_file_path = path.join(self.data_dir, 'experiments_data', self.experiment_file)
        self.pathway_file_dir = path.join(self.data_dir, 'pathways', self.pathway_file)
        if self.interesting_pathway_file:
            self.interesting_pathway_file_dir = path.join(self.data_dir, 'pathways', self.interesting_pathway_file)
        self.random_networks_dir = path.join(self.root_folder, self.random_network_file)
        if self.genes_names_file:
            self.genes_names_file_path = path.join(self.data_dir, 'genes_names', self.genes_names_file)
        else:
            self.genes_names_file_path = None
        self.propagation_scores_path = path.join(self.root_folder, self.propagation_folder)
        return self.experiment_name, self.data_dir, self.output_folder, self.experiment_reader, self.network_file,\
               self.experiment_file_path, self.pathway_file_dir, self.interesting_pathway_file_dir,\
               self.genes_names_file_path, self.propagation_scores_path

    def set_experiment_reader(self):
        self.experiment_reader = get_experiment_reader(self.experiment_reader_name)


class DeltaArgs(Args):
    def __init__(self, test_name=None, is_create_output_folder=True):
        super().__init__(is_create_output_folder=False)
        self.experiment_file = 'delta_scores.xlsx'
        self.network_file = 'H_sapiens.net'
        self.sheet_name = 'India2_24h-IC19_24h'
        self.propagation_input_type = 'Score'
        self.experiment_reader_name = 'cov_data'
        self.interesting_pathway_file = 'interesting_pathways.txt'
        self.genes_names_file = 'H_sapiens_uniprot'
        self.test_name = test_name
        self.remove_self_propagation = True
        self.get_derived_parameters(is_create_output_folder=is_create_output_folder)

class ProstateArgs(Args):
    def __init__(self, test_name=None, is_create_output_folder=True):
        super().__init__(is_create_output_folder=False)
        self.experiment_file = 'KEGG_TGF_BETA_SIGNALING_PATHWAY_10_25.xlsx'
        self.propagation_folder = 'propagation_scores_prostate'
        self.network_file = 'H_sapiens.net'
        self.sheet_name = 'Table_A'
        self.propagation_input_type = 'bootstrap_1'
        self.experiment_reader_name = 'prostate_data'
        self.interesting_pathway_file = 'interesting_pathways.txt'
        self.genes_names_file = 'H_sapiens_symbol'
        self.test_name = test_name
        self.remove_self_propagation = True
        self.get_derived_parameters(is_create_output_folder=is_create_output_folder)


class KentArgs(Args):
    def __init__(self, test_name=None, is_create_output_folder=True):
        super().__init__(is_create_output_folder=False)
        self.experiment_file = 'kent_scores.xlsx'
        self.network_file = 'H_sapiens.net'
        self.sheet_name = 'Kent_24h-VIC_24h'
        self.propagation_input_type = 'Score'
        self.experiment_reader_name = 'cov_data'
        self.interesting_pathway_file = 'interesting_pathways.txt'
        self.genes_names_file = 'H_sapiens_symbol'
        self.test_name = test_name
        self.remove_self_propagation = True
        self.propagation_folder = 'propagation_scores'
        self.get_derived_parameters(is_create_output_folder=is_create_output_folder)

class MockArgs(Args):
    def __init__(self, test_name='test', is_create_output_folder=False):
        super().__init__(is_create_output_folder=False)
        self.experiment_file = 'mock_scores.xlsx'
        self.network_file = 'H_sapiens.net'
        self.sheet_name = 'Table_A'
        self.propagation_input_type = 'Score'
        self.experiment_reader_name = 'cov_data'
        self.interesting_pathway_file = 'mock_interesting_pathways.txt'
        self.genes_names_file = 'H_sapiens_symbol'
        self.test_name = test_name
        self.remove_self_propagation = True
        self.get_derived_parameters(is_create_output_folder=is_create_output_folder)


class HDArgs(Args):
    def __init__(self, test_name=None, is_create_output_folder=True):
        super().__init__(is_create_output_folder=False)
        self.data_file = 'HD_data'
        self.network_file = 'HD.net'
        self.experiment_reader_name = 'huntington_DDA'
        self.sheet_name = 'Table_A'
        self.interesting_pathway_file = 'None'
        self.genes_names_file = 'HD_symbol'
        self.get_derived_parameters(is_create_output_folder=True)
        self.alpha = 0.9

class ColorectalArgs(Args):
    def __init__(self, test_name=None, is_create_output_folder=True):
        super().__init__(is_create_output_folder=False)
        self.data_file = 'colorectal_data'
        self.network_file = 'H_sapiens.net'
        self.experiment_file = 'scores.xlsx'
        self.experiment_reader_name = 'colorectal_cancer'
        self.propagation_input_type = 'abs_Score'
        self.sheet_name = 'EV_positive'
        self.interesting_pathway_file = None
        self.get_derived_parameters(is_create_output_folder=True)
        self.genes_names_file = 'H_sapiens_symbol'
        self.get_derived_parameters(is_create_output_folder=True)

if __name__ == '__main__':
    # instantiate an argument object with default parameter
    default_args = Args()

    # instantiate an argument object with covid preset
    cov_args = Args()
