import pandas as pd
import networkx as nx
import numpy as np
import os
import pickle as pl
from os import path
from datetime import datetime
import mygene
import json

mg = mygene.MyGeneInfo()


# Convert list of prior symbols to ids
def convert_symbols_to_ids(prior_symbols=None, genes_names_file_path=None):
    if genes_names_file_path:
        assert path.isfile(genes_names_file_path), 'Could not find {}'.format(genes_names_file_path)
        genes_names_to_ids = load_genes_ids_from_file(prior_symbols, genes_names_file_path)
    else:
        genes_names_to_ids = load_genes_ids_from_net(prior_symbols)
    return genes_names_to_ids

def load_genes_ids_from_file(genes_names, names_file):
    with open(names_file, 'r') as f:
        all_genes_names_to_ids = json.loads(f.read())
    if genes_names is not None:
        genes_names_to_ids = {name: int(all_genes_names_to_ids[name]) for name in genes_names if name in all_genes_names_to_ids}
    else:
        genes_names_to_ids = {name: int(all_genes_names_to_ids[name]) for name in all_genes_names_to_ids}
    return genes_names_to_ids


def load_genes_ids_from_net(prior_symbols):
    missing_names = []
    prior_gene_dict = {}

    ncbi_query = mg.querymany(prior_symbols, scopes="symbol", fields=["entrezgene", "symbol"], species="human")
    for result in ncbi_query:
        if "entrezgene" in result.keys():
            if int(result["entrezgene"] not in [7795, 100187828]):
                prior_gene_dict[result["symbol"]] = int(result["entrezgene"])
        else:
            missing_names.append(result['query'])

    return prior_gene_dict


def read_network(network_filename):
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)

def read_network_create_graph(network_filename, network_sheetname):
    network_table = pd.read_excel(network_filename, sheet_name=network_sheetname)
    df = pd.DataFrame(network_table, columns=['gene id a', 'gene id b'])
    network_graph = nx.Graph()
    for index, row in df.iterrows():
        network_graph.add_edge(row['gene id a'], row['gene id b'], weight=1)
    return network_graph


def save_file(obj, save_dir=None, compress=True):
    import pickle, zlib
    obj = pickle.dumps(obj)
    if compress:
        obj = zlib.compress(obj)
    with open(save_dir, 'wb') as f:
        pickle.dump(obj, f)
    print('File was saved in {}'.format(save_dir))


def load_file(load_dir, decompress=True):
    import pickle, zlib
    with open(load_dir, 'rb') as f:
        file = pickle.load(f)
    if decompress:
        file = zlib.decompress(file)
        file = pickle.loads(file)
    return file


def read_prior_set(condition_fucntion, excel_dir, sheet_name,):
    # xls = pd.ExcelFile(excel_dir)
    data_frame = pd.read_excel(excel_dir, engine='openpyxl', sheet_name=sheet_name)
    prior_set, prior_data, reference_data = condition_fucntion(data_frame)
    return prior_set, prior_data, reference_data


def create_output_folder(test_name):
    time = get_time()
    output_folder = path.join('output', test_name, time)
    os.makedirs(output_folder)
    return output_folder, time


def listdir_nohidden(path):
    file_list = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            file_list.append(f)
    return file_list


def load_interesting_pathways(pathways_dir):
    with open(pathways_dir, 'r') as f:
        pathways = [str.replace(str.upper(x.strip()), ' ', '_') for x in f]
    pathways = list(dict.fromkeys(pathways))
    return pathways


def load_pathways_genes(pathways_dir, interesting_pathways=None):
    with open(pathways_dir, 'r') as f:
        lines = [str.upper(x.strip()).split('\t') for x in f]
    pathways = {x[0]: [int(y) for y in x[2:]] for x in lines}

    if interesting_pathways is not None:
        interesting_pathways = [str.upper(x) for x in interesting_pathways]
        filtered_pathways = {pathway: pathways[pathway] if pathway in pathways else [] for pathway in interesting_pathways}
    else:
        filtered_pathways=pathways

    return filtered_pathways


def get_propagation_input(prior_gene_dict, prior_data, input_type):
    if input_type == 'ones':
        inputs = {x: 1 for x in prior_gene_dict.values()}
    elif input_type is None:
        inputs = {x: 1 for x in prior_gene_dict.values()}
    elif input_type == 'abs_log2FC':
        inputs = {id: np.abs(float(prior_data[prior_data.Gene_Name == name]['log2FC'])) for name, id in prior_gene_dict.items()}
    elif input_type == 'log2FC':
        inputs = {id: float(prior_data[prior_data.Gene_Name == name]['log2FC']) for name, id in prior_gene_dict.items()}
    elif input_type == 'Absolute AVG Log2 Ratio':
        inputs = {id: float(prior_data[prior_data.Gene == name]['Absolute AVG Log2 Ratio']) for name, id in prior_gene_dict.items()}
    elif input_type == 'AVG Log2 Ratio':
        inputs = {id: float(prior_data[prior_data.Gene == name]['AVG Log2 Ratio']) for name, id in prior_gene_dict.items()}
    elif input_type == 'Absolute Log2FC (HD/C116)':
        inputs = {id: float(prior_data[prior_data.Gene == name]['Absolute Log2FC (HD/C116)']) for name, id in prior_gene_dict.items()}
    elif input_type == 'Log2FC (HD/C116)':
        inputs = {id: float(prior_data[prior_data.Gene == name]['Log2FC (HD/C116)']) for name, id in prior_gene_dict.items()}
    else:
        assert 0, '{} is not a valid input type'.format(input_type)
    return inputs


def get_time():
    return datetime.today().strftime('%d_%m_%Y__%H_%M_%S')


def save_propagation_score(file_name, propagation_scores, prior_set, propagation_input, genes_idx_to_id, args, date=None):
    if date is None:
        date = args.date

    # save propagation score
    os.makedirs(args.propagation_scores_path, exist_ok=True)
    propagation_results_path = path.join(args.propagation_scores_path, file_name)

    args_dict = {'alpha': args.alpha, 'n_max_iterations': args.n_max_iterations, 'convergence_th': args.convergence_th,
                 'date': date}
    save_dict = {'propagation_args': args_dict, 'prior_set': prior_set, 'propagation_input': propagation_input,
                 'gene_idx_to_id': genes_idx_to_id, 'gene_prop_scores': propagation_scores}
    save_file(save_dict, propagation_results_path)

