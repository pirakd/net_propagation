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
        try:
                file = zlib.decompress(file)
                file = pickle.loads(file)
                return file
        except:
            print('entered an uncompressed file but asked decompress it')
            return file
    else:
        return file



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


def get_propagation_input(prior_gene_dict, prior_data, input_type, network):
    """
    :param prior_gene_dict: dictionary contains gene_name:gene_id pairs
    :param prior_data: all excel file
    :param input_type:
    :return:
    """

    if input_type == 'ones':
        inputs = {x: 1 for x in prior_gene_dict.values() if x in network.nodes}
    elif input_type is None:
        inputs = {x: 1 for x in prior_gene_dict.values() if x in network.nodes}
    elif input_type == 'abs_Score':
        inputs = {id: np.abs(float(prior_data[prior_data.Gene_Name == name]['Score'])) for name, id in prior_gene_dict.items()}
    elif input_type == 'Score':
        inputs = {id: float(prior_data[prior_data.Gene_Name == name]['Score']) for name, id in prior_gene_dict.items()}
    elif input_type == 'Score_all':
        inputs = {id: float(prior_data[prior_data.Gene_Name == name]['Score'])
                  for name, id in prior_gene_dict.items()}
        mean_input = np.mean([x for x in inputs.values()])
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = mean_input

    elif input_type == 'abs_Score_all':
        inputs = {id: np.abs(float(prior_data[prior_data.Gene_Name == name]['Score']))
                  for name, id in prior_gene_dict.items()}
        mean_input = np.mean([x for x in inputs.values()])
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = mean_input
    elif input_type == 'ones_all':
        inputs = dict()
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = 1
    else:
        try:
            inputs = {id: float(prior_data[prior_data.Gene_Name == name][input_type]) for name, id in prior_gene_dict.items()}
        except KeyError:
            assert 0, '{} is not a valid input type'.format(input_type)

    inputs = {id: input_score for id, input_score in inputs.items() if id in network.nodes}
    return inputs


def get_time():
    return datetime.today().strftime('%d_%m_%Y__%H_%M_%S')


def save_propagation_score(propagation_scores, prior_set, propagation_input, genes_idx_to_id, args,
                           self_propagation=None, randomization_ranks=None, n_randomizations=None,
                           scores_p_values=None, date=None, random_networks_prop_score=None,file_name=None,
                           save_dir = None):
    if date is None:
        date = args.date
    if file_name is None:
        file_name = '{}_{}_{}_{}_{}'.format(args.propagation_input_type, args.experiment_name, args.sheet_name,
                                         args.experiment_reader_name, str(args.alpha))

    if save_dir is None:
        os.makedirs(args.propagation_scores_path, exist_ok=True)
        propagation_results_path = path.join(args.propagation_scores_path, file_name)
    else:
        propagation_results_path = path.join(save_dir, file_name)

    save_dict = {'args': args, 'prior_set': prior_set, 'propagation_input': propagation_input,
                 'gene_idx_to_id': genes_idx_to_id, 'gene_prop_scores': propagation_scores,
                 'self_propagation': self_propagation, 'randomization_ranks':randomization_ranks,
                 'n_randomizations':n_randomizations, 'scores_p_values': scores_p_values}
    if random_networks_prop_score:
        save_dict['random_prop_scores'] = random_networks_prop_score

    save_file(save_dict, propagation_results_path)
    return save_dict

def load_propagation_scores(args, add_self_propagation=False, normalize_score = True, propagation_file_name=None,
                            normalization_file_name=None):

    if propagation_file_name is None:
        propagation_file_name = '{}_{}_{}_{}_{}'.format(args.propagation_input_type, args.experiment_name, args.sheet_name, args.experiment_reader_name,
                                                     str(args.alpha))

    propagation_results_path = path.join(args.propagation_scores_path, propagation_file_name)
    propagation_res_dict: dict = load_file(propagation_results_path, decompress=True)
    genes_scores = np.array(propagation_res_dict['gene_prop_scores'])
    propagation_input = propagation_res_dict['propagation_input']
    genes_idx_to_id = propagation_res_dict['gene_idx_to_id']
    genes_id_to_idx = {xx: x for x, xx in genes_idx_to_id.items()}
    if add_self_propagation:
        genes_scores[[genes_id_to_idx[id] for id in propagation_input]] += propagation_res_dict['self_propagation']
        propagation_res_dict['gene_prop_score'] = genes_scores

    if normalize_score:
        genes_normalized_scores, genes_idx_to_id, genes_id_to_idx = normalize_propagation_scores(genes_scores, genes_idx_to_id,                                                                           args, normalization_file_name)
        propagation_res_dict['gene_prop_score'] = genes_normalized_scores
    return propagation_res_dict


def normalize_propagation_scores(gene_scores, genes_idx_to_id, args, normalization_file_name=None):
    genes_id_to_idx = {xx: x for x, xx in genes_idx_to_id.items()}
    if normalization_file_name is None:
        if args.normalization_method == 'EC':
            normalization_file_name = '{}_{}_{}_{}_1'.format(args.propagation_input_type, args.experiment_name,
                                                             args.sheet_name, args.experiment_reader_name)
        else:
            normalization_file_name = 'ones_{}_{}_{}_{}'.format(args.experiment_name, args.sheet_name,
                                                                args.experiment_reader_name, args.alpha)

    propagation_norm_res_path = path.join(args.propagation_scores_path, normalization_file_name)
    norm_propagation_res_dict = load_file(propagation_norm_res_path, decompress=True)

    norm_genes_idx_to_id = norm_propagation_res_dict['gene_idx_to_id']
    assert norm_genes_idx_to_id == genes_idx_to_id, 'Normalization scores did not come from the same network'

    propagation_input = norm_propagation_res_dict['propagation_input']

    norm_scores = np.array(norm_propagation_res_dict['gene_prop_scores'])
    if args.normalization_method != 'EC' and args.add_self_prop_to_norm_factor:
        norm_scores[[genes_id_to_idx[id] for id in propagation_input]] += norm_propagation_res_dict['self_propagation']

    zero_normalization_genes = np.nonzero(norm_scores == 0)[0]
    zero_prop_genes = np.nonzero(gene_scores == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_prop_genes))
    norm_scores[genes_to_delete] = 1
    gene_scores[gene_scores != 0] = np.array(gene_scores[gene_scores != 0] / np.abs(norm_scores[gene_scores != 0]))

    for gene_idx in genes_to_delete:
        gene_id = genes_idx_to_id[gene_idx]
        genes_idx_to_id.pop(gene_idx)
        genes_id_to_idx.pop(gene_id)

    return gene_scores, genes_idx_to_id, genes_id_to_idx

def get_root_path():
    return path.dirname(path.realpath(__file__))