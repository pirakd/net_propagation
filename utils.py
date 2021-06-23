import pandas as pd
import networkx as nx
import numpy as np
import scipy.stats
import os
from os import path
from datetime import datetime
import mygene

mg = mygene.MyGeneInfo()
# Global Variables
PROPAGATE_ALPHA = 0.9
PROPAGATE_ITERATIONS = 200
PROPAGATE_EPSILON = 10 ** (-4)


# Convert list of prior symbols to ids
def convert_symbols_to_ids(prior_symbols):
    missing_ids = []
    prior_gene_dict = {}

    ncbi_query = mg.querymany(prior_symbols, scopes="symbol", fields=["entrezgene", "symbol"], species="human")
    for result in ncbi_query:
        if "entrezgene" in result.keys():
            prior_gene_dict[int(result["entrezgene"])] = result["symbol"]
        else:
            missing_ids.append(result)

    # delete some unrelated redundant genes names
    if 7795 in prior_gene_dict:
        del prior_gene_dict[7795]
    if 100187828 in prior_gene_dict:
        del prior_gene_dict[100187828]

    return prior_gene_dict


def read_network(network_filename):
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


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
    prior_set, prior_data = condition_fucntion(data_frame)
    return prior_set, prior_data


def bh_correction(p_values):
    p_vals_rank = scipy.stats.rankdata(p_values, 'max') - 1
    p_vals_rank_ord = scipy.stats.rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) /(p_vals_rank+1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_vals_rank]
    return adj_p_vals


def two_sample_z_test(mu_1, mu_2, mu_diff, sd_1, sd_2, n_1, n_2):
    from numpy import sqrt, abs, round
    from scipy.stats import norm
    pooledSE = sqrt(sd_1**2/n_1 + sd_2**2/n_2)
    z = ((mu_1 - mu_2) - mu_diff)/pooledSE
    pval = 2*(1 - norm.cdf(abs(z)))
    return z, pval

def create_output_folder(test_name):
    time = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
    output_folder = path.join('output', test_name, time)
    os.makedirs(output_folder)
    return output_folder


def listdir_nohidden(path):
    file_list = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            file_list.append(f)
    return file_list


def load_interesting_pathways(pathways_dir):
    with open(pathways_dir, 'r') as f:
        pathways = [str.replace(str.upper(x.strip()), ' ', '_') for x in f]
    return pathways


def load_pathways_genes(pathways_dir, interesting_pathways=None):
    with open(pathways_dir, 'r') as f:
        lines = [str.upper(x.strip()).split('\t') for x in f]
    pathways = {x[0]: [int(y) for y in x[2:]] for x in lines}

    if interesting_pathways:
        interesting_pathways = [str.upper(x) for x in interesting_pathways]
        filtered_pathways = {x: xx for x, xx in pathways.items() if x in interesting_pathways}

    return filtered_pathways


def get_propagation_input(prior_gene_dict, prior_data, input_type):
    if input_type == 'ones':
        inputs = {x:1 for x in prior_gene_dict.keys()}
    elif input_type is None:
        inputs = {x:1 for x in prior_gene_dict.keys()}
    elif input_type == 'abs_log2FC':
        inputs = {x: np.abs(float(prior_data[prior_data.Gene_Name == xx]['log2FC'])) for x, xx in prior_gene_dict.items()}
    else:
        assert 0, '{} is not a valid input type'.format(input_type)
    return inputs