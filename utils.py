import pandas as pd
import networkx as nx
import numpy as np
import scipy as sp
import math
import os
from time import time
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

    return prior_gene_dict


def read_network(network_filename):
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


def save_file(file, save_dir=None, compress=True):
    import pickle, zlib
    file = pickle.dumps(file)
    if compress:
        file = zlib.compress(file)
    with open(save_dir, 'wb') as f:
        pickle.dump(file, f)
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
    prior_set = condition_fucntion(data_frame)
    return prior_set

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