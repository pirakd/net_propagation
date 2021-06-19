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

#Convert list of prior symbols to ids
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


def generate_similarity_matrix(network):
    genes = sorted(network.nodes)
    matrix = nx.to_scipy_sparse_matrix(network, genes, weight=2)
    norm_matrix = sp.sparse.diags(1 / sp.sqrt(matrix.sum(0).A1), format="csr")
    matrix = norm_matrix * matrix * norm_matrix
    return matrix, genes


def propagate(seeds, matrix, gene_indexes, num_genes):
    F_t = np.zeros(num_genes)
    F_t[[gene_indexes[seed] for seed in seeds if seed in gene_indexes]] = 1

    Y = PROPAGATE_ALPHA * F_t
    matrix = (1-PROPAGATE_ALPHA) * matrix
    for _ in range(PROPAGATE_ITERATIONS):
        F_t_1 = F_t
        F_t = matrix.dot(F_t_1) + Y

        if math.sqrt(sp.linalg.norm(F_t_1 - F_t)) < PROPAGATE_EPSILON:
            break
    return F_t


def generate_propagate_data(network, prior_set=None, interactors=None):
    matrix, genes = generate_similarity_matrix(network)
    num_genes = len(genes)
    gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])

    if interactors:
        gene_scores = {gene: propagate(
            [gene], matrix, gene_indexes, num_genes) for gene in interactors}
    if prior_set:
        gene_scores = {gene: propagate(
            [gene], matrix, gene_indexes, num_genes) for gene in prior_set}
    else:
        gene_scores = {gene: propagate(
            [gene], matrix, gene_indexes, num_genes) for gene in genes}

    return matrix, num_genes, gene_indexes, gene_scores


def load_network_scores(g, prior_set):
    print('start propagating network')
    W, num_genes, gene_indexes, gene_scores = generate_propagate_data(g, prior_set)
    network_scores = {"W": W, "num_genes": num_genes, "gene_indexes": gene_indexes, "gene_scores": gene_scores}
    return W, num_genes, gene_indexes, gene_scores


def load_random_networks(g, prior_set, random_networks_dir=None, n_networks=100):
    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = os.listdir(random_networks_dir)

    E = g.number_of_edges()
    Q = 10
    random_networks = {}

    if random_networks_files:
        networks_to_process =np.min([len(random_networks_files), n_networks])
        for n, network_file_name in enumerate(random_networks_files[:networks_to_process]):
            print('propagating network {}'.format(n))

            H = load_file(os.path.join(random_networks_dir, network_file_name))
            W_temp, num_genes_temp, gene_indexes_temp, gene_scores_temp = generate_propagate_data(H, prior_set)
            random_networks[n] = gene_scores_temp
    else:
        for n in range(n_networks):
            H = g.copy()
            nx.swap.double_edge_swap(H, nswap=Q * E, max_tries=Q * E * 2)
            W_temp, num_genes_temp, gene_indexes_temp, gene_scores_temp = generate_propagate_data(H, prior_set)
            random_networks[n] = gene_scores_temp
    return random_networks


def load_and_map_interactions():
    # https://www.nature.com/articles/s41586-020-2286-9#Sec36
    interactions = pd.read_csv("covid_files/data/inputs/interactions.csv")
    # mapping to entrez id
    xli = interactions["PreyGene"].unique().tolist()
    out = pd.DataFrame(mg.querymany(xli, scopes="symbol", fields="entrezgene", species="human"))
    interactions = pd.merge(interactions, out[["query", "entrezgene"]], left_on="PreyGene", right_on="query")
    interactions["entrezgene"] = interactions["entrezgene"].astype(np.float).astype("Int32")
    return interactions


def load_all_targets():
    all_targets = {'enterocytes': pd.read_csv("covid_files/data/inputs/Enterocytes.csv",
                                              dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
        'entrez'],
                   'proximal': pd.read_csv("covid_files/data/inputs/Proximal_tubule_cells.csv",
                                           dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
                       'entrez'],
                   'cardiomyocytes': pd.read_csv("covid_files/data/inputs/cardiomyocytes.csv",
                                                 dtype={'entrez': pd.Int64Dtype()}).dropna().set_index(
                       'gene').to_dict()['entrez'],
                   'bronchial': pd.read_csv("covid_files/data/inputs/human_bronchial_epithelial.csv",
                                            dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
                       'entrez'],
                   'lung': pd.read_csv("covid_files/data/inputs/lung.csv", dtype={'entrez': pd.Int64Dtype()}).set_index(
                       'gene').dropna().to_dict()['entrez'],
                   'lymphocytes': pd.read_csv("covid_files/data/inputs/lymphocytes.csv",
                                              dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
                       'entrez'],
                   'neuronal': pd.read_csv("covid_files/data/inputs/neuronal.csv",
                                           dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
                       'entrez'],
                   'vascular': pd.read_csv("covid_files/data/inputs/vascular.csv",
                                           dtype={'entrez': pd.Int64Dtype()}).dropna().set_index('gene').to_dict()[
                       'entrez']}
    return all_targets


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
