import pandas as pd
import networkx as nx
import numpy as np
import scipy as sp
import math
import os
import tqdm
from utils import listdir_nohidden
from time import time

# Global Variables
PROPAGATE_ALPHA = 0.9
PROPAGATE_ITERATIONS = 200
PROPAGATE_EPSILON = 10 ** (-4)

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

def generate_similarity_matrix(network):
    genes = sorted(network.nodes)
    matrix = nx.to_scipy_sparse_matrix(network, genes, weight=2)
    norm_matrix = sp.sparse.diags(1 / sp.sqrt(matrix.sum(0).A1), format="csr")
    matrix = norm_matrix * matrix * norm_matrix
    return matrix, genes

def propagate_network(network, prior_set=None):
    matrix, genes = generate_similarity_matrix(network)

    num_genes = len(genes)
    gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])
    if prior_set:
        gene_scores = propagate(prior_set, matrix, gene_indexes, num_genes)
    else:
        gene_scores =propagate(genes, matrix, gene_indexes, num_genes)

    return matrix, num_genes, gene_indexes, gene_scores


def propagate_networks(network, prior_set=None, random_networks_dir=None, n_networks=100):
    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = os.listdir(random_networks_dir)
    gene_scores = []

    if random_networks_files:
        networks_to_process = np.min([len(random_networks_files), n_networks])

        for n, network_file_name in enumerate(random_networks_files[:networks_to_process]):
            print('propagating network {}'.format(n))
            H = load_file(os.path.join(random_networks_dir, network_file_name))
            matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, prior_set)
            gene_scores.append(current_gene_scores)

    else:
        E = network.number_of_edges()
        Q = 10
        for n in range(n_networks):
            H = network.copy()
            nx.swap.double_edge_swap(H, nswap=Q * E, max_tries=Q * E * 2)
            matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, prior_set)
            gene_scores.append(current_gene_scores)

    return gene_indexes, np.array(gene_scores)


def parallel_propagate(network_dir, prior_set):
    H = load_file(network_dir)
    matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, prior_set)
    return current_gene_scores


def propagate_networks_parallel(network, prior_set=None, random_networks_dir=None, n_networks=100, n_processes):
    import multiprocessing as mp
    from functools import partial
    mp.set_start_method('fork')

    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = listdir_nohidden(random_networks_dir)
    gene_scores = []
    n_networks_to_process = np.min([len(random_networks_files), n_networks])
    pbar = tqdm.tqdm(total=n_networks_to_process)
    network_dirs = (os.path.join(random_networks_dir, random_networks_files[i]) for i in range(n_networks_to_process))

    network_dirs = (os.path.join(random_networks_dir, random_networks_files[i]) for i in range(n_networks_to_process))

    with mp.Pool(processes=4) as pool:
        gene_scores = \
            [i for i in tqdm.tqdm(pool.imap_unordered(partial(parallel_propagate, prior_set=prior_set), network_dirs),
                           total=n_networks_to_process, desc='propagating networks')]
    return np.array(gene_scores)

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

def read_prior_sets(excel_dir, sheet_name, conditions=None):
    # xls = pd.ExcelFile(excel_dir)
    data_frame = pd.read_excel(excel_dir, engine='openpyxl', sheet_name=sheet_name)

    prior_sets = dict()
    for condition in conditions:
        pos_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                    data_frame.Label == condition.expressed_in)].Gene_Name.unique())
        neg_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                    data_frame.Label == condition.not_expressed_in)].Gene_Name.unique())
        prior_sets[condition.name] = list(set.difference(pos_genes, neg_genes))

    return prior_sets

def get_genes_p_values(original_network_scores, random_networks_scores):

    n_experiments = random_networks_scores.shape[0]
    sorted_scores = np.sort(random_networks_scores, axis=0)
    gene_score_rank = np.array(
        [np.searchsorted(sorted_scores[:, i], original_network_scores[i]) for i in range(random_networks_scores.shape[1])])
    p_values = (n_experiments-gene_score_rank) / n_experiments
    return p_values