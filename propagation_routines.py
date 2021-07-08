import pandas as pd
import networkx as nx
import numpy as np
import scipy as sp
import os
import tqdm
from utils import listdir_nohidden, load_file
import scipy.sparse
from args import Args

def propagate(seeds, propagation_input, matrix, gene_indexes, num_genes, args:Args):
    F_t = np.zeros(num_genes)
    if not propagation_input:
        propagation_input = {x: 1 for x in seeds}
    for seed in seeds:
        if seed in gene_indexes:
            F_t[gene_indexes[seed]] = propagation_input[seed]
    Y = (1-args.alpha) * F_t
    matrix = args.alpha * matrix
    for _ in range(args.n_max_iterations):
        F_t_1 = F_t
        F_t = matrix.dot(F_t_1) + Y

        if sp.linalg.norm(F_t_1 - F_t) < args.convergence_th:
            break
    F_t[[gene_indexes[seed] for seed in seeds if seed in gene_indexes]] = 0
    return F_t


def generate_similarity_matrix(network, genes=None):
    if not sp.sparse.issparse(network):
        genes = sorted(network.nodes())
        matrix = nx.to_scipy_sparse_matrix(network, genes, weight=2)
    else:
        matrix = network
        assert genes, 'must enter genes id mapping to indexes when not loading them from file'
    norm_matrix = sp.sparse.diags(1 / sp.sqrt(matrix.sum(0).A1), format="csr")
    matrix = norm_matrix * matrix * norm_matrix
    return matrix, genes


def propagate_network(network, propagation_input, args:Args, genes=None, prior_set=None):
    matrix, genes = generate_similarity_matrix(network, genes)
    num_genes = len(genes)
    gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])
    effective_prior_set = [x for x in prior_set if x in gene_indexes]
    if prior_set:
        gene_scores = [propagate([x], propagation_input, matrix, gene_indexes, num_genes, args) for x in effective_prior_set]
        gene_scores = np.sum(np.array(gene_scores), axis=0)
        gene_scores[[gene_indexes[x] for x in effective_prior_set]] *= len(effective_prior_set) / (len(effective_prior_set)-1)
    else:
        gene_scores =propagate(genes, propagation_input=None, matrix=matrix, gene_indexes=gene_indexes,
                               num_genes=num_genes, args=args)

    return matrix, num_genes, gene_indexes, gene_scores


def propagate_networks(network,  args:Args, genes=None, prior_set=None, propagation_input=None,
                       random_networks_dir=None, n_networks=100):
    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = listdir_nohidden(random_networks_dir)
    gene_scores = []

    if random_networks_files:
        networks_to_process = np.min([len(random_networks_files), n_networks])

        for n, network_file_name in enumerate(random_networks_files[:networks_to_process]):
            print('propagating network {}'.format(n))
            H = load_file(os.path.join(random_networks_dir, network_file_name))
            matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, propagation_input, args, genes, prior_set)
            gene_scores.append(current_gene_scores)

    else:
        E = network.number_of_edges()
        Q = 10
        for n in range(n_networks):
            H = network.copy()
            nx.swap.double_edge_swap(H, nswap=Q * E, max_tries=Q * E * 2)
            matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, prior_set, args=args)
            gene_scores.append(current_gene_scores)

    return gene_indexes, np.array(gene_scores)


def parallel_propagate(network_dir=None, genes=None,  prior_set=None):
    H = load_file(network_dir)
    matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, genes, prior_set)
    return current_gene_scores


def propagate_networks_parallel(network, genes=None, prior_set=None,
                                random_networks_dir=None, n_networks=100, n_processes=1):
    import multiprocessing as mp
    from functools import partial
    mp.set_start_method('fork')

    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = listdir_nohidden(random_networks_dir)

    n_networks_to_process = np.min([len(random_networks_files), n_networks])
    network_dirs = (os.path.join(random_networks_dir, random_networks_files[i]) for i in range(n_networks_to_process))

    with mp.Pool(processes=4) as pool:
        gene_scores = \
            [i for i in tqdm.tqdm(pool.imap_unordered(partial(parallel_propagate, genes=genes, prior_set=prior_set),
                                                      network_dirs), total=n_networks_to_process, desc='propagating networks')]
    return np.array(gene_scores)


def save_file(file, save_dir=None, compress=True):
    import pickle, zlib
    file = pickle.dumps(file)
    if compress:
        file = zlib.compress(file)
    with open(save_dir, 'wb') as f:
        pickle.dump(file, f)
    print('File was saved in {}'.format(save_dir))


def read_prior_sets(excel_dir, sheet_name, conditions=None):
    # xls = pd.ExcelFile(excel_dir)
    data_frame = pd.read_excel(excel_dir, engine='openpyxl', sheet_name=sheet_name)

    prior_sets = dict()
    for condition in conditions:
        pos_genes = set(data_frame[(data_frame.diffexpressed is True) & (
                    data_frame.Label == condition.expressed_in)].Gene_Name.unique())
        neg_genes = set(data_frame[(data_frame.diffexpressed is True) & (
                    data_frame.Label == condition.not_expressed_in)].Gene_Name.unique())
        prior_sets[condition.name] = list(set.difference(pos_genes, neg_genes))

    return prior_sets