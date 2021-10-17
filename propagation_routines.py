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


def propagate_network(network, propagation_input, args:Args, genes=None):
    matrix, genes = generate_similarity_matrix(network, genes)
    num_genes = len(genes)
    gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])

    if args.remove_self_propagation:
        self_propagation = list()
    else:
        self_propagation = None

    if (not args.remove_self_propagation) or args.alpha==1:
        gene_scores = propagate([x for x in propagation_input.keys()], propagation_input, matrix, gene_indexes, num_genes, args)
    else:
        gene_scores = list()
        for x in list(propagation_input.keys()):
            gene_scores.append(propagate([x], propagation_input, matrix, gene_indexes, num_genes, args))
            self_propagation.append(gene_scores[-1][gene_indexes[x]])
            gene_scores[-1][gene_indexes[x]] = 0

        gene_scores = np.sum(np.array(gene_scores), axis=0)
        correction_factor = len(propagation_input.keys()) / (len(propagation_input.keys())-1)
        gene_scores[[gene_indexes[x] for x in list(propagation_input.keys())]] *= correction_factor

    return matrix, num_genes, gene_indexes, gene_scores, np.array(self_propagation)


def propagate_networks(network,  args:Args, genes=None, prior_set=None, propagation_input=None,
                       random_networks_dir=None, n_networks=100):
    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = listdir_nohidden(random_networks_dir)
    gene_scores = []
    self_propagations = []

    if random_networks_files:
        networks_to_process = np.min([len(random_networks_files), n_networks])

        for n, network_file_name in enumerate(random_networks_files[:networks_to_process]):
            print('propagating network {}'.format(n))
            H = load_file(os.path.join(random_networks_dir, network_file_name))
            matrix, num_genes, gene_indexes, current_gene_scores, current_self_propagation = propagate_network(H, propagation_input, args, genes)
            gene_scores.append(current_gene_scores)
            self_propagations.append(current_self_propagation)

    else:
        E = network.number_of_edges()
        Q = 10
        for n in range(n_networks):
            H = network.copy()
            nx.swap.double_edge_swap(H, nswap=Q * E, max_tries=Q * E * 2)
            matrix, num_genes, gene_indexes, current_gene_scores = propagate_network(H, propagation_input, genes=genes, prior_set=prior_set, args=args)
            gene_scores.append(current_gene_scores)

    return gene_indexes, np.array(gene_scores), self_propagations


def parallel_propagate(network_dir=None, propagation_input=None, args=None, genes=None):
    H = load_file(network_dir)
    matrix, num_genes, gene_indexes, current_gene_scores, self_propagation = propagate_network(H, propagation_input, args, genes)
    return [current_gene_scores, self_propagation]

def propagate_networks_parallel(network, args, genes=None, prior_set=None,
                                propagation_input=None, random_networks_dir=None, n_networks=100, n_processes=1):
    import multiprocessing as mp
    from functools import partial
    mp.set_start_method('fork')

    random_networks_files = None
    if random_networks_dir and os.path.isdir(random_networks_dir):
        random_networks_files = listdir_nohidden(random_networks_dir)

    n_networks_to_process = np.min([len(random_networks_files), n_networks])
    network_dirs = (os.path.join(random_networks_dir, random_networks_files[i]) for i in range(n_networks_to_process))

    with mp.Pool(processes=n_processes) as pool:
        results = \
            [i for i in tqdm.tqdm(pool.imap_unordered(partial(parallel_propagate, args=args, propagation_input=propagation_input, genes=genes),
                                                      network_dirs), total=n_networks_to_process, desc='propagating networks')]
    genes_scores, self_propagations = list(zip(*results))
    return np.array(genes_scores), self_propagations