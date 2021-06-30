import utils as utils
from utils import read_prior_set, load_pathways_genes, load_interesting_pathways, get_propagation_input, bh_correction
from os import path
from propagation_routines import propagate_network, get_genes_p_values
from visualization_tools import plot_enrichment_table
import numpy as np
from args import Args
import pickle as pl
import mygene
import json

def load_names_from_net(gene_ids, genes_file_path):
    mg = mygene.MyGeneInfo()
    # query_output = mg.getgenes(gene_ids)
    query_output = mg.querymany(genes_ids, scopes="entrezgene", fields=["entrezgene", "symbol"], species="human")

    ids_found, ids_not_found = [], []
    for gene in query_output:
        if 'notfound' in gene:
            ids_not_found.append(gene['query'])
        else:
            ids_found.append(gene)

    sorted = np.argsort([int(x['_id']) for x in ids_found if 'notfound' not in x.keys()])
    name_to_id = {ids_found[x]['symbol'] : ids_found[x]['_id'] for x in sorted if int(ids_found[x]['_id'])
                  not in [7795, 100187828]}
    # for id in ids_not_found:
    #     id_to_name[id] = id
    with open(genes_file_path, 'w') as f:
        json.dump(name_to_id, f,indent=4, separators=(',', ': '))
    return name_to_id

args = Args()
network_graph = utils.read_network(args.network_file)
genes_ids = list(network_graph.nodes())
genes_names = load_names_from_net(genes_ids, args.genes_names_file_path)
