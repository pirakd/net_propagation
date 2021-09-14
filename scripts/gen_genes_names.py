"""
This script creates a name-id mapping of all genes in the network and saves it afterwards
"""
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
import utils as utils
import numpy as np
from args import Args
import mygene
import json
import pandas as pd


def load_names_from_net(gene_ids_to_search, name_type):
    mg = mygene.MyGeneInfo()
    query_output = mg.querymany(gene_ids_to_search, scopes="entrezgene", fields=["entrezgene",  name_type], species=4932)

    ids_found, ids_not_found = [], []
    for g,gene in enumerate(query_output):
        a = len(ids_found)
        if 'notfound' in gene or name_type not in gene:
            ids_not_found.append(gene['query'])
        else:
            ids_found.append(gene)

    sorted = np.argsort([int(x['_id']) for x in ids_found if 'notfound' not in x.keys()])
    if name_type == 'uniprot':
        name_to_id = dict()
        for x in sorted:
            if 'Swiss-Prot' in ids_found[x][name_type]:
                if not isinstance(ids_found[x][name_type]['Swiss-Prot'], list):
                    name_to_id[ids_found[x][name_type]['Swiss-Prot']] =  ids_found[x]['_id']
                else:
                    name_to_id[ids_found[x][name_type]['Swiss-Prot'][0]] =  ids_found[x]['_id']
    else:
        name_to_id = {ids_found[x][name_type] : ids_found[x]['_id'] for x in sorted}
    return name_to_id

args = Args(is_create_output_folder=False)
name_type = 'symbol'
network_graph = utils.read_network(args.network_file)
genes_ids = list(network_graph.nodes())
name_to_id = load_names_from_net(genes_ids, name_type)


# save json
with open(args.genes_names_file_path, 'w') as f:
    json.dump(name_to_id, f, indent=4, separators=(',', ': '))
