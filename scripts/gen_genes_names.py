"""
This script creates a name-id mapping of all genes in the network and saves it afterwards
"""

import utils as utils
import numpy as np
from args import Args
import mygene
import json


def load_names_from_net(gene_ids_to_search):
    mg = mygene.MyGeneInfo()
    query_output = mg.querymany(gene_ids_to_search, scopes="entrezgene", fields=["entrezgene", "symbol"], species="human")

    ids_found, ids_not_found = [], []
    for gene in query_output:
        if 'notfound' in gene:
            ids_not_found.append(gene['query'])
        else:
            ids_found.append(gene)

    sorted = np.argsort([int(x['_id']) for x in ids_found if 'notfound' not in x.keys()])
    name_to_id = {ids_found[x]['symbol'] : ids_found[x]['_id'] for x in sorted if int(ids_found[x]['_id'])
                  not in [7795, 100187828]}
    return name_to_id


args = Args(is_create_output_folder=False)
network_graph = utils.read_network(args.network_file)
genes_ids = list(network_graph.nodes())
name_to_id = load_names_from_net(genes_ids)

# save json
with open(args.genes_names_file_path, 'w') as f:
    json.dump(name_to_id, f, indent=4, separators=(',', ': '))
