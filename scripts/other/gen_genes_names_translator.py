"""
This script creates a name-id mapping of all genes in the network and saves it afterwards
"""
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.dirname(path.realpath(__file__)))))
import utils as utils
from args import Args, HAP40

import json
from gene_name_translator.gene_translator import GeneTranslator
gs = GeneTranslator()
# gs.generate_dictionaries()
gs.load_dictionary()
args = Args(is_create_output_folder=False)
name_type = 'symbol'
network_graph = utils.read_network(args.network_file_path)
genes_ids = list(network_graph.nodes())
id_to_name = gs.translate(genes_ids, 'entrez_id', 'symbol')
name_to_id = {xx: x for x, xx in id_to_name.items()}

# save json
with open(args.genes_names_file_path, 'w') as f:
    json.dump(name_to_id, f, indent=4, separators=(',', ': '))
