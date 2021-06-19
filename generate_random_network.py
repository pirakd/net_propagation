import joblib
import old_utils as utils
from old_utils import save_file
import networkx as nx
verbosity = True
import os
from datetime import datetime

time = datetime.today().strftime('%d_%m_%Y__%H_%M_%S_%f')
output_folder_name = 'random_networks'
os.makedirs(output_folder_name, exist_ok=True)
number_of_networks = 100

#Read the h_sapiens network
network_graph = utils.read_network("data/H_sapiens.net")

E = network_graph.number_of_edges()
Q = 10
random_networks = {}
networks_without_propagation = {}

for i in range(number_of_networks):
    network_name = datetime.today().strftime('%d_%m_%Y__%H_%M_%S_%f.pl')
    print('processing network {}'.format(i))
    H = network_graph.copy()
    nx.swap.double_edge_swap(H, nswap=Q*E, max_tries=Q*E*2)
    network_file_path = os.path.join(output_folder_name, network_name)
    save_file(H, network_file_path)
