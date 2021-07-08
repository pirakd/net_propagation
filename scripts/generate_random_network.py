from utils import save_file, read_network, read_network_create_graph
import networkx as nx
verbosity = True
import os
from datetime import datetime

time = datetime.today().strftime('%d_%m_%Y__%H_%M_%S_%f')
output_folder_name = 'random_networks'
os.makedirs(output_folder_name, exist_ok=True)
number_of_networks = 100

#Read the h_sapiens network
network_graph = read_network("data/H_sapiens.net")
#Use this line for htt network
#network_graph = read_network_create_graph(r'C:\Users\USER\Documents\MSc\net_propagation\net_propagation-master\net_propagation-master\data\nn1_network.xlsx', 'NN1.0ext')

E = network_graph.number_of_edges()
Q = 10
random_networks = {}
networks_without_propagation = {}

for i in range(number_of_networks):
    network_name = datetime.today().strftime('%d_%m_%Y__%H_%M_%S_%f.pl')
    print('processing network {}'.format(i))
    H = network_graph.copy()
    nx.swap.double_edge_swap(H, nswap=Q*E, max_tries=Q*E*2)
    sparse_mat = nx.to_scipy_sparse_matrix(H, weight=2)
    network_file_path = os.path.join(output_folder_name, network_name)
    save_file(sparse_mat, network_file_path)
