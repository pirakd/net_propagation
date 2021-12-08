import numpy as np
import re
import utils

def get_statistic_score(data_mat, classes, type='s2n'):
    pos = data_mat[:, classes == 1]
    neg = data_mat[:, classes == 0]
    pos_mean = np.mean(pos,axis=1)
    neg_mean = np.mean(neg, axis=1)
    pos_std = np.std(pos, axis=1)
    neg_std = np.std(neg, axis=1)

    if type == 's2n' or type == 'signal_to_noise':
        return (pos_mean - neg_mean) / (pos_std + neg_std)
    elif type == 't_test':
        return (pos_mean - neg_mean) / np.sqrt(pos_std**2/pos.shape[1] + neg_std**2/neg.shape[1])

def fetch_experiment_data(data_file_path, labels_file_path, args):
    genes_names = list()
    raw_data = list()
    print('Loading data'),
    with open(data_file_path, encoding='utf-8') as f:
        for l_idx, line in enumerate(f):
            if l_idx == 0:  # get protein names
                regions_names = re.split(r'[\t]', line.strip())[2:]
            else:  # get values and regions names
                split_line = re.split(r'[ \t]', line.strip())
                if split_line[0] != 'NA':
                    raw_data.append(split_line[2:])
                    genes_names.append(split_line[0])

    symbol_to_id = utils.convert_symbols_to_ids(genes_names, args.genes_names_file_path)
    id_to_symbol = {xx: x for x, xx in symbol_to_id.items()}

    filtered_raw_data = []
    filtered_genes_names = []
    unique_genes = dict()
    for i in range(len(genes_names)):
        if genes_names[i] in symbol_to_id and genes_names[i] not in unique_genes:
            filtered_genes_names.append(genes_names[i])
            filtered_raw_data.append(raw_data[i])
            unique_genes[genes_names[i]] = 1

    id_to_idx = {symbol_to_id[xx]: x for x, xx in enumerate(filtered_genes_names)}
    class_names_to_idx = {}
    with open(labels_file_path, encoding='utf-8') as f:
        for l_idx, line in enumerate(f):
            if l_idx == 0:
                pass
            elif l_idx == 1:
                split_line = line.strip().split(' ')
                class_names_to_idx[split_line[1]] = 1
                class_names_to_idx[split_line[2]] = 0
            else:
                split_line = line.strip().split('\t')
                sample_classes = [class_names_to_idx[x] for x in split_line]
    return np.array(filtered_raw_data, dtype=np.float), np.array(sample_classes, dtype=int), filtered_genes_names, symbol_to_id, id_to_idx

def change_expression(raw_data, labels, genes_idx_to_change, change_amount):
    pos_samples = np.nonzero(labels == 1)[0]
    samples_std = np.std(raw_data[genes_idx_to_change], axis=1)
    changed_elements_indexes = np.meshgrid(genes_idx_to_change, pos_samples, indexing='ij')
    raw_data[genes_idx_to_change] -= change_amount * samples_std[:, np.newaxis]
    raw_data[tuple(changed_elements_indexes)] += (2 * change_amount) * samples_std[:, np.newaxis]

    return raw_data

def shuffle_data(data_table, labels):
    n_samples = data_table.shape[1]
    suffeled_indexes = np.random.choice(n_samples, replace=True, size=n_samples)
    shuffled_labels = labels[suffeled_indexes]
    shuffled_data = data_table[:, suffeled_indexes]

    return shuffled_data, shuffled_labels