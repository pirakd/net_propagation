import pandas as pd

def get_condition_function(function_name:str):
    try:
        return eval(function_name)
    except:
        assert 0, '{} is not a valid condition'.format(function_name)


def huntington_DDA(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    aliases = list(all_data[gene_with_aliases].Gene_Name)
    first_aliases = [x.split(';')[0] for x in aliases]
    all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene_Name'])

    # data = all_data[(all_data['P_Value (HD/C116)'] <= 0.05) & (all_data['Absolute Log2FC (HD/C116)'] >= 0.58)]
    prior_genes = list(all_data.Gene_Name)
    return prior_genes, all_data, all_data


def huntington_DDA_significant(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    aliases = list(all_data[gene_with_aliases].Gene_Name)
    first_aliases = [x.split(';')[0] for x in aliases]
    all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene_Name'])
    data = all_data[(all_data['P_Value (HD/C116)'] <= 0.05) & (all_data['Absolute Log2FC (HD/C116)'] >= 0.58)]
    prior_genes = list(data.Gene_Name)
    return prior_genes, data, all_data


def colorectal_cancer(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]
    all_data = all_data[~(all_data['Score'] > 300)]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    if gene_with_aliases.any():
        aliases = list(all_data[gene_with_aliases].Gene_Name)
        first_aliases = [x.split(';')[0] for x in aliases]
        all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene_Name'])

    # data = all_data[(all_data['P_Value (HD/C116)'] <= 0.05) & (all_data['Absolute Log2FC (HD/C116)'] >= 0.58)]
    prior_genes = list(all_data.Gene_Name)
    return prior_genes, all_data, all_data


def colorectal_cancer_significant(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    if gene_with_aliases.any():
        aliases = list(all_data[gene_with_aliases].Gene_Name)
        first_aliases = [x.split(';')[0] for x in aliases]
        all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene_Name'])
    data = all_data[(all_data['p_value'] <= 0.001) & ((all_data['Score'] >= 1.5) | (all_data['Score'] < -1.5)) ]
    prior_genes = list(data.Gene_Name)
    return prior_genes, data, all_data


def cov_data(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    if gene_with_aliases.any():
        gene_with_aliases = all_data.Gene_Name.str.contains(';')
        gene_with_aliases_indexes = gene_with_aliases.to_numpy().nonzero()[0]
        if len(gene_with_aliases):
            aliases_list = list(all_data[gene_with_aliases].Gene_Name)
            aliases_list_split = [x.split(';') for x in aliases_list]
            first_aliases = [x[0] for x in aliases_list_split]
            all_data = all_data.replace(aliases_list, first_aliases)

            for i, idx in enumerate(gene_with_aliases_indexes):
                entry = all_data.iloc[[idx], :]
                for alias in aliases_list_split[i][1:]:
                    entry.at[idx, 'Gene_Name'] = alias
                    all_data = all_data.append(entry)

    all_data = all_data.drop_duplicates(subset=['Gene_Name'])
    prior_genes = list(all_data.Gene_Name)
    return prior_genes, all_data, all_data

def cov_data_significant(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    if gene_with_aliases.any():
        gene_with_aliases = all_data.Gene_Name.str.contains(';')
        gene_with_aliases_indexes = gene_with_aliases.to_numpy().nonzero()[0]
        if len(gene_with_aliases):
            aliases_list = list(all_data[gene_with_aliases].Gene_Name)
            aliases_list_split = [x.split(';') for x in aliases_list]
            first_aliases = [x[0] for x in aliases_list_split]
            all_data = all_data.replace(aliases_list, first_aliases)

            for i, idx in enumerate(gene_with_aliases_indexes):
                entry = all_data.iloc[[idx], :]
                for alias in aliases_list_split[i][1:]:
                    entry.at[idx, 'Gene_Name'] = alias
                    all_data = all_data.append(entry)

    all_data = all_data.drop_duplicates(subset=['Gene_Name'])
    data = all_data[(all_data['pvalue'] <= 0.01)]
    prior_genes = list(data.Gene_Name)
    return prior_genes, data, all_data


def yeast_data(data_frame):
    all_data = data_frame[~data_frame['Score'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene_Name.str.contains(';')
    if gene_with_aliases.any():
        gene_with_aliases = all_data.Gene_Name.str.contains(';')
        gene_with_aliases_indexes = gene_with_aliases.to_numpy().nonzero()[0]
        if len(gene_with_aliases):
            aliases_list = list(all_data[gene_with_aliases].Gene_Name)
            aliases_list_split = [x.split(';') for x in aliases_list]
            first_aliases = [x[0] for x in aliases_list_split]
            all_data = all_data.replace(aliases_list, first_aliases)

            for i, idx in enumerate(gene_with_aliases_indexes):
                entry = all_data.iloc[[idx], :]
                for alias in aliases_list_split[i][1:]:
                    entry.at[idx, 'Gene_Name'] = alias
                    all_data = all_data.append(entry)

    all_data = all_data.drop_duplicates(subset=['Gene_Name'])
    prior_genes = list(all_data.Gene_Name)
    return prior_genes, all_data, all_data