def get_condition_function(function_name:str):
    if function_name == 'kent_vic_10h':
        return kent_vic_10h
    elif function_name == 'kent_vic_24h':
        return kent_vic_24h
    elif function_name == 'huntington_DIA':
        return huntington_DIA
    elif function_name == 'huntington_DDA':
        return huntington_DDA
    elif function_name == 'huntington_DDA_significant':
        return huntington_DDA_significant
    elif function_name == 'colorectal_cancer':
        return colorectal_cancer
    else:
        assert 0, 'function_name is not a valid condition'

def kent_vic_10h(data_frame):
    expressed_in = 'Kent_10h-VIC_10h'
    data = data_frame[data_frame.Label == expressed_in]
    data = data.drop_duplicates(subset=['Gene_Name'])
    genes = list(data.Gene_Name)
    return genes, data, data


def kent_vic_24h(data_frame):
    expressed_in = 'Kent_24h-VIC_24h'
    data = data_frame[data_frame.Label == expressed_in]
    data = data.drop_duplicates(subset=['Gene_Name'])
    genes = list(data.Gene_Name)
    return genes, data, data


def huntington_DIA(data_frame):
    all_data = data_frame[~data_frame.Qvalue.isnull()]
    all_data = all_data[~all_data.Gene.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene.str.contains(';')
    aliases = list(all_data[gene_with_aliases].Gene)
    first_aliases = [x.split(';')[0] for x in aliases]
    all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene'])
    prior_genes = list(all_data.Gene)
    # data = all_data[(all_data.Qvalue <= 0.05) & (all_data['Absolute AVG Log2 Ratio'] >= 0.58)]
    # prior_genes = list(data.Gene)
    return prior_genes, all_data, all_data


def huntington_DDA(data_frame):
    all_data = data_frame[~data_frame['Log2FC (HD/C116)'].isnull()]
    all_data = all_data[~all_data.Gene.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene.str.contains(';')
    aliases = list(all_data[gene_with_aliases].Gene)
    first_aliases = [x.split(';')[0] for x in aliases]
    all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene'])

    # data = all_data[(all_data['P_Value (HD/C116)'] <= 0.05) & (all_data['Absolute Log2FC (HD/C116)'] >= 0.58)]
    prior_genes = list(all_data.Gene)
    return prior_genes, all_data, all_data

def huntington_DDA_significant(data_frame):
    all_data = data_frame[~data_frame['Log2FC (HD/C116)'].isnull()]
    all_data = all_data[~all_data.Gene.isnull()]

    # keep only first name of genes with aliases
    gene_with_aliases = all_data.Gene.str.contains(';')
    aliases = list(all_data[gene_with_aliases].Gene)
    first_aliases = [x.split(';')[0] for x in aliases]
    all_data = all_data.replace(aliases, first_aliases)
    all_data = all_data.drop_duplicates(subset=['Gene'])

    data = all_data[(all_data['P_Value (HD/C116)'] <= 0.05) & (all_data['Absolute Log2FC (HD/C116)'] >= 0.58)]
    prior_genes = list(data.Gene)
    return prior_genes, data, all_data


def colorectal_cancer(data_frame):
    all_data = data_frame[~data_frame['log2FC'].isnull()]
    all_data = all_data[~all_data.Gene_Name.isnull()]

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