def get_condition_function(function_name:str):
    if function_name == 'kent_mock_no_vic_mock_24h':
        return kent_mock_no_vic_mock_24h
    elif function_name == 'kent_mock_no_vic_mock_10h':
        return kent_mock_no_vic_mock_10h


def kent_mock_no_vic_mock_24h(data_frame):
    expressed_in = 'Kent_24h-Mock_24h'
    not_expressed_in = 'VIC_24h-Mock_24h'
    pos_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                data_frame.Label == expressed_in)].Gene_Name.unique())
    neg_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                data_frame.Label == not_expressed_in)].Gene_Name.unique())
    return list(set.difference(pos_genes, neg_genes))


def kent_mock_no_vic_mock_10h(data_frame):
    expressed_in = 'Kent_10h-Mock_10h'
    not_expressed_in = 'VIC_10h-Mock_10h'
    pos_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                data_frame.Label == expressed_in)].Gene_Name.unique())
    neg_genes = set(data_frame[(data_frame.diffexpressed == True) & (
                data_frame.Label == not_expressed_in)].Gene_Name.unique())
    return list(set.difference(pos_genes, neg_genes))


