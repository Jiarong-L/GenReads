


def Gen_degenerate_dict(fill_marker_key = ['-'],fill_marker_val = ['']):
    degenerate_dict = {
        'A':'A',
        'T':'T',
        'C':'C',
        'G':'G',
        'R':'AG',
        'Y':'CT',
        'M':'AC',
        'K':'GT',
        'S':'GC',
        'W':'AT',
        'H':'ATC',
        'B':'GTC',
        'V':'GAC',
        'D':'GAT',
        'N':'ATCG'
    }
    for key,val in zip(fill_marker_key,fill_marker_val):
        degenerate_dict[key] = val
    return degenerate_dict