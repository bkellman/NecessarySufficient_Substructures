"""
@author: jtsorren
Conditional Probability Matrix (CPM)
"""

import pandas as pd

def conditional_matrix(glycan_id, substructure_info):
    """
    Generates the substructure conditional probability matrix for 
    an input glycan. Input glycans must be an entry in the original 
    substructure information. Returns a conditional probability matrix
    where columns are conditions; P_ij is P(i|j). *appearance is binary*
    """
    substructure_db = pd.DataFrame(data=substructure_info)
    substructures = substructure_db[substructure_db[glycan_id] != 0].index
    cond_prob_matrix = pd.DataFrame(columns=substructures, index=substructures)
    for index in cond_prob_matrix.index:
        for column in cond_prob_matrix.columns:
            b = substructure_db.T
            b = b[b[column] != 0]
            a = b[b[index] != 0]
            cond_prob_matrix[index][column] = len(a) / len(b)
    return cond_prob_matrix