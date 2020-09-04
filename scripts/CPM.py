"""
@author: jtsorren
Conditional Probability Matrix (CPM)
"""

import pandas as pd
import itertools

def substr_prop(sub,mx,b):
    return (((b[sub]*b[mx])>0).sum())/(b[sub]>0).sum()


def conditional_vect(glycan_id, substructure_info):
    """
    Generates the substructure conditional probability matrix for 
    an input glycan. Input glycans must be an entry in the original 
    substructure information. Returns a conditional probability matrix
    where columns are conditions; P_ij is P(i|j). *appearance is binary*
    """
    
    substructure_db = pd.DataFrame(data=substructure_info)
    substructures = substructure_db[substructure_db[glycan_id] != 0].index
    
    b = substructure_db.T
    mx = substructures[-1]
    
    cond = pd.DataFrame( {glycan_id : [substr_prop(sub,mx,b) for sub in substructures]} ,
                        index = substructures)
    return cond
    
    
def substr_cbx_prop(sub,mx,b):
    return (((b[sub[0]]*b[sub[1]]*b[mx])>0).sum())/(b[sub]>0).sum()


def conditional_cbx_vect(glycan_id, substructure_info):
    """
    Generates the substructure conditional probability matrix for 
    an input glycan. Input glycans must be an entry in the original 
    substructure information. Returns a conditional probability matrix
    where columns are conditions; P_ij is P(i|j). *appearance is binary*
    """
    
    substructure_db = pd.DataFrame(data=substructure_info)
    substructures = substructure_db[substructure_db[glycan_id] != 0].index
    substructures = list(itertools.combinations(substructures,2))
    
    b = substructure_db.T
    mx = substructures[-1]
    
    cond = pd.DataFrame( {glycan_id : [substr_cbx_prop(sub,mx,b) for sub in substructures]} ,
                        index = substructures)
    return cond
    
  




########################

def conditional_matrix(glycan_id, substructure_info,largest_only=True):
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
            cond_prob_matrix[column][index] = len(a) / len(b)
    if largest_only:
        return cond_prob_matrix.loc[max(substructures)]
    else:
        return cond_prob_matrix
    

def conditional_matrix_cbx(glycan_id, substructure_info,largest_only=True):
    """
    Generates the substructure conditional probability matrix for 
    an input glycan. Input glycans must be an entry in the original 
    substructure information. Returns a conditional probability matrix
    where columns are conditions; P_ij is P(i|j). *appearance is binary*
    """
    substructure_db = pd.DataFrame(data=substructure_info)
    substructures = substructure_db[substructure_db[glycan_id] != 0].index
    cond_prob_matrix = pd.DataFrame(columns=list(itertools.combinations(substructures,2)), index=substructures)
    for index in cond_prob_matrix.index:
        for column in cond_prob_matrix.columns:
            b = substructure_db.T
            if not index in column:
                b = b[b[column[0]] != 0]
                b = b[b[column[1]] != 0]
                a = b[b[index] != 0]
                cond_prob_matrix[column][index] = len(a) / len(b)
    if largest_only:
        return cond_prob_matrix.loc[max(substructures)]
    else:
        return cond_prob_matrix
