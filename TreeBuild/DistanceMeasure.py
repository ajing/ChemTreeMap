'''
    Measure and write the distance matrix for the whole molecule library
'''
import numpy as np

SMILE_COLUMNNAME = "Smiles"

def DistanceMatrix(smile_list):
    d_matrix = np.zeros()
    list_len = len(smile_list)
    for i in range(list_len):
        for j in range

def WriteAsPHYLIPFormat(d_matrix, ligand_name):
    list_len = len(ligand_name)
    for i in range(list_len)
    
    return filename

def GenerateDistanceFile(ligand_dict):
    smile_list = [ligand_dict[lig_name][SMILE_COLUMNNAME] for lig_name in ligand_dict.keys()]
    d_matrix   = DistanceMatrix(smile_list)
    filename   = WriteAsPHYLIPFormat(d_matrix, ligand_dict.keys())
    return filename
