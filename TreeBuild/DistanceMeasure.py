'''
    Measure and write the distance matrix for the whole molecule library
'''
import numpy as np
from Model import SMILE_COLUMNNAME

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def ToFPObj(alist):
    newlist = []
    for each in alist:
        smile = each[1]
        m = Chem.MolFromSmiles(smile)
        if m is None:
            continue
        fp = AllChem.GetMorganFingerprint(m, 3)
        newlist.append([each[0], fp])
    return newlist

def getSimilarity(fp1, fp2):
    # generate similarity score for two smiles strings
    if (fp1 is None or fp2 is None):
        return
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def DistanceMatrix(smile_list):
    d_matrix = np.zeros()
    list_len = len(smile_list)
    smile_fp_list = ToFPObj(smile_list)
    for i in range(list_len):
        for j in range(0, i):
            pass

def WriteAsPHYLIPFormat(d_matrix, ligand_name):
    list_len = len(ligand_name)
    for i in range(list_len):
        for j in range(i, list_len):
            pass

    return filename

def GenerateDistanceFile(ligand_dict):
    smile_list = [ligand_dict[lig_name][SMILE_COLUMNNAME] for lig_name in ligand_dict.keys()]
    d_matrix   = DistanceMatrix(smile_list)
    filename   = WriteAsPHYLIPFormat(d_matrix, ligand_dict.keys())
    return filename

if __name__ == "__main__":
    from Util import ParseLigandFile
    ligand_dict = ParseLigandFile("./Data/result_clean_no0_20.txt")
    GenerateDistanceFile(ligand_dict)
