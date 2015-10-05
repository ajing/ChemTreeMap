'''
    Get the size of a fingerprint
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

import numpy as np

# ECFP6
def ECFP4(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

#smi = "c(cc(n1C)C)(c1)c(csc2\N=C(\N)/N)n2"
smi = "[nH]1c(ccc(C(=O)O)c2)c2cc1C(=O)OCC"
mol = Chem.MolFromSmiles(smi)
fp1 = ECFP4(mol)
fp2 = Chem.RDKFingerprint(mol)

fpa = np.zeros((1,len(fp2)),np.double)
DataStructs.ConvertToNumpyArray(fp2, fpa)
print len(fp2)
print dir(fp1)
print fp1.GetLength()
