'''
    Measure and write the distance matrix for the whole molecule library
'''
from Model import SMILE_COLUMNNAME, FILE_FORMAT

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import datetime

def ToFPObj(alist):
    # for alist, the first item is ligand name, the second is smile
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

def DistanceList(smile_list):
    d_list   = []
    list_len = len(smile_list)
    fp_list = ToFPObj(smile_list)
    for i in range(list_len):
        lig1 = fp_list[i]
        lig1list = []
        for j in range(0, i):
            lig2 = fp_list[j]
            sim  = getSimilarity(lig1[1], lig2[1])
            if sim:
                lig1list.append([lig2[0], 1 - sim])
        d_list.append([lig1[0], lig1list])
    return d_list

def WriteAsPHYLIPFormat(d_list):
    newfilename = datetime.datetime.now().strftime(FILE_FORMAT) + ".dist"
    fileobj  = open(newfilename, "w")
    fileobj.write( str(len(d_list)) + "\n")
    for eachrow in d_list:
        sim_values = [ "%.4f" % x[1] for x in eachrow[1]]
        line = "\t".join([eachrow[0], "\t".join(sim_values)]) + "\n"
        fileobj.write(line)
    return newfilename

def GenerateDistanceFile(ligand_dict):
    # smile_list contains ligand name and ligand smile
    smile_list = [ [lig_name, ligand_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in ligand_dict.keys()]
    d_list     = DistanceList(smile_list)
    filename   = WriteAsPHYLIPFormat(d_list)
    return filename

if __name__ == "__main__":
    from Util import ParseLigandFile
    ligand_dict = ParseLigandFile("./Data/result_clean_no0_20.txt")
    print GenerateDistanceFile(ligand_dict)
