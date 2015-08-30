'''
    Measure and write the distance matrix for the whole molecule library
'''
from Model import SMILE_COLUMNNAME, FILE_FORMAT, POTENCY

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import datetime
import math

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

def WriteAsPHYLIPFormat(smile_list):
    fp_list = ToFPObj(smile_list)
    print "finish parsing smile list"
    list_len = len(fp_list)

    newfilename = datetime.datetime.now().strftime(FILE_FORMAT) + ".dist"
    fileobj  = open(newfilename, "w")
    fileobj.write( str(list_len) + "\n")

    for i in range(list_len):
        lig1 = fp_list[i]
        lig1list = []
        for j in range(list_len):
            lig2 = fp_list[j]
            sim  = getSimilarity(lig1[1], lig2[1])
            lig1list.append([lig2[0], 1 - sim])

        sim_values = [ "%.4f" % x[1] for x in lig1list]
        line = "\t".join([lig1[0], "\t".join(sim_values)]) + "\n"
        fileobj.write(line)

    fileobj.close()

    return newfilename

def GenerateDistanceFile(ligand_dict):
    # smile_list contains ligand name and ligand smile
    smile_list = [ [lig_name, ligand_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in ligand_dict.keys()]
    print "finish smile list"
    filename   = WriteAsPHYLIPFormat(smile_list)
    print "finish writing phyli file"
    return filename

def LigEff(smile, ic50):
    m = Chem.MolFromSmiles(smile)
    num_heavy = m.GetNumHeavyAtoms()
    return 1.37 * (9 - math.log10(ic50)) / num_heavy

def AddLigEff(lig_show, liganddict):
    for ligname in lig_show:
        lig_show[ligname]["lig_eff"] = LigEff(liganddict[ligname][SMILE_COLUMNNAME], liganddict[ligname][POTENCY])
    return lig_show

if __name__ == "__main__":
    from Util import ParseLigandFile
    ligand_dict = ParseLigandFile("./Data/result_clean_no0_20.txt")
    print GenerateDistanceFile(ligand_dict)
