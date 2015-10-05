'''
    Joining two molecule libraries by clustering two libraries seperately first and then join them together
'''
from Util import ParseLigandFile
from ClusterMolecules import ToFPObjRDK
from DistanceMeasure import getSimilarity

import argparse
from csv import DictWriter

def IsInList1(ligand):
    return ligand["orig_id"].startswith("CHEMBL")

def MergeTwoDicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def JoinTwoLib(infile, outfile):
    group1      = "blue"
    group2      = "red"
    ligand_dict = ParseLigandFile(infile)
    ligand_list = ligand_dict.values()
    fieldnames  = ligand_list[0].keys() + ["group"]
    lib1_list   = [MergeTwoDicts(ligand, {"group": group1}) for ligand in ligand_list if IsInList1(ligand)]
    lib2_list   = [MergeTwoDicts(ligand, {"group": group2}) for ligand in ligand_list if not IsInList1(ligand)]
    smile1_list = [[x["orig_id"], x["Canonical_Smiles"]] for x in lib1_list]
    smile2_list = [[x["orig_id"], x["Canonical_Smiles"]] for x in lib2_list]
    fp1_list     = ToFPObjRDK(smile1_list)
    fp2_list     = ToFPObjRDK(smile2_list)
    print "size of fp1 list", len(fp1_list)
    print "size of fp2 list", len(fp2_list)
    newlist = list(lib1_list)
    for idx in range(len(lib2_list)):
        fp2 = fp2_list[idx]
        overlap = False
        for idy in range(len(lib1_list)):
            fp1 = fp1_list[idy]
            if getSimilarity(fp1[1], fp2[1]) == 1:
                newlist[idy]["group"] = "black"
                overlap = True
                break
        if not overlap:
            newlist.append(lib2_list[idx])
    writer = DictWriter(open(outfile, "w"), fieldnames=fieldnames, delimiter = "\t")
    writer.writeheader()
    for item in newlist:
        writer.writerow(item)

def main():
    parser = argparse.ArgumentParser(description='join two molecule libraries.')
    parser.add_argument('--infile')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    JoinTwoLib(arg.infile, arg.outfile)


if __name__ == "__main__":
    main()
