'''
    Provide utility functions
'''

from csv import DictReader
import json

def ParseLigandFile(infile):
    '''
       parse ligand file to an dictionary, key is ligand id and valud is a dictionray with properties and property values
    '''
    mol_dict = dict()
    for line in DictReader(open(infile), delimiter = "\t"):
        mol_dict[line["ligandid"]] = line
    return mol_dict


def WriteJSON(dict_obj, outfile):
    fileobj = open(outfile, "w")
    fileobj.write(json.dumps(dict_obj, indent=2))

def SelectColumn(lig_dict, colname):
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
    return lig_new


if __name__ == "__main__":
    ParseLigandFile("./Data/thrombin_clean_ct.txt")

