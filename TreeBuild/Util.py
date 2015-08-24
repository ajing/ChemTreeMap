'''
    Provide utility functions
'''

from csv import DictReader

def ParseLigandFile(infile):
    '''
       parse ligand file to an dictionary, key is ligand id and valud is a dictionray with properties and property values
    '''
    mol_dict = dict()
    for line in DictReader(open(infile), delimiter = "\t"):
        mol_dict[line["ligandid"]] = line
    return mol_dict


def WriteJSON(dict_obj, outfile):
    pass


if __name__ == "__main__":
    ParseLigandFile("./Data/thrombin_clean_ct.txt")

