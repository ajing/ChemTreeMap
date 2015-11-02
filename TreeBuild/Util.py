'''
    Provide utility functions
'''

from csv import DictReader
import json
import math

from Model import IDENTIFIER

def GuestByFirstLine(firstline):
    num_colnam = []
    for key in firstline:
        try:
            float(firstline[key])
            num_colnam.append(key)
        except:
            continue
    return num_colnam

def ConvertToFloat(line, colnam_list):
    for name in colnam_list:
        line[name] = round(float(line[name]), 3)
    return line

def ParseLigandFile(infile):
    '''
       parse ligand file to an dictionary, key is ligand id and valud is a dictionray with properties and property values
       This program will guess the type for each column based on the first row. The program will assume there is only two type of data: number and string.
    '''
    mol_dict = dict()
    flag = 1 # first line flag
    id_count = 0
    for line in DictReader(open(infile), delimiter = "\t"):
        if flag:
            num_colnam = GuestByFirstLine({k:v for k,v in line.items() if not k in [IDENTIFIER] })
        new_id    = "B" + str(id_count)
        id_count += 1
        mol_dict[new_id] = ConvertToFloat({k:v for k,v in line.items() if not k in [IDENTIFIER] }, num_colnam)
        mol_dict[new_id]["orig_id"] = line[IDENTIFIER]
    return mol_dict


def WriteJSON(dict_obj, outfile, write_type):
    fileobj = open(outfile, write_type)
    fileobj.write(json.dumps(dict_obj, indent=2))

def SelectColumn(lig_dict, colname):
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
    return lig_new

def ReArrangeActivity(lig_dict, colname):
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if not sk in colname}
        lig_new[k]["activities"] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
        if "IC50" in lig_new[k]["activities"]:
            # This is for nM
            #lig_new[k]["activities"]["pIC50"] = round(9 - math.log10(float(lig_new[k]["activities"]["IC50"])), 5)
            # This is for uM
            lig_new[k]["activities"]["pIC50"] = round(6 - math.log10(float(lig_new[k]["activities"]["IC50"])), 5)
    return lig_new

def Dict2List(lig_dict):
    lig_list = []
    for idx in range(len(lig_dict.keys())):
        new_entry = lig_dict["B" + str(idx)]
        new_entry["id"] = "B" + str(idx)
        lig_list.append(new_entry)
    return lig_list


if __name__ == "__main__":
    ParseLigandFile("./Data/thrombin_clean_ct.txt")
