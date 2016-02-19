'''
    Update the old JSON file
'''

import json
import argparse

from Util import ParseLigandFile

def EachNodeRecursive(tree_json, newdict, c_min_list):
    if not "children" in tree_json:
        n_size = newdict[tree_json["orig_id"]]["SlogP"]
        if n_size < c_min_list[0]:
            c_min_list[0] = n_size
        tree_json["size"] = n_size
    else:
        for child in tree_json["children"]:
            EachNodeRecursive(child, newdict, c_min_list)


def UpdateProp(infile, propfile, outfile):
    lig_dict = ParseLigandFile(propfile)
    new_dict = dict((lig_dict[ekey]["orig_id"], lig_dict[ekey]) for ekey in lig_dict)
    filename = infile.split(".")
    #outfile = filename[0] + "_new." + filename[1]
    outfile = outfile
    outobj  = open(outfile, "w")
    inobj   = open(infile)
    jsonc   = json.load(inobj)
    c_min_list = [10]
    EachNodeRecursive(jsonc, new_dict, c_min_list)
    print c_min_list

    outobj.write(json.dumps(jsonc, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update properties')
    parser.add_argument('--jsonfile')
    parser.add_argument("--propfile")
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    infile = arg.jsonfile
    propfile = arg.propfile
    outfile  = arg.outfile
    UpdateProp(infile, propfile, outfile)
