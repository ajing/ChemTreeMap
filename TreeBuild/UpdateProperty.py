'''
    Update property of nodes
'''
import json
import argparse

from Util import ParseLigandFile
from Model import POTENCY

def UpdateProp(infile, propfile, outfile):
    lig_dict = ParseLigandFile(propfile)
    filename = infile.split(".")
    #outfile = filename[0] + "_new." + filename[1]
    outfile = outfile
    outobj  = open(outfile, "w")
    inobj   = open(infile)
    jsonc   = json.load(inobj)
    for compound in jsonc["compounds"]:
        compound["activities"][POTENCY] = lig_dict[compound["id"]][POTENCY]

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
