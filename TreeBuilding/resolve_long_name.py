'''
 resolve long name problem
'''

from Dot2JSON import Dot2JSON, Root2JSON, RecursiveNode2Dict
from TreeRebuilder import RewriteDot
from CreateGraph import MoleculeDictionary
from TreefromSmile import RecursiveChangeName
from SFDPLayOut import SFDPonDot
import argparse

import json


INTERESTED = ["IC50"]

def first():
    # run the sfdp command, and then run the following code
    root = Dot2JSON("test_sfdp.gv")
    Root2JSON(root, "test.json")

def ReMapName(new2old, moldict, newmap):
    moldict_new = dict()
    new2old_new = dict()
    for key in moldict.keys():
        moldict_new["B" + str(newmap[key])] = moldict[key]
    for key in new2old.keys():
        new2old_new["B" + str(newmap[key])] = new2old[key]
    return new2old_new, moldict_new

def Convert2NJmoldict(moldict):
    # convert from MoleculeDictionary from CreateGraph.py to what can work in TreeConstruction nj
    newdict = dict()
    old2new = dict()
    new2old = dict()
    i = 0
    for eachkey in moldict:
        ligandname = moldict[eachkey]["ligandid"]
        newname = "B" + str(i)
        newdict[newname] = dict()
        old2new[ligandname] = newname
        new2old[newname]    = ligandname
        i = i + 1
        if "size" in moldict[eachkey]:
            clustersize = moldict[eachkey]["size"]
            newdict[newname]["size"] = clustersize
        else:
            newdict[newname]["size"] = 1
        for dict_name in moldict[eachkey].keys():
            if dict_name in INTERESTED:
                newdict[newname][dict_name] = moldict[eachkey][dict_name]
    return newdict, old2new, new2old

def second():
    # rewrite the file
    parser = argparse.ArgumentParser(description='Simplify the graph.')
    parser.add_argument('--incsv')
    parser.add_argument('--indot')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    moldict, old2new, new2old = Convert2NJmoldict(MoleculeDictionary(arg.incsv))

    newfilename, newmapdict = RewriteDot(arg.indot)
    print newfilename
    sfdp_dot  = SFDPonDot(newfilename, 10)
    print sfdp_dot
    root = Dot2JSON(sfdp_dot)

    new2old, moldict = ReMapName(new2old, moldict, newmapdict)
    rootdict = RecursiveNode2Dict(root, moldict)
    RecursiveChangeName(rootdict, new2old)
    fileobj  = open(arg.outfile, "w")
    fileobj.write(json.dumps(rootdict, indent=2))

if __name__ == "__main__":
    second()

