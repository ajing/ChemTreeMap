'''
    Construct a tree from smile strings
'''

import sys
from ligandGraphall import NewLigandFile, parseLigandFile, similarityMatrix, getSimilarity
from TreeConstruction import nj, DistanceMatrix
from TreeRebuilder import RewriteDot
from CreateGraph import MoleculeDictionary
from Dot2JSON import Dot2JSON, RecursiveNode2Dict
import random
from SFDPLayOut import SFDPonDot
from Model import NAME_MAP_FILE
from TreeParser import GetLevelFromName
import json

import argparse

NOT_INTERESTED = ["ligandid", "Canonical_Smiles"]

def SamplingLigandFile(infile, num_allo, num_comp):
    # return new ligand filename
    newfile  = infile + "_sampled"
    newobj   = open(newfile, "w")
    bindingtypes = ["allosteric", "competitive"]
    allolist     = []
    complist     = []
    # header flag
    flag         = 1
    for line in open(infile):
        if flag:
            header = line
            newobj.write(header)
            flag = 0
            continue
        if bindingtypes[0] in line:
            allolist.append(line)
        elif bindingtypes[1] in line:
            complist.append(line)
    new_allo     = random.sample(allolist, num_allo)
    new_comp     = random.sample(complist, num_allo)
    for each in new_allo + new_comp:
        newobj.write(each)
    newobj.close()
    return newfile

def Convert2NJmoldict(moldict):
    # convert from MoleculeDictionary from CreateGraph.py to what can work in TreeConstruction nj
    newdict = dict()
    for eachkey in moldict:
        ligandname = moldict[eachkey]["ligandid"]
        newdict[ligandname] = dict()
        if "size" in moldict[eachkey]:
            clustersize = moldict[eachkey]["size"]
            newdict[ligandname]["size"] = clustersize
        else:
            newdict[ligandname]["size"] = 1
        for dict_name in moldict[eachkey].keys():
            if not dict_name in NOT_INTERESTED:
                newdict[ligandname][dict_name] = moldict[eachkey][dict_name]
    return newdict

def RecursiveChangeName(node, namedict):
    if not "children" in node.keys():
        node["name"] = namedict[node["name"]]
    else:
        if not "_" in namedict[node["name"]]:
            node["name"] = namedict[node["name"]]
        for c_node in node["children"]:
            RecursiveChangeName(c_node, namedict)

def GetRootName(namedict):
    maxlevel = 0
    rootname = ""
    for name in namedict:
        namelevel = GetLevelFromName(name)
        if maxlevel < namelevel:
            maxlevel = namelevel
            rootname = name
    return rootname

def Matrix2JSON(smatrix, liganddict, newfile, filename):
    moldict  = Convert2NJmoldict(MoleculeDictionary(newfile))
    dmatrix  = DistanceMatrix(liganddict.keys(), 1 - smatrix)
    # so write dot language to file
    dotfile  = nj(dmatrix, moldict, True)
    # dotfile with sfdp layout
    # write to JSON file
    newdotfile, newmapdict = RewriteDot(dotfile)

    rootname = GetRootName(newmapdict)

    sfdp_dot  = SFDPonDot(newdotfile, 10)
    root = Dot2JSON(sfdp_dot, newmapdict[rootname])
    # update keys for moldict
    moldict_new = dict()
    for key in moldict.keys():
        moldict_new[newmapdict[key]] = moldict[key]
    rootdict    = RecursiveNode2Dict(root, moldict_new)
    inv_mapdict = { v : k for k, v in newmapdict.items()}
    RecursiveChangeName(rootdict, inv_mapdict )

    fileobj  = open(filename, "w")
    fileobj.write(json.dumps(rootdict, indent=2))

def TreefromSmile(infile, outfile, sample = False):
    if sample:
        newfile = SamplingLigandFile(infile, 10, 10)
    else:
        newfile    = infile
    liganddict = parseLigandFile(newfile)
    NewLigandFile(liganddict, newfile)
    smatrix  = similarityMatrix(liganddict, getSimilarity)
    Matrix2JSON(smatrix, liganddict, newfile, outfile)

def test():
    #parse argument
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')
    parser.add_argument('--outfile')
    arg = parser.parse_args()
    print arg.infile
    print arg.outfile

    samplefile = "Data/ligand_5_7_ppilot.txt"
    samplefile = "Data/result_clean.txt"
    samplefile = "Data/result_clean_no0.txt"
    samplefile = "Data/result_clean_no0_20.txt"
    #samplefile = "Data/result_clean_no0_100.txt"
    # test for Sampling
    #SamplingLigandFile(samplefile, 100, 100)
    # test for TreefromSmile
    TreefromSmile(arg.infile, arg.outfile, False)

if __name__ == "__main__":
    test()
