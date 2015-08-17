'''
    Construct a tree from smile strings
'''

import sys
from ligandGraphall import NewLigandFile, parseLigandFile, similarityMatrix, getSimilarity
from TreeConstruction import nj, DistanceMatrix
from CreateGraph import MoleculeDictionary
from Dot2JSON import Dot2JSON, RecursiveNode2Dict
import random
from SFDPLayOut import SFDPonDot
from Model import NAME_MAP_FILE
import json

INTERESTED = ["IC50"]

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

def RecursiveChangeName(node, namedict):
    if not "children" in node.keys():
        node["name"] = namedict[node["name"]]
    else:
        for c_node in node["children"]:
            RecursiveChangeName(c_node, namedict)

def Matrix2JSON(smatrix, liganddict, newfile, filename):
    moldict, old2new, new2old  = Convert2NJmoldict(MoleculeDictionary(newfile))
    dmatrix  = DistanceMatrix([old2new[x] for x in liganddict.keys()], 1 - smatrix)
    # so write dot language to file
    dotfile  = nj(dmatrix, moldict, True)
    # dotfile with sfdp layout
    # write to JSON file
    sfdp_dot  = SFDPonDot(dotfile, 10)
    root = Dot2JSON(sfdp_dot)
    rootdict = RecursiveNode2Dict(root, moldict)
    RecursiveChangeName(rootdict, new2old)
    fileobj  = open(filename, "w")
    fileobj.write(json.dumps(rootdict, indent=2))

def TreefromSmile(infile, sample = False):
    if sample:
        newfile = SamplingLigandFile(infile, 10, 10)
    else:
        newfile    = infile
    liganddict = parseLigandFile(newfile)
    NewLigandFile(liganddict, newfile)
    smatrix  = similarityMatrix(liganddict, getSimilarity)
    Matrix2JSON(smatrix, liganddict, newfile, "test.json")

def test():
    samplefile = "Data/ligand_5_7_ppilot.txt"
    samplefile = "Data/result_clean.txt"
    samplefile = "Data/result_clean_no0.txt"
    #samplefile = "Data/result_clean_no0_100.txt"
    # test for Sampling
    #SamplingLigandFile(samplefile, 100, 100)
    # test for TreefromSmile
    TreefromSmile(samplefile, False)

if __name__ == "__main__":
    test()
