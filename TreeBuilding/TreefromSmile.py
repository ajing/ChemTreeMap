'''
    Construct a tree from smile strings
'''

import sys
from ligandGraphall import NewLigandFile, parseLigandFile, similarityMatrix, getSimilarity
from TreeConstruction import nj, DistanceMatrix
from CreateGraph import MoleculeDictionary
from Dot2JSON import Dot2JSON, Root2JSON
import random
from SFDPLayOut import SFDPonDot

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
        if "size" in moldict[eachkey]:
            clustersize = moldict[eachkey]["size"]
            newdict[ligandname] = [clustersize]
        else:
            newdict[ligandname] = [1]
        for dict_name in moldict[eachkey].keys():
            if dict_name != "ligandid":
                newdict[ligandname].append(moldict[eachkey][dict_name])
    return newdict

def Matrix2JSON(smatrix, liganddict, newfile, filename):
    dmatrix  = DistanceMatrix(liganddict.keys(), 1 - smatrix)
    moldict  = Convert2NJmoldict(MoleculeDictionary(newfile))
    # so write dot language to file
    dotfile  = nj(dmatrix, moldict, True)
    # dotfile with sfdp layout
    # write to JSON file
    sfdp_dot  = SFDPonDot(dotfile, 10)
    root = Dot2JSON(sfdp_dot)
    Root2JSON(root, filename)

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
    samplefile = "Data/result_clean_10.txt"
    # test for Sampling
    #SamplingLigandFile(samplefile, 100, 100)
    # test for TreefromSmile
    TreefromSmile(samplefile, False)

if __name__ == "__main__":
    test()
