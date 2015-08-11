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
        if "size" in moldict[eachkey]:
            ligandname = moldict[eachkey]["ligandid"]
            ligandtype = moldict[eachkey]["typeofbinding"]
            clustersize = moldict[eachkey]["size"]
            newdict[ligandname] = [clustersize, ligandtype]
        else:
            ligandname = moldict[eachkey]["ligandid"]
            ligandtype = moldict[eachkey]["typeofbinding"]
            newdict[ligandname] = [1, ligandtype]
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

def TreefromSmile(infile):
    newfile = SamplingLigandFile(infile, 50, 50)
    liganddict = parseLigandFile(newfile)
    NewLigandFile(liganddict, newfile)
    smatrix  = similarityMatrix(liganddict, getSimilarity)
    Matrix2JSON(smatrix, liganddict, newfile, "test.json")

def test():
    samplefile = "Data/ligand_5_7_ppilot.txt"
    # test for Sampling
    #SamplingLigandFile(samplefile, 100, 100)
    # test for TreefromSmile
    TreefromSmile(samplefile)

if __name__ == "__main__":
    test()
