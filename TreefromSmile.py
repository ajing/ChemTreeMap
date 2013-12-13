'''
    Construct a tree from smile strings
'''

import sys
sys.path.append("../clusterVis")
from ligandGraphall import NewLigandFile, parseLigandFile, similarityMatrix, getSimilarity
from TreeConstruction import nj
from CreateGraph import MoleculeDictionary
from Dot2JSON import Dot2JSON, Root2JSON
import random

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

def TreefromSmile(infile):
    newfile = SamplingLigandFile(infile, 300, 300)
    liganddict = parseLigandFile(newfile)
    NewLigandFile(liganddict, newfile)
    smatrix  = similarityMatrix(liganddict, getSimilarity)
    moldict  = MoleculeDictionary(newfile)
    # so write dot language to file
    dotfile  = nj(smatrx, moldict, True)
    # write to JSON file
    root = Dot2JSON(dotfile)
    Root2JSON(root)

def test():
    samplefile = "Data/ligand_5_7_ppilot.txt"
    # test for Sampling
    #SamplingLigandFile(samplefile, 100, 100)
    # test for TreefromSmile
    TreefromSmile(samplefile)

if __name__ == "__main__":
    test()
