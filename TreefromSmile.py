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
    bindingtypes = ["Allosteric", "Competitive"]
    allolist     = []
    complist     = []
    for line in open(infile):
        if bindingtypes[0] in line:
            allolist.append(line)
        elif bindingtypes[1] in line:
            complist.append(line)
    new_allo     = random.sample(allolist, num_allo)
    new_comp     = random.sample(complist, num_allo)
    newfile  = infile + "_sampled"
    newobj   = open(newfile, "w")
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
