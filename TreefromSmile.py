'''
    Construct a tree from smile strings
'''

import sys
sys.path.append("../clusterVis")
from ligandGraphall import NewLigandFile, parseLigandFile, similarityMatrix, getSimilarity
from TreeConstruction import nj
from CreateGraph import MoleculeDictionary
from Dot2JSON import Dot2JSON, Root2JSON

def SamplingLigandFile(infile, num_allo, num_comp):
    # return new ligand filename
    pass

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
