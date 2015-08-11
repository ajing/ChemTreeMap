'''
    Construct a tree from affinity file
'''

import sys
sys.path.append("..")
sys.path.append("../../clusterVis")
from ParseFile import ParseAffinity
from ligandGraphall import similarityMatrix, getSimilarity
from TreeConstruction import nj, DistanceMatrix
from TreefromSmile import Matrix2JSON
from SFDPLayOut import SFDPonDot
from Dot2JSON import Dot2JSON, Root2JSON

def NewLigandFile(ligand_dict, filename):
    newfilename = filename + "_new"
    out_obj     = open( newfilename, "w" )
    in_lines    = open( filename, "r" ).readlines()
    for eachLigandID in ligand_dict:
        #print eachLigandID
        for eachline in in_lines:
            content = eachline.split("\t")
            #print "content:", content[0]
            if content[0] == eachLigandID:
                out_obj.write(eachline)
                break
    out_obj.close()
    print "finish writing to new file"

def Matrix2JSON(smatrix, affinity_dict, outfile):
    dmatrix  = DistanceMatrix(affinity_dict.keys(), 1 - smatrix)
    moldict  = dict( (key, [ 1, affinity_dict[key]]) for key in affinity_dict.keys())
    # so write dot language to file
    dotfile  = nj(dmatrix, moldict, True)
    # dotfile with sfdp layout
    # write to JSON file
    sfdp_dot  = SFDPonDot(dotfile, 10)
    print sfdp_dot
    root = Dot2JSON(sfdp_dot)
    Root2JSON(root, outfile)

def TreefromAffinity(infile):
    smile_dict, affin_dict = ParseAffinity(infile)
    NewLigandFile(smile_dict, infile)
    smatrix  = similarityMatrix(smile_dict, getSimilarity)
    Matrix2JSON(smatrix, affin_dict, "test.json")

if __name__ == "__main__":
    infile = "affinity.txt"
    TreefromAffinity(infile)
