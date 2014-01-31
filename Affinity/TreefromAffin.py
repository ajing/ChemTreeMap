'''
    Construct a tree from affinity file
'''

import sys
sys.path.append("..")
sys.path.append("../../clusterVis")
from ParseFile import ParseAffinity
from ligandGraphall import similarityMatrix, getSimilarity
from TreefromSmile import Matrix2JSON

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

def TreefromAffinity(infile):
    smile_dict, affin_dict = ParseAffinity(infile)
    NewLigandFile(smile_dict, infile)
    smatrix  = similarityMatrix(smile_dict, getSimilarity)
    print smatrix

if __name__ == "__main__":
    infile = "affinity.txt"
    TreefromAffinity(infile)
