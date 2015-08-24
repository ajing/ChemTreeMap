'''
    The main function for tree building
'''

import argparse
from DistanceMeasure import GenerateDistanceFile
from Util import ParseLigandFile, WriteJSON
from RunRapidNJ import RunRapidNJ
from RunGraphViz import WriteDotFile, RunGraphViz, Dot2JSON


def main():
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    infile = arg.infile
    outfile= arg.outfile

    liganddict   = ParseLigandFile(newfile)  # keep all ligand related information, key is ligand id, value is also a dictionray
    distant_file = GenerateDistanceFile(liganddict)
    newick_cont  = RunRapidNJ(distant_file)
    dotfilename  = WriteDotFile(newick_cont)
    outdot_file  = RunGraphViz(dotfilename)
    tree_dict    = Dot2Dict(outdotfile, liganddict)
    WriteJSON(tree_dict, outfile)

if __name__ == "__main__":
    main()
