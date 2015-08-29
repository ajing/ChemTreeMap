'''
    The main function for tree building
'''

import argparse
from DistanceMeasure import GenerateDistanceFile
from Util import ParseLigandFile, WriteJSON, SelectColumn
from RunRapidNJ import RunRapidNJ
from RunGraphViz import WriteDotFile, SFDPonDot, Dot2Dict

from Model import INTEREST_COLUMN

def main():
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    infile = arg.infile
    outfile= arg.outfile

    liganddict   = ParseLigandFile(infile)  # keep all ligand related information, key is ligand id, value is also a dictionray
    distant_file = GenerateDistanceFile(liganddict)
    print "finish generating distance file..."
    newick_cont  = RunRapidNJ(distant_file)
    print "finish running RapidNJ..."
    dotfilename  = WriteDotFile(newick_cont)
    outdot_file  = SFDPonDot(dotfilename, 10)
    print "finish GraphViz..."
    lig_show     = SelectColumn(liganddict, INTEREST_COLUMN)
    tree_dict    = Dot2Dict(outdot_file, lig_show)
    WriteJSON(tree_dict, outfile)

if __name__ == "__main__":
    main()
