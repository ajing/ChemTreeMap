'''
    The main function for tree building (clustering)
'''
import argparse
from rdkit.Chem import AllChem

from DistanceMeasure import GenerateDistanceFile

from RunGraphViz import WriteDotFile, SFDPonDot, Dot2Dict
from RunRapidNJ import RunRapidNJ
from Util import ParseLigandFile, WriteJSON, SelectColumn
from treebuild.Model import INTEREST_COLUMN


# ECFP6
def ECFP6(mol):
    return AllChem.GetMorganFingerprint(mol, 3)

def main():
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    infile = arg.infile
    outfile= arg.outfile

    liganddict   = ParseLigandFile(infile)  # keep all ligand related information, key is ligand id, value is also a dictionray
    #print liganddict
    print "finish parsing ligand file"
    distant_file = GenerateDistanceFile(liganddict, ECFP6)
    print "finish generating distance file..."
    newick_cont  = RunRapidNJ(distant_file)
    print "finish running RapidNJ..."
    dotfilename  = WriteDotFile(newick_cont)
    outdot_file  = SFDPonDot(dotfilename, 10)
    print "finish GraphViz..."
    lig_show     = SelectColumn(liganddict, INTEREST_COLUMN)
    # add ligand efficiency
    #lig_show     = AddLigEff(lig_show, liganddict)
    #print outdot_file
    tree_dict    = Dot2Dict(outdot_file, lig_show)
    WriteJSON(tree_dict, outfile, "w")

if __name__ == "__main__":
    main()
