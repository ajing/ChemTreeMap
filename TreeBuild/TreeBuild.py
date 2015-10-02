'''
    The main function for tree building
'''

import argparse
from DistanceMeasure import GenerateDistanceFile, AddLigEff, AddSLogP
from Util import ParseLigandFile, WriteJSON, SelectColumn, ReArrangeActivity, Dict2List
from RunRapidNJ import RunRapidNJ
from RunGraphViz import WriteDotFile, SFDPonDot, Dot2Dict
from MakeStructuresForSmiles import MakeStructuresForSmiles

from Model import INTEREST_COLUMN, ACTIVITY_COLUMN

# two kinds of fingerprint
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem

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
    ################## for ECFP 6 ######################
    print "finish parsing ligand file"
    distant_file = GenerateDistanceFile(liganddict, ECFP6)
    print "finish generating distance file..."
    newick_cont  = RunRapidNJ(distant_file)
    print "finish running RapidNJ..."
    dotfilename  = WriteDotFile(newick_cont)
    outdot_file  = SFDPonDot(dotfilename, 10)
    print "finish GraphViz..."
    tree_dict_ecfp = Dot2Dict(outdot_file, None)
    ################## for Atom pair ######################
    print "finish parsing ligand file"
    distant_file = GenerateDistanceFile(liganddict, Pairs.GetAtomPairFingerprint)
    print "finish generating distance file..."
    newick_cont  = RunRapidNJ(distant_file)
    print "finish running RapidNJ..."
    dotfilename  = WriteDotFile(newick_cont)
    outdot_file  = SFDPonDot(dotfilename, 10)
    print "finish GraphViz..."
    tree_dict_atom = Dot2Dict(outdot_file, None)

    lig_show     = SelectColumn(liganddict, INTEREST_COLUMN)
    lig_show     = ReArrangeActivity(lig_show, ACTIVITY_COLUMN)

    # add ligand efficiency
    lig_show     = AddLigEff(lig_show, liganddict)
    # add ligand efficiency
    lig_show     = AddSLogP(lig_show, liganddict)
    lig_list     = Dict2List(lig_show)

    trees = {"ECFP": tree_dict_ecfp, "AtomPair": tree_dict_atom}

    if not "pIC50" in ACTIVITY_COLUMN:
        ACTIVITY_COLUMN.append("pIC50")

    WriteJSON({"metadata": {"activityTypes": [{"name": x, "metadata": "nothing"} for x in ACTIVITY_COLUMN], "treeTypes": [{"name":x, "metadata": "nothing"} for x in trees.keys()], "circleSizeTypes": [{"name": x, "metadata": "nothing"} for x in lig_list[0]["properties"].keys()], "circleBorderTypes": [{"name": x, "metadata": "nothing"} for x in lig_list[0]["properties"].keys()]}, "trees": trees, "compounds":lig_list}, outfile, "w")
    print "finish writing to JSON..."

    # make image file
    MakeStructuresForSmiles(liganddict)

    #print outdot_file
    # the current output structure is
    # {trees:{ECFP:root, AtomPair:root}, compounds:[{id:, name:, activities: {IC50:}, properties:},]}

if __name__ == "__main__":
    main()
