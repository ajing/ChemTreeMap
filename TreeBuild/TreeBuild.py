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

    ################## Adding properties ######################
    lig_show     = SelectColumn(liganddict, INTEREST_COLUMN)
    # add pIC50, here we assume the unit for IC50 is nM
    lig_show     = ReArrangeActivity(lig_show, ACTIVITY_COLUMN)

    # add ligand efficiency
    lig_show     = AddLigEff(lig_show, liganddict)
    # add ligand efficiency
    lig_show     = AddSLogP(lig_show, liganddict)
    lig_list     = Dict2List(lig_show)

    ################## Tree structure info  ######################
    tree_meta = [{ "name":     "ECFP",
                  "metadata": """Extended Connectivity fingerprint, """
                              """implemented in <a href="http://www.rdkit.org">RDKit</a>. <br/>"""
                              """Parameters used: Radius = 4"""},
                 { "name": "AtomPair",
                   "metadata": """Atom Pairs as Molecular Features, describe in """
                         """ R.E. Carhart, D.H. Smith, R. Venkataraghavan;"""
                         """ "Atom Pairs as Molecular Features in Structure-Activity Studies:"""
                         """ Definition and Applications" JCICS 25, 64-73 (1985)."""
                         """implemented in <a href="http://www.rdkit.org">RDKit</a>. <br/>"""}]

    ################## Activity info  ######################
    activity_meta = [{ "name": "IC50",
                      "metadata" : "User input activity"},
                     { "name": "pIC50",
                       "metadata" : "This number assumes IC50 in nM unit, so please change your data or the code to make it appropriate." }]

    ################## Property info  ######################
    property_meta = [{"name": "SLogP",
                      "metadata": """SLogP, the coefficients are a measure of the difference in solubility of the compound in water and octanol. describe in """
                         """   S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)"""
                         """ R.E. Carhart, D.H. Smith, R. Venkataraghavan;"""
                         """ "Atom Pairs as Molecular Features in Structure-Activity Studies:"""},
                     {"name": "Lig_Eff",
                      "metadata": "Ligand efficiency. The value is calculated by function 1.37 * pIC50 / a_heavy"}]


    WriteJSON({"metadata": {"activityTypes": activity_meta,
                            "treeTypes": tree_meta,
                            "circleSizeTypes": property_meta,
                            "circleBorderTypes": property_meta},
               "trees": trees,
               "compounds":lig_list},
               outfile, "w")

    print "finish writing to JSON..."

    # make image file
    MakeStructuresForSmiles(liganddict)


if __name__ == "__main__":
    main()
