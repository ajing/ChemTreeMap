#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

import argparse
import os

from DistanceMeasure import GenerateDistanceFile, AddLigEff, AddSLogP
from model import INTEREST_COLUMN, ACTIVITY_COLUMN

from RunGraphViz import WriteDotFile, SFDPonDot, Dot2Dict
from RunRapidNJ import RunRapidNJ
from util import ParseLigandFile, WriteJSON, SelectColumn, ReArrangeActivity, Dict2List
from treebuild.MakeStructuresForSmiles import MakeStructuresForSmiles

# two kinds of fingerprint
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem



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
DIRNAME = os.path.dirname(os.path.abspath(__file__))

from types import FingerPrintType


class TreeBuild:
    """
    There are a few assumptions for the input file:
        1. potency unit is nM
    """
    def __init__(self, inputfile, fps, activities, other_properties):
        lig_dict = self.parse_lig_file(inputfile)
        trees = dict()
        for fp in fps:
            assert isinstance(fp, FingerPrintType)
            trees[fp.name] = self.build_single_tree(lig_dict, fp)
        metadata = dict()
        metadata["activityType"] = [act.to_dict() for act in activities]
        metadata["treeTypes"] = [fp.to_dict() for fp in fps]
        metadata["circleSizeTypes"] = [prop.to_dict() for prop in other_properties]
        metadata["circleBorderTypes"] = [prop.to_dict() for prop in other_properties]

        comp_info = self.add_properties()
        final_dict = {"metadata": metadata, "trees": trees, compounds: comp_info}

    def build_single_tree(self, lig_dict, fp):
        distfile = self.gen_dist_file(lig_dict, fp)
        newick_f = self.run_rapidnj(distfile)
        dot_inf = self.write_dotfile(newick_f)
        dot_out = self.sfdp_dot(dot_inf, 10)
        dot_dict = self.dot2dict(dot_out, None)
        return dot_dict



    def gen_dist_file(self, liganddict, fp):
        pass



if __name__ == "__main__":
    main()
