#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

# building tree for allosteric and competitive compounds
# assume the file is in the format
# ligandid[tab]smile[tab]allosteric

from treebuild import TreeBuild
from treebuild import DEFAULT_FINGERPRINT_TYPES
from sklearn.cluster import KMeans

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def SMILE2Matrix(smile_list):
    # To ECFP6
    def ToECFP(smile):
        mol = Chem.MolFromSmiles(smile)
        return AllChem.GetMorganFingerprint(mol, 3)

    fps = map(ToECFP, smile_list)

    np_fps = []
    for fp in fps:
      arr = numpy.zeros((1,))
      DataStructs.ConvertToNumpyArray(fp, arr)
      np_fps.append(arr)

    print np_fps



def LigandClusteringByClass(lig_dict, class_col = "allosteric", num_clusters = {"allo": 5, "comp" : 3}):
    """
    Ligand clustering by ligand class

    :param lig_dict: all ligand information
    :param class_col: the column name with ligand class
    :return: new lig_dict with cluster size
    """
    smile_list = [ [lig_name, lig_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in lig_dict.keys()]

    # to fp matrix
    fp_mat = SMILE2Matrix(smile_list)
    print fp_mat.shape
    print fp_mat
    return fp_mat

    # KMean clustering
    # allo_c <- KMeans(n_clusters = num_clusters["allo"]).fit()
    #
    # print "finish parsing smile list"
    # list_len = len(fp_list)
    #
    # newfilename = datetime.datetime.now().strftime(FILE_FORMAT) + ".dist"
    # fileobj  = open(newfilename, "w")
    # fileobj.write( str(list_len) + "\n")

if __name__ == "__main__":
    input_file = "allo.txt"
    lig_dict = TreeBuild.parse_lig_file(input_file, "ligandid")

    print lig_dict

    LigandClusteringByClass(lig_dict)

    # distfile = TreeBuild.gen_dist_file(lig_dict, fp.fp_func)
    # newick_o = TreeBuild.run_rapidnj(distfile)
    # dot_inf = TreeBuild.write_dotfile(newick_o)
    # dot_out = TreeBuild.sfdp_dot(dot_inf, 10)
    # dot_dict = TreeBuild.dot2dict(dot_out)
