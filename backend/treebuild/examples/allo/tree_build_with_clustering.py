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

import numpy as np
import os
import json

def SMILE2Matrix(smile_list):
    # To ECFP6
    def ToECFP(id_smile):
        cid = id_smile[0]
        smile = id_smile[1]
        mol = Chem.MolFromSmiles(smile)
        return [cid, AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)]

    fps = map(ToECFP, smile_list)

    np_fps = []
    ids = []
    for fp in fps:
      arr = np.zeros((1,))
      vfp = DataStructs.FoldFingerprint(fp[1], 4)
      DataStructs.ConvertToNumpyArray(vfp, arr)
      ids.append(fp[0])
      np_fps.append(arr)

    return  ids, np_fps



def LigandClusteringByClass(lig_dict, class_col = "allosteric", smile_col = "Canonical_Smiles", num_clusters = {"allosteric": 5, "competitive" : 3}):
    """
    Ligand clustering by ligand class

    :param lig_dict: all ligand information
    :param class_col: the column name with ligand class
    :param smile_col: the column name with SMILE string
    :param num_clusters: the number of clusters for each class
    :return: new lig_dict with center only and cluster size
    """
    all_classes = set([lig_dict[lig_name][class_col] for lig_name in lig_dict.keys()])
    lig_dict_center = dict()

    for e_class in num_clusters.keys():
        smile_list = [ [lig_name, lig_dict[lig_name][smile_col]] for lig_name in lig_dict.keys() if lig_dict[lig_name][class_col] == e_class]
        ids, fp_mat = SMILE2Matrix(smile_list)
        kcluster = KMeans(n_clusters = num_clusters[e_class]).fit(fp_mat)
        print "cluster labels:", kcluster.labels_
        print "cluster ids:", ids

        for c_idx in range(kcluster.cluster_centers_.shape[0]):
            dist = (fp_mat - kcluster.cluster_centers_[c_idx,])**2
            dist = np.sum(dist, axis=1)
            print "dist:", dist
            m_idx = np.argmin(dist)
            print "m_idx", m_idx
            lig_dict_center[ids[m_idx]] = lig_dict[ids[m_idx]]
            lig_dict_center[ids[m_idx]]["cluster_size"] = sum(kcluster.labels_ == c_idx)

    return lig_dict_center

def RecursiveUpdate(node, info_dict):
    '''
    Recursively populate information to the tree object with info_dict.
    :param node: tree object with all info
    :param info_dict: information for each ligand.
    :return: a tree dictionary
    '''
    children = None
    if not "children" in node:
        name = node["name"]
        node.update(info_dict[name])
    else:
        children = [RecursiveUpdate(c, info_dict) for c in node["children"]]


if __name__ == "__main__":
    TMP_FOLDER  = "./.tmp"
    IMG_DIR = "./images/"

    input_file = "allo.txt"
    output_file = "allo.json"
    lig_dict = TreeBuild.parse_lig_file(input_file, "ligandid")

    lig_dict_center = LigandClusteringByClass(lig_dict, num_clusters = {"allosteric": 5, "competitive" : 3})

    if not os.path.exists(TMP_FOLDER):
        os.makedirs(TMP_FOLDER)
    if not os.path.exists(IMG_DIR):
        os.makedirs(IMG_DIR)

    distfile = TreeBuild.gen_dist_file(lig_dict, lambda mol: AllChem.GetMorganFingerprint(mol, 3))
    newick_o = TreeBuild.run_rapidnj(distfile)
    dot_inf = TreeBuild.write_dotfile(newick_o)
    dot_out = TreeBuild.sfdp_dot(dot_inf, 10)
    dot_dict = TreeBuild.dot2dict(dot_out)

    print lig_dict
    print dot_dict
    RecursiveUpdate(dot_dict, lig_dict)
    print "after:", dot_dict

    fileobj = open(output_file, "w")
    fileobj.write(json.dumps(dot_dict))

    # make image file
    TreeBuild.make_structures_for_smiles(lig_dict)
