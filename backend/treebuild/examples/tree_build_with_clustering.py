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
        print "cluster center:", dir(kcluster.cluster_centers_)
        print "cluster labels:", kcluster.labels_

        for c_idx in range(kcluster.cluster_centers_.shape[0]):
            min_dist = float("inf")
            dist = (fp_mat - kcluster.cluster_centers_[c_idx,])**2
            dist = np.sum(dist, axis=1)
            print "dist:", dist
            m_idx = np.argmin(dist)
            print "m_idx", m_idx
            lig_dict_center[ids[m_idx]] = lig_dict[ids[m_idx]]
            lig_dict_center[ids[m_idx]]["cluster_size"] = sum(kcluster.labels_ == c_idx)

    return lig_dict_center










    # to fp matrix
    # KMean clustering
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

    LigandClusteringByClass(lig_dict)

    # distfile = TreeBuild.gen_dist_file(lig_dict, fp.fp_func)
    # newick_o = TreeBuild.run_rapidnj(distfile)
    # dot_inf = TreeBuild.write_dotfile(newick_o)
    # dot_out = TreeBuild.sfdp_dot(dot_inf, 10)
    # dot_dict = TreeBuild.dot2dict(dot_out)
