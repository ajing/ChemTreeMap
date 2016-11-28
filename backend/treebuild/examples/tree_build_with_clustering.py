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



def LigandClusteringByClass(lig_dict, class_col = "allosteric"):
    """
    Ligand clustering by ligand class

    :param lig_dict: all ligand information
    :param class_col: the column name with ligand class
    :return: new lig_dict with cluster size
    """


def MakeAlloFile(filename):
    """
        
    """

if __name__ == "__main__":
    input_file = "allo.txt"
    lig_dict = TreeBuild.parse_lig_file(input_file, "ligandid")

    # distfile = TreeBuild.gen_dist_file(lig_dict, fp.fp_func)
    # newick_o = TreeBuild.run_rapidnj(distfile)
    # dot_inf = TreeBuild.write_dotfile(newick_o)
    # dot_out = TreeBuild.sfdp_dot(dot_inf, 10)
    # dot_dict = TreeBuild.dot2dict(dot_out)
