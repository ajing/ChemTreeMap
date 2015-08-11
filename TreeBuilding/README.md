#clusterVis

clusterVis is a set of tools for tree building from molecules and preparing input file for TreeViz.

## Function for each file

***dependency.sh*** is a simple way to install all dependencies for tree building. This is for Ubuntu operating system. This file may not be applicable to other operating system.

***CreateGraphNumpy.py*** are main function for creating a tree from a molecular library. The molecular library is a tab delimited file with canonical smile string. You can check the file under TreeViz/ligand_5_7.txt as an example.

***ligandGraphall.py*** and ***ligandGraphnames*** are for creating distance matrix from a list of SMILES strings. You can change fingerprints and similarity functions inside.

***BuildTree.py*** is for building tree from distance matrix. The output file will saved in __Data__ directory.

***MakeStructuresForSmiles.py*** can generate images for molecules (smile string).

***TreeParser*** is for helping analyze the tree structure. You can use the function to find which branch a molecule belongs to.
