"""
    For picking ligand names for a branch
"""

from TreeRebuilder import *

def GetLevelFromName(name):
    return len(name.strip().split("_"))

def GetLigandBranch(ligandname, level, treefile):
    for eachline in open(treefile):
        if ligandname in eachline and not IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            if GetLevelFromName(name) > level:
                return name.split("_")

class Node:
    def __init__(self, name, width):
        self.name = name
        self.width = width

def SigNodeParser(treefile, cri_width):
    # input treefile, output node dict class
    node_list = []
    node_size = []
    for eachline in open(treefile):
        if NodeNameExist(eachline) and not IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            width = GetAttributeValue("width", attr)
            if float(width) > cri_width:
                anode = Node(name, width)
                node_list.append(anode)
                node_size.append(int(GetSize(float(width))))
    return node_list, node_size

def SignificantClusters(ligands, nodelist, nodesize):
    sig_ligands = []
    sig_nodesize = []
    for each in ligands:
        for i in range(len(nodesize)):
            if each == nodelist[i].name:
                sig_ligands.append(each)
                sig_nodesize.append(nodesize[i])
    return sig_ligands, sig_nodesize

def GetBranchLargeCluster(ligandname, treefile):
    level = 200
    width_cut = 0.15
    all_ligands_in_cluster = GetLigandBranch(ligandname, level, treefile)
    node_list, node_size = SigNodeParser(treefile, width_cut)
    s_ligands, sig_nodesize = SignificantClusters(all_ligands_in_cluster, node_list, node_size)
    return s_ligands, sig_nodesize

if __name__ == "__main__":
    tree_file = "./Data/all_0.9.gv"
    ligandname = "ASD01910452"
    GetBranchLargeCluster(ligandname, tree_file)
