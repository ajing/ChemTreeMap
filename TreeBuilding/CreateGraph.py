"""
    Create graph from dmatrix
!!!!!!!!!!!!!!!! Please use CreateGraphNumpy.py, this file is obsolete !!!!!!!!!!!!!!!
"""
import networkx as nx
import numpy as np
import os
import pickle
from random import choice
from FileParser import GetAllinfo
from BuildTree import BuildTree
# Global variable so I can change easily
__SAVEDIR__ = "./Data/"

def createGraph( dmatrix, criteria, moldict, typeofbinding ):
    newGraph = nx.Graph()
    row, col = dmatrix.shape
    for eachrow in range(row):
        # in case, isolated nodes exist in graph
        if not IsTypeofBinding(eachrow, moldict, typeofbinding):
            continue
        newGraph.add_node( eachrow )
        for eachcol in range( eachrow + 1, col ):
            if not IsTypeofBinding(eachcol, moldict, typeofbinding):
                continue
            eachweight = dmatrix[eachrow, eachcol]
            if eachweight < criteria:
                newGraph.add_edge( eachrow, eachcol, weight = eachweight )
    return newGraph

def RandomPickFromList( alist ):
    return choice( alist )

def BindingTypeFilter( alist, moldict, bindingType = None ):
    if bindingType is None:
        return alist
    newlist = []
    for eachIndex in alist:
        if IsTypeofBinding( eachIndex, moldict, bindingType):
            newlist.append( eachIndex )
    if len(newlist) == 0:
        raise RuntimeError("Empty list after filtering!!")
    return newlist

def LeaderInCluster( graphObj, moldict ):
    leaderList = []
    for eachSubGraph in nx.connected_component_subgraphs( graphObj ):
        centerlist = nx.center(eachSubGraph)
        if len(centerlist) < 1:
            raise "A subgraph with less than 1 center"
        leaderID   = RandomPickFromList(centerlist)
        leaderList.append( leaderID )
        graphSize  = len(eachSubGraph)
        moldict[ leaderID ][ "size" ] = graphSize
    return leaderList

def IsTypeofBinding( index, moldict, typeofbinding ):
    if typeofbinding == "all":
        return True
    return moldict[index]["typeofbinding"] == typeofbinding

def MoleculeDictionary( infile ):
    all_info = GetAllinfo( infile )
    totalrow = len( all_info[ "ligandid" ] )
    molDict  = dict()
    for index in range( totalrow ):
        molDict[index] = dict()
        for dict_name in all_info.keys():
            molDict[index][dict_name] = all_info[dict_name][index]
    return molDict

def CheckExistingLeaderlist( typeofbinding, criteria ):
    criString = str(criteria)
    expectFile = typeofbinding + "_" + criString
    numpyarray = None
    moldictfile = None
    for eachfile in os.listdir(__SAVEDIR__):
        if typeofbinding in eachfile and criString in eachfile:
            filedir = __SAVEDIR__ + eachfile
            if eachfile == expectFile + ".npz":
                numpyarray = filedir
            if eachfile == expectFile + ".p":
                moldictfile = filedir
    if not numpyarray is None and not moldictfile is None:
        return (numpyarray, moldictfile)
    return None

def SizeHistogram( moldict ):
    sizeList = []
    for eachID in moldict:
        try:
            sizeList.append(moldict[eachID]["size"])
        except:
            continue
    value, binedge = np.histogram(sizeList)
    print "number:", value
    print "value:", binedge

def SaveLeaderAndMolDict(leaderlist, moldict, typeofbinding, criteria):
    criString = str(criteria)
    filename = "_".join([ typeofbinding, criString ])
    filedir  = __SAVEDIR__ + filename
    np.savez( filedir, leaderlist )
    pickle.dump( moldict, open( filedir + ".p", "wb" ))

def main( bindingtype, minDistance, dmatrix ):
    #minDistance = 0.75
    infile      = "./Data/ligand_5_7_ppilot_modified.txt"
    #leaderAndmol = CheckExistingLeaderlist( bindingtype, minDistance )
    leaderAndmol = None
    if leaderAndmol:
        leaderfile, moldictfile = leaderAndmol
        with np.load(leaderfile) as leader_moldict:
            print leader_moldict.files
            leaderlist     = leader_moldict["arr_0"]
            print "leaderlist length:", len(leaderlist)
        with open(moldictfile, "rb") as moldictobj:
            moldict        = pickle.load(moldictobj)
        SizeHistogram( moldict )
    else:
        moldict = MoleculeDictionary( infile )
        SizeHistogram( moldict )
        newgraph = createGraph( dmatrix, minDistance, moldict, bindingtype)
        leaderlist = LeaderInCluster( newgraph, moldict )
        SaveLeaderAndMolDict(leaderlist, moldict, bindingtype, minDistance)
    leaderlist = BindingTypeFilter( leaderlist, moldict, bindingtype)
    BuildTree( leaderlist, dmatrix, moldict, str(minDistance) + "_" + bindingtype )

if __name__ == "__main__":
    smatrixfile = "./Data/similarityMatrix.npy"
    dmatrix = 1 - np.load(smatrixfile)
    print dmatrix.shape
    #distanceList = [ 0.6, 0.65, 0.7, 0.8 ]
    distanceList = [ 0.1 ]
    #for each in ["allosteric", "competitive"]:
    for each in ["all"]:
        for distance in distanceList:
            main(each, distance, dmatrix)
