'''
    Dot language to JSON
'''
import sys
sys.path.append("../clusterVis")
from TreeRebuilder import *
from TreeParser import GetLevelFromName

def GetRoot(dotfile):
    maxlevel = 0
    maxnode  = ""
    for eachline in contents:
        if NodeNameExist(eachline):
            name, attr = NameAndAttribute(eachline)
            nodelevel = GetLevelFromName(name)
            if nodelevel > maxlevel:
                maxnode  = name
                maxlevel = nodelevel
    return maxnode


def Dot2JSON(dotfile):
    # dotfile is a dot file
    contents = open(dotfile).readlines()
    # get the root of the network
    for eachline in contents:
        if NodeNameExist(eachline):

