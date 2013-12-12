'''
    Dot language to JSON
'''
import sys
sys.path.append("../clusterVis")
from TreeRebuilder import *
from TreeParser import GetLevelFromName
from GetSize import GetSize

def GetRoot(dotfile):
    # return root name with most levels
    maxlevel = 0
    maxnode  = ""
    for eachline in open(dotfile):
        if NodeNameExist(eachline) and not IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            nodelevel = GetLevelFromName(name)
            if nodelevel > maxlevel:
                maxnode  = eachline
                maxlevel = nodelevel
    name, width, group = GetNodeProperty(maxnode)
    size = GetSize(width)
    return Node(name, size = size, group = group)

def GetNodeProperty(line):
    name, attr = NameAndAttribute(line)
    name = ProcessName(name, False)
    attr = CleanAttribute(attr)
    width = GetAttributeValue("width", attr)
    group = GetAttributeValue("color", attr)
    return name, width, group

class Node(dict):
    # class for node of tree, each node can only have one parent
    def __init__(self, name, **attr):
        self.name = name
        self.parent = None
        self.children = []
        self.dist   = 0
        if not attr is None:
            self.update(attr)

    def __str__(self):
        return "a node with name:" + self.name

    def get_dist(self, a_node):
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        if a_node == self.parent:
            return self.dist
        if a_node in self.children:
            return a_node.dist
        else:
            return None

    def add_child(self, a_node):
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        self.children.append(a_node)
        a_node.set_parent(self)

    def set_parent(self, a_node):
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        self.parent = a_node

def Dot2JSON(dotfile):
    # dotfile is a dot file
    contents = open(dotfile).readlines()
    # get the root of the network
    root = GetRoot(dotfile)
    for eachline in contents:
        if NodeNameExist(eachline):
            if not IsEdge(eachline):
                pass

def test():
    testfile = "test.gv"
    print GetRoot(testfile)

if __name__ == "__main__":
    test()
