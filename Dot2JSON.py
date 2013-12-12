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

    def set_dist(self, dist):
        self.dist = dist

def NodeByName(name, contents):
    for eachline in contents:
        if NodeNameExist(eachline) and not IsEdge(eachline) and name in eachline:
            name, width, group = GetNodeProperty(eachline)
            size = GetSize(width)
            return Node(name, size = size, group = group)

def AddNewChild(contents, a_node, new_node_name, edge_length, childrens, currentlist):
    # return a node object
    newnode = NodeByName(new_node_name, contents)
    newnode.set_dist(edge_length)
    childrens.append(newnode)
    currentlist.append(new_node_name)

def ExtendChildren(a_node, contents, cur_list):
    children_list = []
    for eachline in contents:
        if IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            fnode, snode = ProcessName(name, True)
            if fnode == a_node.name and not snode in cur_list:
                edge_len  = GetAttributeValue("len", eachline)
                AddNewChild(contents, a_node, snode, edge_len, children_list, cur_list)
            if snode == a_node.name and not fnode in cur_list:
                edge_len  = GetAttributeValue("len", eachline)
                AddNewChild(contents, a_node, fnode, edge_len, children_list, cur_list)
    return children_list

def Root2JSON(root):
    cur_node = root
    for eachchild in

def Dot2JSON(dotfile):
    # dotfile is a dot file
    contents = open(dotfile).readlines()
    # get the root of the network
    root = GetRoot(dotfile)
    curr_nodes = [root]
    curr_name_list = [root.name]
    next_nodes     = 1
    while next_nodes:
        next_nodes = []
        for each_node in curr_nodes:
            next_nodes += ExtendChildren(each_node, contents, curr_name_list)
        curr_nodes = next_nodes

def test():
    testfile = "test.gv"
    root = GetRoot(testfile)
    print root.__dict__
    print root["size"]
    Dot2JSON(testfile)

if __name__ == "__main__":
    test()
