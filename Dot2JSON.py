'''
    Dot language to JSON
'''
import sys
sys.path.append("../clusterVis")
from TreeRebuilder import *
from TreeParser import GetLevelFromName
from GetSize import GetSize
import json

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
    name, size, group, position = GetNodeProperty(maxnode)
    return Node(name, size = size, group = group, position = position)

def SizeScale(size):
    # size is a string
    return 1000*float(size)

def GetNodeProperty(line):
    name, attr = NameAndAttribute(line)
    name = ProcessName(name, False)
    position = GetAttributeValue("pos", attr)[:-1].replace(",", "-")
    print position
    attr = CleanAttribute(attr)
    width = GetAttributeValue("width", attr)
    group = GetAttributeValue("color", attr)
    size = SizeScale(GetSize(width))
    return name, size, group, position

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
        if not IsEdge(eachline) and name in eachline:
            name, size, group, position = GetNodeProperty(eachline)
            return Node(name, size = size, group = group, position = position)

def AddNewChild(contents, a_node, new_node_name, edge_length, childrens, currentlist):
    # return a node object
    newnode = NodeByName(new_node_name, contents)
    newnode.set_dist(edge_length)
    a_node.add_child(newnode)
    childrens.append(newnode)
    currentlist.append(new_node_name)

def ExtendChildren(a_node, contents, cur_list):
    children_list = []
    for eachline in contents:
        if IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            fnode, snode = ProcessName(name, True)
            eachline = CleanAttribute(eachline)
            if fnode == a_node.name and not snode in cur_list:
                edge_len  = GetAttributeValue("len", eachline)
                AddNewChild(contents, a_node, snode, edge_len, children_list, cur_list)
            if snode == a_node.name and not fnode in cur_list:
                edge_len  = GetAttributeValue("len", eachline)
                AddNewChild(contents, a_node, fnode, edge_len, children_list, cur_list)
    return children_list

def static_var(varname, value):
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate

@static_var("namedict", dict())
@static_var("counter", 0)
def SimpleName(name):
    SimpleName.counter += 1
    if not "_" in name:
        return name
    if not name in SimpleName.namedict:
        SimpleName.namedict[name] = SimpleName.counter
        return "N" + str(SimpleName.counter)
    return SimpleName.namedict[name]

def RecursiveNode2Dict(node):
    if not node.children:
        result = {"name": SimpleName(node.name), "size": node["size"], "group": node["group"], "position": node["position"], "dist": float(node.dist)}
    else:
        result = {"name": SimpleName(node.name), "position": node["position"], "dist": float(node.dist)}
    children = [RecursiveNode2Dict(c) for c in node.children]
    if children:
        result["children"] = children
    return result

def Root2JSON(root, filename):
    fileobj  = open(filename, "w")
    rootdict = RecursiveNode2Dict(root)
    fileobj.write(json.dumps(rootdict, indent=2))

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
    return root

def test():
    testfile = "./Data/test.gv"
    #root = GetRoot(testfile)
    #print root.__dict__
    #print root["size"]
    root = Dot2JSON(testfile)
    Root2JSON(root, "./Data/test.json")
    #aline = "ASD00030001 [label=, width=0.027778, color=red, pos=5.6889,266.3, height=0.041667];\n"
    #print GetNodeProperty(aline)

if __name__ == "__main__":
    test()
    #infile = sys.argv[1]
    #root = Dot2JSON(infile)
    #Root2JSON(root, "test.json")
