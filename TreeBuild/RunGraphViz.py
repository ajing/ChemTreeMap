'''
    Run GraphViz on Newick format
'''
from ete2 import Tree
from Model import FILE_FORMAT
import datetime
import subprocess
import os

def WriteDotFile(newick):
    tree = Tree(newick)

    dot_file_name = datetime.datetime.now().strftime(FILE_FORMAT) + ".gv"
    fileobj = open(dot_file_name, "w")

    # rename internal tree name
    i = 0
    for n in tree.traverse():
        if not n.name:
            n.name = "F" + str(i)
            i = i + 1
        else:
            n.name = n.name.replace("\'", "")

    aline = "graph G{\nnode [shape=circle, style=filled];"
    fileobj.write(aline + "\n")
    filecontent = []
    for n in tree.traverse():
        if n.up:
            filecontent.append(n.name + "--" + n.up.name + "[len=" + "{:f}".format(n.dist).rstrip("0") + "]")
        else:
            filecontent.append(n.name)

    fileobj.write("\n".join(filecontent) + "}")
    return dot_file_name


def SFDPonDot(dotfile, size):
    fmt= FILE_FORMAT + '_sfdp.gv'
    newfilename = datetime.datetime.now().strftime(fmt)
    if os.path.isfile(newfilename):
        os.remove(newfilename)
    command = "sfdp -Gsmoothing=triangle -Gsize={size} {infile} > {outfile}".format(size=size, infile=dotfile, outfile=newfilename)
    #print command
    subprocess.Popen( command, shell = True, stdout = subprocess.PIPE ).communicate()
    RemoveBackSlash(newfilename)
    return newfilename

def RemoveBackSlash(dotfile):
    # remove backslash and replace all " quote sign
    f = open(dotfile, 'r+')
    content = f.readlines()
    newcontent = []
    for line in content:
        line = line.replace("\"", "")
        if line.endswith("\\\n"):
            newcontent.append(line[:-2])
        elif line.endswith("\n") and line[-2] != ";":
            newcontent.append(line[:-1])
        else:
            newcontent.append(line)
    f.seek(0)
    f.write("".join(newcontent))
    f.truncate()
    f.close()


## Functions for parsing DOT file
def NodeNameExist(line):
    if "CHEMBL" in line or "ASD" in line or "Chk1" in line or "B" in line or "F" in line:
        return True
    else:
        return False

def IsEdge(line):
    if "--" in line:
        return True
    else:
        return False

def NameAndAttribute(line):
    #print line
    split_index = line.index("[")
    name   = line[:split_index]
    attr   = line[split_index:]
    return name, attr

def SizeScale(size):
    # size is a string
    return float(size)

def GetNodeProperty(line):
    name, attr = NameAndAttribute(line)
    name = ProcessName(name, False)
    position = GetAttributeValue("pos", attr)[:-1].replace(",", "-")
    attr = CleanAttribute(attr)
    width = GetAttributeValue("width", attr)
    #group = GetAttributeValue("color", attr)
    size = SizeScale(GetSize(width))
    return name, size, position


def GetRoot(dotfile, rootname):
    # return root name with most levels
    rootnode  = ""
    for eachline in open(dotfile):
        if NodeNameExist(eachline) and not IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            name = name.strip()
            if name == rootname:
                name, size, position = GetNodeProperty(eachline)
                return Node(name, size = size, position = position)

def ProcessName(name, isedge):
    if isedge:
        firstnode, secondnode = name.split("--")
        firstnode = firstnode.strip()
        secondnode = secondnode.strip()
        return firstnode, secondnode
    else:
        return name.strip()

def GetAttributeValue(attrname, attr):
    left = attr.index("[") + 1
    right = attr.index("]")
    attr  = attr[left:right]
    attrlist = attr.split()
    for each in attrlist:
        if attrname in each:
            value = each.split("=")[1]
            if value.endswith("!"):
                return value[:-1]
            else:
                return value

def CleanAttribute(attr):
    new_attr = attr.replace(",", "")
    return new_attr

def GetSize(width):
    if isinstance(width, str):
        width = float(width)
    return width

def NodeByName(name, contents):
    for eachline in contents:
        if not IsEdge(eachline) and NodeNameExist(eachline):
            nodename, attr = NameAndAttribute(eachline)
            if name == nodename.strip():
                name, size, position = GetNodeProperty(eachline)
                return Node(name, size = size, position = position)

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


def RecursiveNode2Dict(node, info_dict):
    if not node.children:
        result = {"name": node.name, "position": node["position"], "dist": abs(float(node.dist)) / 10}
        result.update(info_dict[node.name])
    else:
        result = {"name": node.name, "position": node["position"], "dist": abs(float(node.dist)) / 10}
        if node.name in info_dict:
            result.update(info_dict[node.name])
    children = [RecursiveNode2Dict(c, info_dict) for c in node.children]
    if children:
        result["children"] = children
    return result

def Dot2Dict(dotfile, moldict):
    rootname = "F0"
    # dotfile is a dot file
    contents = open(dotfile).readlines()
    # get the root of the network
    root = GetRoot(dotfile, rootname)
    curr_nodes = [root]
    curr_name_list = [root.name]
    next_nodes     = 1
    while next_nodes:
        next_nodes = []
        for each_node in curr_nodes:
            next_nodes += ExtendChildren(each_node, contents, curr_name_list)
        curr_nodes = next_nodes
    rootdict = RecursiveNode2Dict(root, moldict)
    return rootdict

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




def test():
    #dot = "./Data/test.gv"
    #SFDPonDot(dot, 10)
    dot = "./Data/sfdp.gv"
    dot = "./test_sfdp.gv"
    RemoveBackSlash(dot)

if __name__ == "__main__":
    sdfile = "./Data/2015-08-24-15h-11m_dot_sfdp.gv"
    filename = WriteDotFile("('B50266893':0.13092,(((((((((((((((('B50266922':0.26611,'B50443857':-0.12731):0.0030516,'B50266920':0.17615):-0.064354,'B50266923':0.13973):0.12317,'B50142090':0.00073348):0.029559,'B50266921':0.13694):0.15116,'B50328697':-0.18416):0.029298,'B50266924':0.14146):0.15875,'B50103666':-0.10715):0.09859,'B50266891':0.017935):0.093078,'B50266743':-0.068619):0.35558,'B50113591':-0.22354):0.2073,'B50142163':-0.18887):0.26968,'B50266919':-0.2017):0.57275,'B50328726':-0.36233):0.55194,'B50266775':-0.53153):1.3155,'B50142161':-1.1255):2.3503,'B50266773':-2.2315);")
    sdfile = SFDPonDot(filename, 10)
    print sdfile
    from Util import ParseLigandFile
    liganddict   = ParseLigandFile("./Data/result_clean_no0_20.txt")  # keep all ligand related information, key is ligand id, value is also a dictionray
    Dot2Dict(sdfile, liganddict)
