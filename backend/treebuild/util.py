'''
    Provide utility functions
'''
import datetime
import json
import math
from csv import DictReader
from rdkit import Chem, DataStructs

from ete2 import Tree

from treebuild.model import FILE_FORMAT


def GuestByFirstLine(firstline):
    num_colnam = []
    for key in firstline:
        try:
            float(firstline[key])
            num_colnam.append(key)
        except:
            continue
    return num_colnam

def ConvertToFloat(line, colnam_list):
    for name in colnam_list:
        line[name] = round(float(line[name]), 3)
    return line

def ParseLigandFile(infile, identifier):
    '''
       parse ligand file to an dictionary, key is ligand id and valud is a dictionray with properties and property values
       This program will guess the type for each column based on the first row. The program will assume there is only two type of data: number and string.
    '''
    mol_dict = dict()
    flag = 1 # first line flag
    id_count = 0
    for line in DictReader(open(infile), delimiter = "\t"):
        if flag:
            num_colnam = GuestByFirstLine({k:v for k,v in line.items() })
        new_id = "B" + str(id_count)
        id_count += 1
        mol_dict[new_id] = ConvertToFloat({k:v for k,v in line.items()}, num_colnam)
        mol_dict[new_id]["orig_id"] = line[identifier]
    return mol_dict


def WriteJSON(dict_obj, outfile, write_type):
    fileobj = open(outfile, write_type)
    fileobj.write(json.dumps(dict_obj))

def SelectColumn(lig_dict, colname):
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
    return lig_new

def ReArrangeActivity(lig_dict, colname):
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if not sk in colname}
        lig_new[k]["activities"] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
        if "IC50" in lig_new[k]["activities"]:
            lig_new[k]["activities"]["pIC50"] = round(9 - math.log10(float(lig_new[k]["activities"]["IC50"])), 5)
            #lig_new[k]["activities"]["pIC50"] = round(6 - math.log10(float(lig_new[k]["activities"]["IC50"])), 5)
    return lig_new

def Dict2List(lig_dict):
    lig_list = []
    for idx in range(len(lig_dict.keys())):
        new_entry = lig_dict["B" + str(idx)]
        new_entry["id"] = "B" + str(idx)
        lig_list.append(new_entry)
    return lig_list


def WriteAsPHYLIPFormat(smile_list, fp_func):
    fp_list = ToFPObj(smile_list, fp_func)
    print "finish parsing smile list"
    list_len = len(fp_list)

    newfilename = datetime.datetime.now().strftime(FILE_FORMAT) + ".dist"
    fileobj  = open(newfilename, "w")
    fileobj.write( str(list_len) + "\n")

    for i in range(list_len):
        lig1 = fp_list[i]
        lig1list = []
        for j in range(list_len):
            lig2 = fp_list[j]
            sim  = getSimilarity(lig1[1], lig2[1])
            lig1list.append([lig2[0], 1 - sim])

        sim_values = [ "%.4f" % x[1] for x in lig1list]
        line = "\t".join([lig1[0], "\t".join(sim_values)]) + "\n"
        fileobj.write(line)

    fileobj.close()

    return newfilename


def ToFPObj(alist, fp_func):
    # for alist, the first item is ligand name, the second is smile
    newlist = []
    for each in alist:
        smile = each[1]
        m = Chem.MolFromSmiles(smile)
        if m is None:
            continue
        fp = fp_func(m)
        newlist.append([each[0], fp])
    return newlist


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
            next_nodes += extendChildren(each_node, contents, curr_name_list)
        curr_nodes = next_nodes
    rootdict = RecursiveNode2Dict(root, moldict)
    return rootdict


def getSimilarity(fp1, fp2):
    # generate similarity score for two smiles strings
    if (fp1 is None or fp2 is None):
        return
    return DataStructs.TanimotoSimilarity(fp1, fp2)


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


def extendChildren(a_node, contents, cur_list):
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


def IsEdge(line):
    if "--" in line:
        return True
    else:
        return False


def RecursiveNode2Dict(node, info_dict):
    '''
    :param node: tree object with all info
    :param info_dict:
    :return:
    '''
    if not node.children:
        x, y   = map(float, node["position"].split("-"))
        result = {"name": node.name, "size": 1, "x": x, "y": y, "dist": abs(float(node.dist))}
        if info_dict:
            result.update(info_dict[node.name])
    else:
        x, y   = map(float, node["position"].split("-"))
        result = {"name": node.name, "x": x, "y": y, "dist": abs(float(node.dist))}
        if info_dict and node.name in info_dict:
            result.update(info_dict[node.name])
    children = [RecursiveNode2Dict(c, info_dict) for c in node.children]
    if children:
        result["children"] = children
    return result


def NodeNameExist(line):
    ## Functions for parsing DOT file
    if "CHEMBL" in line or "ASD" in line or "Chk1" in line or "B" in line or "F" in line:
        return True
    else:
        return False


def NameAndAttribute(line):
    #print line
    split_index = line.index("[")
    name   = line[:split_index]
    attr   = line[split_index:]
    return name, attr


def AddNewChild(contents, a_node, new_node_name, edge_length, childrens, currentlist):
    # return a node object
    newnode = NodeByName(new_node_name, contents)
    newnode.set_dist(edge_length)
    a_node.add_child(newnode)
    childrens.append(newnode)
    currentlist.append(new_node_name)


def GetNodeProperty(line):
    name, attr = NameAndAttribute(line)
    name = ProcessName(name, False)
    position = GetAttributeValue("pos", attr)[:-1].replace(",", "-")
    attr = CleanAttribute(attr)
    width = GetAttributeValue("width", attr)
    #group = GetAttributeValue("color", attr)
    size = SizeScale(GetSize(width))
    return name, size, position


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


def NodeByName(name, contents):
    for eachline in contents:
        if not IsEdge(eachline) and NodeNameExist(eachline):
            nodename, attr = NameAndAttribute(eachline)
            if name == nodename.strip():
                name, size, position = GetNodeProperty(eachline)
                return Node(name, size = size, position = position)


def SizeScale(size):
    # size is a string
    return float(size)


def GetSize(width):
    if isinstance(width, str):
        width = float(width)
    return width


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


if __name__ == "__main__":
    ParseLigandFile("./Data/thrombin_clean_ct.txt")