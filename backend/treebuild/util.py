'''
    Provide utility functions
'''
import datetime
import json
from csv import DictReader
from rdkit import Chem, DataStructs

from ete2 import Tree

from .model import FILE_FORMAT


def GuessByFirstLine(firstline):
    """
    Guess the number of columns with floats by the first line of the file
    :param firstline:
    :return:
    """
    num_colnam = []
    for key in firstline:
        try:
            float(firstline[key])
            num_colnam.append(key)
        except:
            continue
    return num_colnam

def ConvertToFloat(line, colnam_list):
    """
    Convert some columns (in colnam_list) to float, and round by 3 decimal.
    :param line: a dictionary from DictReader.
    :param colnam_list: float columns
    :return: a new dictionary
    """
    for name in colnam_list:
        line[name] = round(float(line[name]), 3)
    return line

def ParseLigandFile(infile, identifier):
    """
    Parse ligand file to an dictionary, key is ligand id and valud is a dictionary with properties and property values.
    This program will guess the type for each column based on the first row. The program will assume there is only two types of data: number and string.
    :param infile: input filename
    :param identifier: the identifier column name
    :return: a dictionray
    """
    '''

    '''
    mol_dict = dict()
    flag = 1 # first line flag
    id_count = 0
    for line in DictReader(open(infile), delimiter = "\t"):
        if flag:
            num_colnam = GuessByFirstLine({k:v for k,v in line.items() })
        new_id = "B" + str(id_count)
        id_count += 1
        mol_dict[new_id] = ConvertToFloat({k:v for k,v in line.items()}, num_colnam)
        mol_dict[new_id]["orig_id"] = line[identifier]
    return mol_dict


def WriteJSON(dict_obj, outfile, write_type):
    """
    Dump json object to a file
    :param dict_obj: dictionary object
    :param outfile: output file name
    :param write_type: append or rewrite ('a' or 'w')
    :return: void
    """
    fileobj = open(outfile, write_type)
    fileobj.write(json.dumps(dict_obj))


def SelectColumn(lig_dict, colname):
    """
    Prune the dictionary, only attribute in colname will be left.
    :param lig_dict: a tree like dictionary
    :param colname: what attribute you want to keep.
    :return: a new dictionary
    """
    lig_new = dict()
    for k in lig_dict:
        lig_new[k] = {sk:v for sk, v in lig_dict[k].items() if sk in colname}
    return lig_new


def WriteAsPHYLIPFormat(smile_list, fp_func):
    """
    Prepare the input for RapidNJ.
    :param smile_list: a list of smiles string
    :param fp_func: the fingerprint function
    :return: tje filename with PHYLIP format (input for rapidnj)
    """
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
    """
    :param alist: a list of two element list, the first item is ligand name, the second is smile
    :param fp_func: the fingerprint function
    :return: a new list of two element list, with first item as ligand name, second item as a fingerprint object.
    """
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
    """
    Write newick string to a DOT file
    :param newick: a string with newick tree structure
    :return: DOT file name
    """
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
    """
    Rewrite dot file, with removing back slash of dot file
    :param dotfile: DOT file name
    :return: void
    """
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
    """
    Read a DOT file to generate a tree and save it to a dictionary
    :param dotfile: DOT file name
    :param moldict: a dictionary with ligand information
    :return: a dictionary with the tree
    """
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
    """
    Generate similarity score for two smiles strings
    :param fp1: fingerprint object (rdkit)
    :param fp2: fingerprint object (rdkit)
    :return: Tanimoto similarity
    """
    if (fp1 is None or fp2 is None):
        return
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def GetRoot(dotfile, rootname):
    """
    Return root name with rootname
    :param dotfile: DOT file
    :param rootname: the name of the root
    :return: the object of the root
    """
    for eachline in open(dotfile):
        if NodeNameExist(eachline) and not IsEdge(eachline):
            name, attr = NameAndAttribute(eachline)
            name = name.strip()
            if name == rootname:
                name, size, position = GetNodeProperty(eachline)
                return Node(name, size = size, position = position)


def extendChildren(a_node, contents, cur_list):
    """
    Find all children of a node in a tree
    :param a_node: a node in a tree
    :param contents: contents from DOT file
    :param cur_list: current children
    :return: a list of node objects (children)
    """
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
    """
    Whether this line in DOT file is an edge
    :param line: a string line in DOT file
    :return: True or False
    """
    if "--" in line:
        return True
    else:
        return False


def RecursiveNode2Dict(node, info_dict):
    '''
    Recursively populate information to the tree object with info_dict
    :param node: tree object with all info
    :param info_dict: information for each ligand.
    :return: a tree dictionary
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
    """
    Functions for parsing DOT file
    :param line: a line from DOT file
    :return: whether there is a node name in this line
    """
    if "CHEMBL" in line or "ASD" in line or "Chk1" in line or "B" in line or "F" in line:
        return True
    else:
        return False


def NameAndAttribute(line):
    """
    Split name and attribute
    :param line: DOT file name
    :return: name string and attribute string
    """
    split_index = line.index("[")
    name   = line[:split_index]
    attr   = line[split_index:]
    return name, attr


def AddNewChild(contents, a_node, new_node_name, edge_length, children, currentlist):
    """
    Add a new child to a node
    :param contents: a string, a line from DOT
    :param a_node: a node object
    :param new_node_name: new node name
    :param edge_length: the length of edge
    :param children: existing children
    :param currentlist: current list of node name
    :return: void
    """
    # return a node object
    newnode = NodeByName(new_node_name, contents)
    newnode.set_dist(edge_length)
    a_node.add_child(newnode)
    children.append(newnode)
    currentlist.append(new_node_name)


def GetNodeProperty(line):
    """
    Get node property from a string
    :param line: a string
    :return: name, size, and position of the node
    """
    name, attr = NameAndAttribute(line)
    name = ProcessName(name, False)
    position = GetAttributeValue("pos", attr)[:-1].replace(",", "-")
    attr = CleanAttribute(attr)
    width = GetAttributeValue("width", attr)
    #group = GetAttributeValue("color", attr)
    size = SizeScale(GetSize(width))
    return name, size, position


def ProcessName(name, isedge):
    """
    Process the name of the node
    :param name: name of the node
    :param isedge: whether this is a edge
    :return: new name
    """
    if isedge:
        firstnode, secondnode = name.split("--")
        firstnode = firstnode.strip()
        secondnode = secondnode.strip()
        return firstnode, secondnode
    else:
        return name.strip()


def GetAttributeValue(attrname, attr):
    """
    Get node attribute
    :param attrname: name of the attribute
    :param attr: the attribute string
    :return: the value for the attribute
    """
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
    """
    Clean attribute, remove ','.
    :param attr: old attribute string
    :return: new string
    """
    new_attr = attr.replace(",", "")
    return new_attr


def NodeByName(name, contents):
    """
    Create node with name name
    :param name: a string with node name
    :param contents: a list of string from DOT file
    :return: node object
    """
    for eachline in contents:
        if not IsEdge(eachline) and NodeNameExist(eachline):
            nodename, attr = NameAndAttribute(eachline)
            if name == nodename.strip():
                name, size, position = GetNodeProperty(eachline)
                return Node(name, size = size, position = position)


def SizeScale(size):
    """
    Rescale the size (currently only convert to float)
    :param size: a string
    :return: a float
    """
    return float(size)


def GetSize(width):
    """
    Get the size.
    :param width:
    :return:
    """
    if isinstance(width, str):
        width = float(width)
    return width


class Node(dict):
    """
    class for node of tree, each node can only have one parent
    """
    def __init__(self, name, **attr):
        self.name = name
        self.parent = None
        self.children = []
        self.dist   = 0
        self.update(attr)

    def __str__(self):
        return "a node with name:" + self.name

    def get_dist(self, a_node):
        """
        get the node as a dictionary
        :param a_node: Node object
        :return: a dictionary
        """
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        if a_node == self.parent:
            return self.dist
        if a_node in self.children:
            return a_node.dist
        else:
            return None

    def add_child(self, a_node):
        """
        Add child to the node.
        :param a_node: Node object
        :return: void
        """
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        self.children.append(a_node)
        a_node.set_parent(self)

    def set_parent(self, a_node):
        """
        Set the parent for a node
        :param a_node: Node object
        :return: void
        """
        if not isinstance(a_node, Node):
            raise TypeError("argument should be Node class")
        self.parent = a_node

    def set_dist(self, dist):
        """
        set the dictionary attribute for the Node object
        :param dist:
        :return:
        """
        self.dist = dist


if __name__ == "__main__":
    ParseLigandFile("./Data/thrombin_clean_ct.txt")