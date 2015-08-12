"""
    Rewrite file to make the node name simple, so graphviz can understand
    Assumption: 1. node names are started with CHEMBL or ASD
                2. only two kinds of statements contain node name: (1) node declaration (2) edge between nodes.
"""

import sys

def NodeNameExist(line):
    if "CHEMBL" in line or "ASD" in line or "Chk1" in line or "B" in line:
        return True
    else:
        return False

def IsEdge(line):
    if "--" in line:
        return True
    else:
        return False

def NameAndAttribute(line):
    print line
    split_index = line.index("[")
    name   = line[:split_index]
    attr   = line[split_index:]
    return name, attr

def ProcessName(name, isedge):
    if isedge:
        firstnode, secondnode = name.split("--")
        firstnode = firstnode.strip()
        secondnode = secondnode.strip()
        return firstnode, secondnode
    else:
        return name.strip()

def static_var(varname,value):
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate

@static_var("counter", 0)
@static_var("hashtable", dict())
def HashANode(nodename):
    if nodename in HashANode.hashtable:
        return "N" + str(HashANode.hashtable[nodename])
    else:
        HashANode.hashtable[nodename] = HashANode.counter
        HashANode.counter += 1
        return "N" + str(HashANode.hashtable[nodename])

def CleanAttribute(attr):
    new_attr = attr.replace(",", "")
    return new_attr

def AddAttributeLabel(attr, label):
    if "label" not in attr:
        return attr
    idx = attr.index("\"")
    return attr[:(idx + 1)] + label + attr[(idx + 1):]

def AddMoreAttribute(attr, labelname, labelvalue):
    right = attr.index("]")
    new_attr = attr[:right] + " " + str(labelname) + "=" + str(labelvalue) + "]\n"
    return new_attr

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

def GetSize(width):
    return 100 ** (width/0.3)

def SimplifyName(aname):
    if "CHEMBL" in aname:
        return aname.replace("CHEMBL", "C")
    if "ASD" in aname:
        return aname.replace("ASD", "A")

@static_var("counter", 0)
@static_var("hashtable", dict())
def HashAName(nodename):
    if nodename in HashAName.hashtable:
        return str(HashAName.hashtable[nodename])
    else:
        HashAName.hashtable[nodename] = HashAName.counter
        HashAName.counter += 1
        return str(HashAName.hashtable[nodename])

def PrintHash(hashtable):
    newdict = dict((y,x) for x,y in hashtable.iteritems())
    for each in newdict:
        print each, ":", newdict[each]

def RewriteDot(infile):
    nodename = dict()
    newfilename = infile + "_simple"
    newfileobj  = open(newfilename, "w")
    for eachline in open(infile):
        if NodeNameExist(eachline):
            name, attr = NameAndAttribute(eachline)
            attr = CleanAttribute(attr)
            if IsEdge(eachline):
                fnode, snode = ProcessName(name, True)
                fnode_new = HashANode(fnode)
                snode_new = HashANode(snode)
                new_line  = "--".join([fnode_new, snode_new]) + attr
            else:
                node = ProcessName(name, False)
                node_new = HashANode(node)
                #if not "_" in node and float(GetAttributeValue("width", attr)) > 0.15:
                ##for display large node
                #    new_name = HashAName(node)
                #    attr = AddAttributeLabel(attr, new_name)
                #    attr = AddMoreAttribute(attr, "fixedsize", "true")
                new_line = node_new + attr
            newfileobj.write(new_line)
        else:
            newfileobj.write(eachline)
    newfileobj.close()
    PrintHash(HashAName.hashtable)

def GetMaxWidth(infile):
    widthlist = []
    for eachline in open(infile):
        if NodeNameExist(eachline):
            name, attr = NameAndAttribute(eachline)
            if not IsEdge(eachline):
                widthlist.append(float(GetAttributeValue("width", attr)))
    print max(widthlist)

if __name__ == "__main__":
    infile = sys.argv[1]
    print infile
    RewriteDot(infile)
    #GetMaxWidth(infile)
