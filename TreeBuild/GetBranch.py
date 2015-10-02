'''
    Get the molecules within one branch
'''

import json
import sys

def FindClosestRoot(tree_json, root_id):
    if not u"children" in tree_json:
        return None
    else:
        for child in tree_json["children"]:
            if child["name"] == root_id:
                return child
            sub_tree = FindClosestRoot(child, root_id)
            if sub_tree:
                return sub_tree

def GetIDRecursive(tree_json, id_list):
    if not "children" in tree_json:
        id_list.append(tree_json["name"])
    else:
        for child in tree_json["children"]:
            GetIDRecursive(child, id_list)

def GetIDs(tree_json, compound_list):
    #  the seed id is the node which is closest to the root of the branch
    id_list = []
    GetIDRecursive(tree_json, id_list)
    org_list = []
    for a_id in id_list:
        for compound in compound_list:
            if compound["id"] == a_id:
                org_list.append(compound["orig_id"])
                break
    return org_list

def GenFile(filename, id_list, num):
    namecontent = filename.split(".")
    newfilename = namecontent[0] + "_small" + num + "." + namecontent[1]
    newobj  = open(newfilename, "w")
    newcontent = []
    fline = 1
    new_id_list = []
    for x in id_list:
        if x[0] == "B":
            new_id_list.append(x[1:])
        else:
            new_id_list.append(x)
    for line in open(filename, "r"):
        content = line.strip().split("\t")
        if fline:
            fline = 0
            content.append("operator")
            if "Canonical_Smiles" in content:
                content[content.index("Canonical_Smiles")] = "Smiles"
            newcontent.append("\t".join(content))
        if content[0] in new_id_list:
            content.append("=")
            newcontent.append("\t".join(content))
    newobj.write("\n".join(newcontent))

if __name__ == "__main__":
    jsonfile = sys.argv[1]
    root_name= sys.argv[2]
    filename = sys.argv[3]
    fnumber  = sys.argv[4]
    jsonobj = json.load(open(jsonfile))
    treeobj = jsonobj["trees"]["ECFP"]
    root_json = FindClosestRoot(treeobj, root_name)
    print "root json is found", root_json["name"]
    org_list   = GetIDs(root_json, jsonobj["compounds"])
    print "org list is", org_list
    GenFile(filename, org_list, fnumber)
