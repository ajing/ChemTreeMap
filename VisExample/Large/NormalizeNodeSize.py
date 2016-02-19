'''
    Normalize node size, so the size is proportion to the percentage in each data set
'''
import json
from collections import defaultdict

def GetPropRecursive(tree_json, id_dict):
    if not "children" in tree_json:
        id_dict[tree_json["group"]].append(tree_json["size"])
    else:
        for child in tree_json["children"]:
            GetPropRecursive(child, id_dict)

def ChangePropRecursive(tree_json, sum_dict):
    if not "children" in tree_json:
        tree_json["size"] = tree_json["size"] / sum_dict[tree_json["group"]] * 40000
    else:
        for child in tree_json["children"]:
            ChangePropRecursive(child, sum_dict)

def GetSum(json_obj):
    id_dict = defaultdict(list)
    GetPropRecursive(json_obj, id_dict)
    return {key : sum(id_dict[key]) for key in id_dict}


if __name__ == "__main__":
    jsonobj = json.load(open("test.json"))
    sumdict = GetSum(jsonobj)
    ChangePropRecursive(jsonobj, sumdict)
    json.dump(jsonobj, open("test_norm.json", "w"), indent = 2)
