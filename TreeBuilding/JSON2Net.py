'''
   JSON to Edge and Node files
   The assumption is that each json node has "name", "dist" attribute. Parent node has "children" attribute.
'''

import json
import csv

def RecursiveJSON2Net(node, node_info, edge_info):
    if not "children" in node:
        one_node = dict()
        for one_key in node:
            if one_key is not "children":
                one_node[one_key] = node[one_key]
        node_info.append(one_node)
    else:
        for c_node in node["children"]:
            a_edge = {"orgin": node["name"], "dest": c_node["name"], "dist": c_node["dist"]}
            edge_info.append(a_edge)
            RecursiveJSON2Net(c_node, node_info, edge_info)

def WriteDictToFile(filename, dict_list):
    with open(filename, "wb") as fileobj:
        w = csv.DictWriter(fileobj, dict_list[0].keys())
        w.writeheader()
        w.writerows(dict_list)

if __name__ == "__main__":
    jsonfile = "test.json"
    jsonobj  = json.load(open(jsonfile))
    node_info = []
    edge_info = []
    RecursiveJSON2Net(jsonobj, node_info, edge_info)
    WriteDictToFile("nodeinfo.csv", node_info)
    WriteDictToFile("edgeinfo.csv", edge_info)
