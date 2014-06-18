'''
  Internal node to end node mapping
'''

import json
import shutil
import os

def FindEndNode(root):
  if u'children' in root:
    next_node = root[u'children'][0]
    return FindEndNode(next_node)
  else:
    return root

def MapToAllEndNode(root):
  acc = dict()
  def MapAux(root, accu):
    if u'children' in root:
      acc[root[u'name']] = FindEndNode(root)[u'name']
      for eachchild in root[u'children']:
        MapAux(eachchild, acc)
  MapAux(root, acc)
  return acc

def CopyNewFiles(map_dict, filedir):
  os.chdir(filedir)
  for eachkey in map_dict:
    shutil.copyfile(map_dict[eachkey], eachkey)

def Dict2JSON(dictionary, filename):
  fileobj  = open(filename, "w")
  fileobj.write(json.dumps(dictionary, indent=2))

if __name__ == "__main__":
  root_dict = json.load(open("./Large/all_0.9.json"))
  acc = MapToAllEndNode(root_dict)
  #Dict2JSON(acc, "mapping.json")
  CopyNewFiles(acc, "./Image")
