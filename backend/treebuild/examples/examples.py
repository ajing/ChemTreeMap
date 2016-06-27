#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods
import os

import treebuild
import shutil
from distutils.dir_util import copy_tree
from treebuild import DEFAULT_ACTIVITY_TYPES, DEFAULT_PROPERTY_TYPES, DEFAULT_FINGERPRINT_TYPES, DEFAULT_EXTERNAL
from treebuild.types import bindingdb, pubchem, pic50

# for the static server
import SimpleHTTPServer
import SocketServer

# directory for chemtree
CHEMTREE_DIR = "../../../"


input_file = "./.txt"
out_file = "./aff.json"
properties = {"activities": [pic50], "properties": DEFAULT_PROPERTY_TYPES, "ext_links": []}
treebuild.TreeBuild(input_file, out_file, id_column="ligandid", fps=DEFAULT_FINGERPRINT_TYPES, properties=properties)

# affinity data with chk1 kinase
#input_file = "./aff.txt"
#out_file = "./aff.json"
#properties = {"activities": [pic50], "properties": DEFAULT_PROPERTY_TYPES, "ext_links": []}
#treebuild.TreeBuild(input_file, out_file, id_column="ligandid", fps=DEFAULT_FINGERPRINT_TYPES, properties=properties)

# # factor xa
#input_file = "./factorxa.txt"
#out_file = "./factorxa.json"
#properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": DEFAULT_EXTERNAL}
#treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)
#
# # cdk2
# input_file = "./cdk2.txt"
# out_file = "./cdk2.json"
# properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": DEFAULT_EXTERNAL}
# treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)
#
# # MAP P38
# input_file = "./map_p38.txt"
# out_file = "./map_p38.json"
# properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": [bindingdb]}
# treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)
#
# # cytochrome P450
# input_file = "./cyto.txt"
# out_file = "./cyto.json"
# properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": [pubchem]}
# treebuild.TreeBuild(input_file, out_file, id_column ="PubChem", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)
#
# # move the image directory to frontend
for fname in ["./aff.json", "./factorxa.json", "./cdk2.json", "./map_p38.json", "./cyto.json"]:
    try:
        shutil.copy(fname, os.path.join(CHEMTREE_DIR, "frontend/app/data/"))
    except:
        continue
if os.path.exists(os.path.join(CHEMTREE_DIR, "frontend/app/images")):
    shutil.rmtree(os.path.join(CHEMTREE_DIR, "frontend/app/images"))
copy_tree("./images", os.path.join(CHEMTREE_DIR, "frontend/app/images"))

##
# setup the server
##
# change working directory
os.chdir(os.path.join(CHEMTREE_DIR, "frontend"))

PORT = 8000

Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

httpd = SocketServer.TCPServer(("", PORT), Handler)

print "Serving at port", PORT
print "Please open the following link in your browser:"
print "localhost:8000/dist/#/" + out_file.split("/")[-1].split(".")[0]

httpd.serve_forever()
