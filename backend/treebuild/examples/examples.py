#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

import treebuild
from treebuild import DEFAULT_ACTIVITY_TYPES, DEFAULT_PROPERTY_TYPES, DEFAULT_FINGERPRINT_TYPES, DEFAULT_EXTERNAL
from treebuild.types import bindingdb, pubchem

# affinity data with chk1 kinase


# factor xa
input_file = "./factorxa.txt"
out_file = "./factorxa.json"
properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": DEFAULT_EXTERNAL}
treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)

# cdk2
input_file = "./cdk2.txt"
out_file = "./cdk2.json"
properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": DEFAULT_EXTERNAL}
treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)

# MAP P38
input_file = "./map_p38.txt"
out_file = "./map_p38.json"
properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": [bindingdb]}
treebuild.TreeBuild(input_file, out_file, id_column ="BindingDB", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)

# cytochrome P450
input_file = "./cyto.txt"
out_file = "./cyto.json"
properties = {"activities": DEFAULT_ACTIVITY_TYPES, "properties": DEFAULT_PROPERTY_TYPES, "ext_links": [pubchem]}
treebuild.TreeBuild(input_file, out_file, id_column ="PubChem", fps = DEFAULT_FINGERPRINT_TYPES, properties=properties)

