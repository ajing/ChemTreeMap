#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

"""
    Keep consistent among all files
"""
import os


SMILE_COLUMNNAME = "Canonical_Smiles"

RAPIDNJ_COMMAND = os.path.join(os.path.dirname(__file__), "lib/rapidnj-linux-64")

# temperary files
TMP_FOLDER  = "./.tmp"
FILE_FORMAT = './.tmp/%Y-%m-%d-%Hh-%Mm-%Ss'
## image directory
IMG_DIR = "./images/"

# Potency
POTENCY = "IC50"