#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

from distutils.core import setup

setup(
    name='treebuild',
    version='0.1.0',
    packages=['treebuild'],
    url='https://github.com/ajing/ChemTree',
    license='Apache 2.0',
    author='ajing',
    author_email='ajingnk@gmail.com',
    description='Generate Tree Structures for Biochemical Similarity in Molecular Datasets. Please install rdkit and graphviz first, because they are not pip installable.',
    install_requires=[
        #'rdkit',  ## currently rdkit is not pip installable
        #               ## so this dependency must be met by the installer
        'ete2'
    ],
    package_data={'treebuild': ['data/*.txt', 'data/*.py', 'lib/rapidnj-linux-64']}
)
