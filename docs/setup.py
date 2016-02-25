#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

from distutils.core import setup, Extension

from setuptools.command.install import install
from distutils import log # needed for outputting information messages

import os
import sys
os.chdir('../backend')


class OverrideInstall(install):

    def run(self):
        uid, gid = 0, 0
        mode = 0700
        install.run(self) # calling install.run(self) insures that everything that happened previously still happens, so the installation does not break!
        # here we start with doing our overriding and private magic ..
        print self.install_scripts
        for filepath in self.get_outputs():
            if filepath.endswith("rapidnj-linux-64"):
                log.info("Overriding setuptools mode of scripts ...")
                log.info("Changing ownership of %s to uid:%s gid %s" %
                         (filepath, uid, gid))
                os.chown(filepath, uid, gid)
                log.info("Changing permissions of %s to %s" %
                         (filepath, oct(mode)))
                os.chmod(filepath, mode)

setup(
    name='treebuild',
    version='0.1.0',
    packages=['treebuild'],
    url='https://github.com/ajing/ChemTreeMap',
    license='Apache 2.0',
    author='ajing',
    author_email='ajingnk@gmail.com',
    description='Generate Tree Structures for Biochemical Similarity in Molecular Datasets. Please install rdkit and graphviz first, because they are not pip installable.',
    install_requires=[
        #'rdkit',  ## currently rdkit is not pip installable
        #               ## so this dependency must be met by the installer
        'ete2'
    ],
    package_data={'treebuild': ['data/*.txt', 'data/*.py', 'lib/rapidnj-linux-64']},
    cmdclass={'install': OverrideInstall}
)
