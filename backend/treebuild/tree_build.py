#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

import datetime
import os
import subprocess
import shutil
from rdkit import Chem
from rdkit.Chem.Draw import MolToFile

from .types import FingerPrintType
from .util import ParseLigandFile, WriteJSON, WriteAsPHYLIPFormat, Dot2Dict, \
    WriteDotFile, RemoveBackSlash
from .model import IMG_DIR, SMILE_COLUMNNAME, RAPIDNJ_COMMAND, FILE_FORMAT, TMP_FOLDER

class TreeBuild:
    """
    There are assumptions for the data format of the input file.
    It is very important to understand these assumptions:
        1. potency (e.g. IC50/Ka/Ki) unit is nM
        2. the file must have a id column, you can set the column name with id_column
        3. the file must have a SMILES column, with 'Canonical_Smiles' as column name
        4. the file must have at least one potency column (IC50/Ka/Ki).
    To build the tree
        1. the identity column needs to be specified with id_column
        2. a list of fingerprints and a list of properties need to be specified with rdkit
        3. the directories for input and output file need to be specified
    """
    def __init__(self, input_file, output_file, id_column, fps, properties):
        """Setting parameters to build the tree.

        :param input_file: input file is a tab delimited text file.
        :param output_file: output file is a json file
        :param id_column: the id for each column, which will shown as the identifier in the visualization.
        :param fps: a list of FingperPrintType
        :param properties: a list of PropertyType
        :return: void, the program will generate input file for the visualization.
        """
        # initial setting
        self._RAPIDNJ_COMMAND = RAPIDNJ_COMMAND
        self._FILE_FORMAT = FILE_FORMAT

        # creating folders
        if not os.path.exists(TMP_FOLDER):
            os.makedirs(TMP_FOLDER)
        if not os.path.exists(IMG_DIR):
            os.makedirs(IMG_DIR)

        activities = properties["activities"]
        other_properties = properties["properties"]
        ext_links = properties["ext_links"]
        lig_dict = self.parse_lig_file(input_file, id_column)
        trees = dict()
        for fp in fps:
            assert isinstance(fp, FingerPrintType)
            trees[fp.name] = self._build_single_tree(lig_dict, fp)
        metadata = dict()
        metadata["activityTypes"] = [act.to_dict() for act in activities]
        metadata["treeTypes"] = [fp.to_dict() for fp in fps]
        metadata["circleSizeTypes"] = [prop.to_dict() for prop in other_properties]
        metadata["circleBorderTypes"] = [prop.to_dict() for prop in other_properties]
        metadata["external"] = ext_links

        ext_names = [ext["name"] for ext in ext_links]

        comp_info = self.gen_properties(lig_dict, activities, other_properties, ext_names)
        final_dict = {"metadata": metadata, "trees": trees, "compounds": comp_info}

        WriteJSON(final_dict, outfile=output_file, write_type="w")
         # make image file
        self.make_structures_for_smiles(lig_dict)

        # delete tmp folder
        shutil.rmtree(TMP_FOLDER)


    def _build_single_tree(self, lig_dict, fp):
        """
        Build a single tree with fingerprint function

        :param lig_dict: all ligand information
        :param fp: fingerprint object
        :return: dot filename
        """
        distfile = self.gen_dist_file(lig_dict, fp.fp_func)
        newick_o = self.run_rapidnj(distfile)
        dot_inf = self.write_dotfile(newick_o)
        dot_out = self.sfdp_dot(dot_inf, 10)
        dot_dict = self.dot2dict(dot_out)
        return dot_dict

    @staticmethod
    def parse_lig_file(in_file, identifier):
        """
        parse ligand file and return a dictionary with identifier as IDs

        :param in_file: input file directory
        :param identifier: name for the identifier
        :return: a dictionray with ligand information
        """
        return ParseLigandFile(in_file, identifier)

    @staticmethod
    def gen_dist_file(liganddict, fp_func):
        """
        generate distance file which is the input of rapidnj program.

        :param liganddict: ligand information
        :param fp_func: fingerprint function
        :return: filename for distance file
        """
        smile_list = [ [lig_name, liganddict[lig_name][SMILE_COLUMNNAME]] for lig_name in liganddict.keys()]
        print "finish smile list"
        filename   = WriteAsPHYLIPFormat(smile_list, fp_func)
        print "finish writing phyli file"
        return filename

    def run_rapidnj(self, distance_file):
        """
        run rapidnj program on distance_file

        :param distance_file: directory of distance file
        :return: newick string
        """
        proc = subprocess.Popen([self._RAPIDNJ_COMMAND, distance_file, "-i", "pd"], stdout=subprocess.PIPE)
        newick = proc.stdout.read()
        return newick

    @staticmethod
    def write_dotfile(newick):
        """
        write newick string as dot file

        :param newick: newick string
        :return: dot file
        """
        return WriteDotFile(newick)

    def sfdp_dot(self, dot_infile, size):
        """
        run sdfp on dot file

        :param dot_infile: directory for dot file
        :param size: parameter for the sfdp
        :return: new filename
        """
        fmt= self._FILE_FORMAT + '_sfdp.gv'
        newfilename = datetime.datetime.now().strftime(fmt)
        if os.path.isfile(newfilename):
            os.remove(newfilename)
        command = "sfdp -Gsmoothing=triangle -Gsize={size} {infile} > {outfile}".format(size=size, infile=dot_infile, outfile=newfilename)
        subprocess.Popen( command, shell = True, stdout = subprocess.PIPE ).communicate()
        RemoveBackSlash(newfilename)
        return newfilename

    @staticmethod
    def dot2dict(dot_outfile):
        return Dot2Dict(dot_outfile, None)

    @staticmethod
    def gen_properties(ligand_dict, activities, properties, ext_cols):
        """
        Generate properties for each molecule.

        :param ligand_dict: ligand dictionary which keep all ligand information
        :param activities: a list of PropertyType objects
        :param properties: a list of PropertyType objects
        :param ext_cols: the column name for external links
        :return:
        """
        compounds = []
        for idx in range(len(ligand_dict)):
            lid = "B" + str(idx)
            comp = dict()
            comp["id"] = lid
            comp["orig_id"] = ligand_dict[lid]["orig_id"]
            comp["activities"] = dict()
            comp["properties"] = dict()
            comp["external"] = dict()
            for act in activities:
                comp["activities"][act.name] = act.gen_property(ligand_dict[lid])
            for prop in properties:
                comp["properties"][prop.name] = prop.gen_property(ligand_dict[lid])
            for col in ext_cols:
                ext_val = ligand_dict[lid][col]
                if isinstance(ext_val, float):
                    comp[col] = str(int(ext_val))
                else:
                    comp[col] = str(ext_val)
            compounds.append(comp)

        return compounds

    @staticmethod
    def make_structures_for_smiles( ligand_dict ):
        """
        Make structure figures from smile strings. All image files will be in the IMG_DIR

        :param ligand_dict: ligand dictionary which keep all ligand information
        :return:
        """
        relative_dir = IMG_DIR
        for key in ligand_dict:
            smile = ligand_dict[key][ SMILE_COLUMNNAME ]
            filename = ligand_dict[ key ][ "orig_id" ]
            mol = Chem.MolFromSmiles(smile)
            MolToFile( mol, os.path.join(relative_dir, '{}.svg'.format(filename)) )
