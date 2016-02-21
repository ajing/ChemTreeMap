#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods
import datetime
import os
from rdkit import Chem
from rdkit.Chem.Draw import MolToFile

from types import FingerPrintType
from util import ParseLigandFile, WriteJSON, WriteAsPHYLIPFormat, Dot2Dict, \
    WriteDotFile, RemoveBackSlash
from model import IMG_DIR, SMILE_COLUMNNAME, RAPIDNJ_COMMAND, FILE_FORMAT
import subprocess

class TreeBuild:
    """
    There are a few assumptions for the input file:
        1. potency unit is nM
    """
    def __init__(self, input_file, output_file, id_column, fps, properties):

        self._RAPIDNJ_COMMAND = RAPIDNJ_COMMAND
        self._FILE_FORMAT = FILE_FORMAT

        activities = properties["activities"]
        other_properties = properties["others"]
        lig_dict = self.parse_lig_file(input_file, id_column)
        trees = dict()
        for fp in fps:
            assert isinstance(fp, FingerPrintType)
            trees[fp.name] = self.build_single_tree(lig_dict, fp)
        metadata = dict()
        metadata["activityType"] = [act.to_dict() for act in activities]
        metadata["treeTypes"] = [fp.to_dict() for fp in fps]
        metadata["circleSizeTypes"] = [prop.to_dict() for prop in other_properties]
        metadata["circleBorderTypes"] = [prop.to_dict() for prop in other_properties]

        comp_info = self.gen_properties(lig_dict, activities, other_properties)
        final_dict = {"metadata": metadata, "trees": trees, "compounds": comp_info}

        WriteJSON(final_dict, outfile=output_file, write_type="w")
         # make image file
        self.make_structures_for_smiles(lig_dict)

    def build_single_tree(self, lig_dict, fp):
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
    def gen_properties(ligand_dict, activities, properties):
        """
        Generate properties for each molecule.
        :param ligand_dict: ligand dictionary which keep all ligand information
        :param activities: a list of PropertyType objects
        :param properties: a list of PropertyType objects
        :return:
        """
        compounds = []
        for lid in ligand_dict:
            comp = dict()
            comp["id"] = lid
            comp["orig_id"] = ligand_dict[lid]["orig_id"]
            comp["activities"] = dict()
            for act in activities:
                comp["activities"][act.name] = act.gen_property(ligand_dict[lid])
            for prop in properties:
                comp["properties"][prop.name] = prop.gen_property(ligand_dict[lid])
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
            print filename
            try:
                MolToFile( mol, os.path.join(relative_dir, '{}.svg'.format(filename)) )
                print "successfully print to file"
            except:
                raise Exception("cannot write to file: " + os.path.join(relative_dir, '{}.svg'.format(filename)))