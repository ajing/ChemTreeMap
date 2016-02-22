#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods

class GenerateIDs:
    """
    Retrieve other external ids for BindingDB IDs.
    """
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

        self.pubchem_dict = self.parse_dict_file("./data/BindingDB_CID.txt")
        self.chebi_dict = self.parse_dict_file("./data/BindingDB_CHEBI_ID.txt")

        self.cannot_find_pubchem = 0
        self.cannot_find_chebi   = 0
        self.run_for_file()

        print("# PubChem ID cannot find " + str(self.cannot_find_pubchem))
        print("# CHEBI ID cannot find " + str(self.cannot_find_chebi))


    @staticmethod
    def parse_dict_file(filename):
        id_dict = dict()
        with open(filename) as fobj:
            for line in fobj:
                cont = line.split()
                id_dict[cont[0]] = cont[1]
        return id_dict

    def get_pubchemid_from_bdid(self, bdid):
        try:
            return self.pubchem_dict[bdid]
        except:
            self.cannot_find_pubchem += 1
            return None

    def get_chebiid_from_bdid(self, bdid):
        try:
            return self.chebi_dict[bdid]
        except:
            self.cannot_find_chebi += 1
            return None

    def run_for_file(self):
        import csv
        out_obj = open(self.outfile, "w")
        with open(self.infile) as in_obj:
            reader = csv.DictReader(in_obj, delimiter = '\t')
            writer = csv.DictWriter(out_obj, reader.fieldnames + ["PubChem"], delimiter = '\t')
            writer.writeheader()
            for line in reader:
                line["PubChem"] = self.get_pubchemid_from_bdid(line["BindingDB"])
                #line["CHEBI"] = self.get_chebiid_from_bdid(line["BindingDB"])
                writer.writerow(line)

if __name__ == "__main__":
    import sys
    GenerateIDs(sys.argv[1], sys.argv[2])