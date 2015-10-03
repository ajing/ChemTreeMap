"""
  Using pybel to create image for each smile string
"""

import pybel
from Util import ParseLigandFile
import argparse

from rdkit import Chem
from rdkit.Chem.Draw import MolToFile

import os

from Model import SMILE_COLUMNNAME, IDENTIFIER, IMG_DIR

def MakeStructuresForSmiles( ligand_dict ):
    relativedir = IMG_DIR
    for key in ligand_dict:
        smile = ligand_dict[key][ SMILE_COLUMNNAME ]
        filename = ligand_dict[ key ][ "orig_id" ]
        mol = Chem.MolFromSmiles(smile)
        print filename
        try:
            MolToFile( mol, os.path.join(relativedir, '{}.svg'.format(filename)) )
            print "successfully print to file"
        except:
            continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot images for each molecule.')
    parser.add_argument('--infile')
    arg = parser.parse_args()
    all_info = ParseLigandFile( arg.infile )
    MakeStructuresForSmiles( all_info )
