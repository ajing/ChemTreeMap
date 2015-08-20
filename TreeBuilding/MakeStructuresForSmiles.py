"""
  Using pybel to create image for each smile string
"""

import pybel
from  FileParser import GetAllinfo
import argparse

def MakeStructuresForSmiles( all_info ):
    totalrow = len( all_info[ "ligandid" ] )
    relativedir = "./Data/Image/"
    for index in range( totalrow ):
        smile = all_info[ "Canonical_Smiles" ][ index ]
        pngfile = all_info[ "ligandid" ][ index ]
        mol = pybel.readstring( "smi", smile )
        mol.draw( show=False, filename = relativedir + pngfile )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot images for each molecule.')
    parser.add_argument('--infile')
    arg = parser.parse_args()
    all_info = GetAllinfo( arg.infile )
    MakeStructuresForSmiles( all_info )
