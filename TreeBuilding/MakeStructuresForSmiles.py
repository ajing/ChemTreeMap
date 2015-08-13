"""
  Using pybel to create image for each smile string
"""

import pybel
from  FileParser import GetAllinfo

def MakeStructuresForSmiles( all_info ):
    totalrow = len( all_info[ "ligandid" ] )
    relativedir = "./Data/"
    for index in range( totalrow ):
        smile = all_info[ "Canonical_Smiles" ][ index ]
        pngfile = all_info[ "ligandid" ][ index ]
        mol = pybel.readstring( "smi", smile )
        mol.draw( show=False, filename = relativedir + pngfile )


if __name__ == "__main__":
    infile = "./Data/result_clean.txt"
    all_info = GetAllinfo( infile )
    MakeStructuresForSmiles( all_info )
