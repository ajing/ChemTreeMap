"""
  Using pybel to create image for each smile string
"""

import pybel
from ParseFile import ParseAffinity

def MakeStructuresForSmiles( smile_dict ):
    relativedir = "./Image/"
    for name in smile_dict:
        pngfile = name
        smile   = smile_dict[name]
        filename = relativedir + pngfile
        print name, smile, filename
        mol = pybel.readstring( "smi", smile )
        mol.draw( show = False, filename = filename )


if __name__ == "__main__":
    infile = "./affinity.txt_new"
    s_dict, _ = ParseAffinity(infile)
    MakeStructuresForSmiles( s_dict )
