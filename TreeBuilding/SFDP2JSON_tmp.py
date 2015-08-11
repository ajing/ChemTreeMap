'''
    raw Dot file without position info to JSON file
'''

from SFDPLayOut import SFDPonDot
from Dot2JSON import Dot2JSON, Root2JSON

def main():
    infile = "./Data/all_0.9.gv"
    outfile = "./Data/all_0.9.json"
    sfdp_dot  = SFDPonDot(infile, 10)
    root = Dot2JSON(sfdp_dot)
    Root2JSON(root, outfile)

if __name__ == "__main__":
    main()
