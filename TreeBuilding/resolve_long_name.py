'''
 resolve long name problem
'''

from Dot2JSON import Dot2JSON, Root2JSON
from TreeRebuilder import RewriteDot
from CreateGraph import MoleculeDictionary
from TreefromSmile import Convert2NJmoldict

def first():
    # run the sfdp command, and then run the following code
    root = Dot2JSON("test_sfdp.gv")
    Root2JSON(root, "test.json")

def second():
    # rewrite the file
    parser = argparse.ArgumentParser(description='Simplify the graph.')
    parser.add_argument('--infile')
    arg = parser.parse_args()

    moldict = Convert2NJmoldict(MoleculeDictionary(arg.infile))

    newfilename, newmapdict = RewriteDot(arg.infile)
    sfdp_dot  = SFDPonDot(newfilename, 10)
    root = Dot2JSON(sfdp_dot)
    rootdict = RecursiveNode2Dict(root, moldict)
    RecursiveChangeName(rootdict, new2old)
    fileobj  = open(filename, "w")
    fileobj.write(json.dumps(rootdict, indent=2))

if __name__ == "__main__":
    second()

