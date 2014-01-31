'''
    Parse affinity file, which contains
'''

def ParseAffinity(infile):
    smile_dict = dict()
    affin_dict = dict()
    for line in open(infile):
        name, smile, affinity = line.split()
        affinity= float(affinity)
        smile_dict[name] = smile
        affin_dict[name] = affinity
    return smile_dict, affin_dict

if __name__ == "__main__":
    infile = "affinity.txt"
    smile, affin  = ParseAffinity(infile)
    print smile, affin

