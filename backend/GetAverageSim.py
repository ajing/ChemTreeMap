"""
    Get average similarity
"""

import argparse

def GetAverageSim(infile):
    inobj = open(infile)
    firstline = True
    line_num  = 0
    all_sim   = []
    for line in inobj:
        content = line.strip()
        if firstline:
            firstline = False
            N = int(content)
            continue
        line_num += 1
        distlist = content.split("\t")[line_num + 1 : ]
        distlist = [1 - float(x) for x in distlist]
        all_sim  += distlist
    print sum(all_sim) / len(all_sim)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get average similarity.')
    parser.add_argument('--infile')
    arg = parser.parse_args()
    GetAverageSim(arg.infile)
