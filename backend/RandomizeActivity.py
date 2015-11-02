'''
    Randomlize the activity column
'''

from Model import POTENCY
import argparse
from csv import DictReader, DictWriter

import random

def RandomizeActivity(infile):
    act_list = []
    content = []
    reader = DictReader(open(infile), delimiter = "\t")
    for line in reader:
        content.append(line)
        act_list.append(line[POTENCY])
    random.shuffle(act_list)
    w = DictWriter(open(infile + "_rand", "w"), fieldnames = reader.fieldnames, delimiter = "\t" )
    w.writeheader()
    for idx, line in enumerate(content):
        line[POTENCY] = act_list[idx]
        w.writerow(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')
    arg = parser.parse_args()

    infile = arg.infile
    RandomizeActivity(infile)
