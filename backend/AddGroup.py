'''
    Add color group for each node
'''


def AddColorGroup(infile):
    outfile = infile + "_new"
    outobj  = open(outfile, "w")
    flag = 1
    for line in open(infile):
        if flag:
            newline = line.strip() + "\tgroup\n"
            flag = 0
        else:
            if not line.strip().split("\t")[0].startswith("CHEMBL"):
                newline = line.strip() + "\tred\n"
            else:
                newline = line.strip() + "\tblue\n"
        outobj.write(newline)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Add group to file.')
    parser.add_argument('--infile')
    arg = parser.parse_args()

    infile = arg.infile
    AddColorGroup(infile)
