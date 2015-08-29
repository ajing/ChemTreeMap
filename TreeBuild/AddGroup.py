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
            if line.strip().split("\t")[0].startswith("B"):
                newline = line.strip() + "\tred\n"
            else:
                newline = line.strip() + "\tblue\n"
        outobj.write(newline)


if __name__ == "__main__":
    AddColorGroup("./Data/combine_chembl_bindingdb.txt")
