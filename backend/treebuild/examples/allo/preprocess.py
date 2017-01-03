
infile = "allo_comp.txt"
outfile = "allo_comp_mod.txt"
outobj = open(outfile, "w")

import csv

f_line = True

for line in open(infile):
    if f_line:
        f_line = False
        outobj.write(line)
    else:
        content = line.split("\t")
        content[2] = "allosteric" if content[2] == "1" else "competitive"
        outobj.write("\t".join(content))

outobj.close()
