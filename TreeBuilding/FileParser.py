'''
    Goodies to parse all kinds of files
'''
def GetAllinfo(infile):
    # To understand ligand cluster file
    rownum = 0
    all_info = dict()
    # first_col is only for descriptor file, the first column is a smile string
    first_col = dict()
    for line in open(infile):
        if rownum == 0:
            header = line.strip().split('\t')
            #print header
            for each in header:
                all_info[each] = []
        else:
            content = line.strip().split('\t')
            first_col[content[0]] = content[1:]
            for i in range(len(header)):
                try:
                    all_info[header[i]].append(content[i])
                except:
                    print content
        rownum += 1
    return all_info


