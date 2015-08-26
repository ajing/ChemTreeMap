'''
    Process PipelinePilot cluster output file
'''
from csv import DictReader

def GetCenterList(infile):
    # get the list of ligand name which are the center and the cluster size
    lig_name_col = "Column1"
    center_list  = []
    ig_count     = 0
    for line in DictReader(open(infile), delimiter = "\t"):
        if ig_count != 0:
            ig_count = ig_count - 1
            continue
        print line
        if float(line["DistanceToClosest"]) == 0:
            cluster_size= int(line["ClusterSize"])
            center_list.append([line[lig_name_col], cluster_size])
            ig_count = cluster_size - 1
    print center_list


if __name__ == "__main__":
    import  argparse
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--infile')

    arg = parser.parse_args()
    infile = arg.infile
    GetCenterList(infile)
