'''
    Joining two molecule libraries by finding unique elements first, clustering three libraries seperately and then join them together
'''

from Util import ParseLigandFile
from DistanceMeasure import getSimilarity
from JoinTwoLib2 import MergeTwoDicts
from ClusterMolecules import GetSmallSetAfterClusteringKMeans, WriteToFile

import pickle
import argparse
from csv import DictWriter

def JoinTwoLib(infile1, infile2, fpfile1, fpfile2, outfile):
    group1 = "blue"
    group2 = "red"
    group12= "black"
    ligand_dict1 = ParseLigandFile(infile1)
    ligand_dict2 = ParseLigandFile(infile2)
    ligand_dict_k1 = { value["orig_id"]: value for key, value in ligand_dict1.iteritems()}
    ligand_dict_k2 = { value["orig_id"]: value for key, value in ligand_dict2.iteritems()}
    fp1_list     = pickle.load(open(fpfile1))
    fp2_list     = pickle.load(open(fpfile2))

    fp1_dict     = dict()
    fp2_dict     = dict()
    fp12_dict    = dict()

    for idx in range(len(fp1_list)):
        #print "fp1_list, ", idx
        fp1 = fp1_list[idx]
        try:
            fp1_dict[fp1[1].ToBitString()]
        except:
            fp1_dict[fp1[1].ToBitString()] = fp1

    for idx in range(len(fp2_list)):
        #print "fp1_list, ", idx
        fp2 = fp2_list[idx]
        try:
            fp2_dict[fp2[1].ToBitString()]
        except:
            fp2_dict[fp2[1].ToBitString()] = fp2

    for key2 in fp2_dict.keys():
        try:
            fp1_dict[key2]
            fp12_dict[key2] = fp2_dict[key2]
            del fp1_dict[key2]
            del fp2_dict[key2]
        except:
            continue

    #writer = DictWriter(open(outfile, "w"), fieldnames=fieldnames, delimiter = "\t")
    #writer.writeheader()
    #for item in newlist:
    #    writer.writerow(item)

    shrinkage = 0.005
    fp1_new_list = fp1_dict.values()
    fp2_new_list = fp2_dict.values()
    fp12_new_list = fp12_dict.values()
    print len(fp1_new_list)
    print len(fp2_new_list)
    print len(fp12_new_list)
    num_cluster1 = int(len(fp1_new_list) * shrinkage)
    num_cluster2 = int(len(fp2_new_list) * shrinkage)
    num_cluster12 = int(len(fp12_new_list) * shrinkage)
    print num_cluster1, num_cluster2, num_cluster12
    center_list1 = GetSmallSetAfterClusteringKMeans(fp1_new_list, num_cluster1)
    center_list2 = GetSmallSetAfterClusteringKMeans(fp2_new_list, num_cluster2)
    center_list12 = GetSmallSetAfterClusteringKMeans(fp12_new_list, num_cluster12)
    print len(center_list1)
    print len(center_list2)
    print len(center_list12)
    WriteToFile(center_list1, ligand_dict_k1, outfile, "w", group1)
    WriteToFile(center_list2, ligand_dict_k2, outfile, "a", group2)
    WriteToFile(center_list12, ligand_dict_k2, outfile, "a", group12)
    #WriteToFile(center_list12, ligand_dict_k1, outfile, "w", group12)

def main():
    parser = argparse.ArgumentParser(description='join two molecule libraries.')
    parser.add_argument('--raw1')
    parser.add_argument('--raw2')
    parser.add_argument('--fp1')
    parser.add_argument('--fp2')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    JoinTwoLib(arg.raw1, arg.raw2, arg.fp1, arg.fp2, arg.outfile)

if __name__ == "__main__":
    main()
