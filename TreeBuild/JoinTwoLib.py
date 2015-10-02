'''
    Joining two molecule libraries
'''
import pickle
import numpy as np

from ClusterMolecules import Convert2Numpy, GetCenterOfClusterKMeans
from sklearn.cluster import MiniBatchKMeans

def AssignGroup(class_member_idx, fp_list1_size):
    group1 = np.any(class_member_idx < fp_list1_size)
    group2 = np.any(class_member_idx >= fp_list1_size)
    if group1 and group2:
        #return "Mix"
        return "black"
    elif group1:
        #return "Group1"
        return "blue"
    else:
        #return "Group2"
        return "red"

def WriteToFile(center_list, gp_list, ligand_dict, filename):
    fileobj = open(filename, "w")
    # write header
    first_e = ligand_dict.itervalues().next()
    fileobj.write("\t".join(["mapped_id", "size", "group"] + first_e.keys()) + "\n")
    for idx in range(len(center_list)):
        lig = center_list[idx]
        gp  = gp_list[idx]
        content = "\t".join(map(str, lig + [gp] + ligand_dict[lig[0]].values()))
        fileobj.write(content + "\n")
    fileobj.close()

def GetSmallSetAfterClusteringKMeans(fp_list1, fp_list2, cluster_n):
    fp1_len   = len(fp_list1)
    fp_matrix = Convert2Numpy(fp_list1 + fp_list2)
    mnb = MiniBatchKMeans(n_clusters = cluster_n)
    mnb.fit(fp_matrix)
    centers = mnb.cluster_centers_
    labels   = mnb.labels_

    gp_name  = set(mnb.labels_)
    center_list= []
    gp_list    = []
    for each_gp in gp_name:
        if each_gp == -1:
            continue
        class_member_idx = np.where(labels == each_gp)[0]
        center = GetCenterOfClusterKMeans(fp_list, fp_matrix, centers[each_gp, :], class_member_idx)
        center_list.append(center)
        gp_list.append(AssignGroup(class_member_idx, fp1_len))
    return center_list, gp_list

def MergeTwoDicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

if __name__ == "__main__":
    import argparse
    from Util import ParseLigandFile
    parser = argparse.ArgumentParser(description='Join two libs.')
    parser.add_argument('--raw1')
    parser.add_argument('--raw2')
    parser.add_argument('--fp1')
    parser.add_argument('--fp2')
    parser.add_argument('--cluster_n')
    parser.add_argument('--outfile')
    arg = parser.parse_args()

    ligand_dict1 = ParseLigandFile(arg.raw1)
    ligand_dict2 = ParseLigandFile(arg.raw2)
    fp1list = pickle.load(open(arg.fp1))
    fp2list = pickle.load(open(arg.fp2))

    center_list, gp_list = GetSmallSetAfterClusteringKMeans(fp1list, fp2list, arg.cluster_n)
    WriteToFile(center_list, gp_list, MergeTwoDicts(ligand_dict1, ligand_dict2), arg.outfile)
