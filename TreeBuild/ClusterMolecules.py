'''
    Clustering Large Scale Molecule Library
'''

from sklearn.cluster import DBSCAN, MiniBatchKMeans
from numpy import arange
from DistanceMeasure import getSimilarity
from functools import partial

import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy
import pickle

from rdkit import Chem
from rdkit import DataStructs

def GetDistance(fplist, x, y):
    fp1 = fplist[int(x[0])][1]
    fp2 = fplist[int(y[0])][1]
    return 1 - getSimilarity(fp1, fp2)

def GetCenterOfCluster(fp_list, idx):
    '''
        Get the center of cluster, currently just return the item with median index
    '''
    clusters = []
    for i in idx:
        clusters.append(fp_list[i])
    return [clusters[len(clusters) / 2][0], len(idx)]

def GetSmallSetAfterClustering(fp_list):
    #dbs = DBSCAN(eps = 0.15, min_samples = len(fp_list) / 100000, metric = partial(GetDistance, fp_list), algorithm='ball_tree')
    dbs = DBSCAN(eps = 0.45, min_samples = 2, metric = partial(GetDistance, fp_list), algorithm='ball_tree')
    dbs.fit(arange(len(fp_list)).reshape(-1, 1))
    gp_name  = set(dbs.labels_)
    labels   = dbs.labels_
    center_list= []
    for each_gp in gp_name:
        if each_gp == -1:
            continue
        class_member_idx = np.where(labels == each_gp)[0]
        center = GetCenterOfCluster(fp_list, class_member_idx)
        center_list.append(center)
    #print labels
    print len(center_list)
    return center_list

def WriteToFile(center_list, ligand_dict, filename):
    fileobj = open(filename, "w")
    # write header
    first_e = ligand_dict.itervalues().next()
    fileobj.write("\t".join(["mapped_id", "size"] + first_e.keys()) + "\n")
    for lig in center_list:
        content = "\t".join(map(str, lig + ligand_dict[lig[0]].values()))
        fileobj.write(content + "\n")
    fileobj.close()

def Convert2Numpy(fp_list):
    numpy_list = []
    fold_times = 8
    print "size of zero matrix", len(fp_list), len(fp_list[0][1]) / fold_times
    print (len(fp_list), len(fp_list[0][1]) / fold_times)
    fp_matrix  = np.zeros((len(fp_list), len(fp_list[0][1]) / fold_times))
    for idx in range(len(fp_list)):
        fp  = fp_list[idx][1]
        fpa = np.zeros((1,len(fp)),np.double)
        DataStructs.ConvertToNumpyArray(fp, fpa)
        fp_matrix[idx,:] = np.sum(np.reshape(fpa, (fold_times, -1)), axis = 0)
    return fp_matrix

def GetCenterOfClusterKMeans(fp_list, fp_matrix, center_fp, idx):
    '''
        Get the center of cluster, currently just return the item with median index
    '''
    min_dist = 100000
    center_idx = None
    for i in idx:
        #print center_fp
        dist = scipy.spatial.distance.euclidean(fp_matrix[i,:], center_fp)
        if dist < min_dist:
            min_dist = dist
            center_idx = i
    return [fp_list[center_idx][0], len(idx)]


def GetSmallSetAfterClusteringKMeans(fp_list, cluster_n):
    print cluster_n
    fp_matrix = Convert2Numpy(fp_list)
    mnb = MiniBatchKMeans(n_clusters = cluster_n)
    #mnb = MiniBatchKMeans(n_clusters = 2)
    mnb.fit(fp_matrix)
    centers = mnb.cluster_centers_
    labels   = mnb.labels_

    gp_name  = set(mnb.labels_)
    center_list= []
    for each_gp in gp_name:
        if each_gp == -1:
            continue
        class_member_idx = np.where(labels == each_gp)[0]
        center = GetCenterOfClusterKMeans(fp_list, fp_matrix, centers[each_gp, :], class_member_idx)
        center_list.append(center)
    return center_list

def ToFPObjRDK(alist):
    # for alist, the first item is ligand name, the second is smile
    newlist = []
    for each in alist:
        smile = each[1]
        m = Chem.MolFromSmiles(smile)
        if m is None:
            continue
        fp =  Chem.RDKFingerprint(m)
        newlist.append([each[0], fp])
    return newlist


if __name__ == "__main__":
    from Util import ParseLigandFile
    from Model import SMILE_COLUMNNAME, FILE_FORMAT
    from DistanceMeasure import ToFPObj
    import argparse
    parser = argparse.ArgumentParser(description='Get a clustered result.')
    parser.add_argument('--infile')
    parser.add_argument('--method')
    parser.add_argument('--outfile')
    parser.add_argument('--num_cluster', type=int)
    parser.add_argument('--use_cache')
    parser.add_argument('--save_cache')
    arg = parser.parse_args()
    infile =  arg.infile

    ligand_dict = ParseLigandFile(infile)
    #ligand_dict = ParseLigandFile("./Data/result_clean_no0.txt")
    #ligand_dict = ParseLigandFile("./Data/chembl_20_10000.txt")
    smile_list = [ [lig_name, ligand_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in ligand_dict.keys()]
    if arg.method == "KMeans":
        if arg.use_cache:
            with open(arg.use_cache, "rb") as inputobj:
                fp_list    = pickle.load(inputobj)
        else:
            fp_list    = ToFPObjRDK(smile_list)
            with open(arg.save_cache, "wb") as outputobj:
                pickle.dump(fp_list, outputobj, pickle.HIGHEST_PROTOCOL)
        center_list= GetSmallSetAfterClusteringKMeans(fp_list, arg.num_cluster)
    else:
        fp_list    = ToFPObj(smile_list)
        center_list= GetSmallSetAfterClustering(fp_list)

    WriteToFile(center_list, ligand_dict, arg.outfile)

