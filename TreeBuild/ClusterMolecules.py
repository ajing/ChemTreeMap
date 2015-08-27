'''
    Clustering Large Scale Molecule Library
'''

from sklearn.cluster import DBSCAN, MiniBatchKMeans
from numpy import arange
import numpy as np
from DistanceMeasure import getSimilarity
from functools import partial

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
    return clusters[len(clusters) / 2]

def GetSmallSetAfterClustering(fp_list):
    dbs = DBSCAN(eps = 0.8, min_samples = 2, metric = partial(GetDistance, fp_list), algorithm='ball_tree')
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
    print labels
    print len(center_list)
    return center_list

def GetSmallSetAfterClusteringKMeans(fp_list):
    #mnb = MiniBatchKMeans(n_clusters = 10000, metric = partial(GetDistance, fp_list))
    #mnb = MiniBatchKMeans(n_clusters = 2, metric = partial(GetDistance, fp_list))
    #mnb.fit(arange(len(fp_list)).reshape(-1, 1))
    #center_idx = mnb.cluster_centers_
    #mnb.fit(arange(len(fp_list)).reshape(-1, 1))
    #center_idx = mnb.cluster_centers_
    pass


def WriteToFile():
    pass


if __name__ == "__main__":
    from Util import ParseLigandFile
    from Model import SMILE_COLUMNNAME, FILE_FORMAT
    from DistanceMeasure import ToFPObj
    #ligand_dict = ParseLigandFile("./Data/chembl_20_new.txt")
    #ligand_dict = ParseLigandFile("./Data/result_clean_no0.txt")
    ligand_dict = ParseLigandFile("./Data/chembl_20_10000.txt")
    smile_list = [ [lig_name, ligand_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in ligand_dict.keys()]
    fp_list    = ToFPObj(smile_list)
    GetSmallSetAfterClustering(fp_list)
