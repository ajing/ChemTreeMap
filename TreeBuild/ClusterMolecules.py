'''
    Clustering Large Scale Molecule Library
'''

from sklearn.cluster import DBSCAN, dbscan
from numpy import arange
from DistanceMeasure import getSimilarity
from functools import partial

def GetDistance(fplist, x, y):
    fp1 = fplist[int(x[0])]
    fp2 = fplist[int(y[0])]
    return 1 - getSimilarity(fp1, fp2)

def GetSmallSetAfterClustering(fp_list):
    dbs = DBSCAN(eps = 0.8, min_samples = 2, metric = partial(GetDistance, fp_list))
    dbs.fit(arange(len(fp_list)).reshape(-1, 1))
    core_idx = db.core_sample_indices_
    return core_idx

def WriteToFile()


if __name__ == "__main__":
    from Util import ParseLigandFile
    from Model import SMILE_COLUMNNAME, FILE_FORMAT
    from DistanceMeasure import ToFPObj
    ligand_dict = ParseLigandFile("./Data/result_clean_no0_20.txt")
    smile_list = [ [lig_name, ligand_dict[lig_name][SMILE_COLUMNNAME]] for lig_name in ligand_dict.keys()]
    fp_list    = ToFPObj(smile_list)
    GetSmallSetAfterClustering(fp_list)
