'''
    Build Tree from dmatrix
'''
import datetime
from TreeConstruction import *
import pickle

class BuildTree():
    def __init__( self, leaderlist, dmatrix, moldict, figurename):
        self.distanceMatrix = dmatrix
        self.leaderList = leaderlist
        self.molDict  = moldict
        self.figuresize = len(leaderlist)
        self.figure   = figurename + ".svg"
        self.savePath = "./Data/"
        self.imgPath = "./Image/"
        self.size    = len(leaderlist)
        self.prepareNJ()

    def LeaderName( self, leaderlist ):
        leadername = []
        for each in leaderlist:
            leadername.append(self.molDict[each]["ligandid"])
        return leadername

    def LeaderMatrix( self, dmatrix, leaderlist ):
        return dmatrix[leaderlist, :][:, leaderlist]

    def prepareNJ(self):
        leader_name = self.LeaderName(self.leaderList)
        d_shrink_matrix = self.LeaderMatrix( self.distanceMatrix, self.leaderList )
        if not d_shrink_matrix.size == len(leader_name)**2:
            raise ValueError("size doesn't match between leaderlist and distance matrix")
        DMatrix = DistanceMatrix(leader_name, d_shrink_matrix)
        root    = nj(DMatrix, self.molDict)
        fmt='%Y-%m-%d-%Hh-%Mm_{fname}_pickle_tree'
        newfilename = datetime.datetime.now().strftime(fmt).format(fname = self.figure)
        pickle.dump(root, open(newfilename, "w"))
        return root

if __name__ == "__main__":
    pass
