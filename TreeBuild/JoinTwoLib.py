'''
    Joining two molecule libraries
'''
import pickle

def JoinTwoFPList(fplist1, fplist2):
    pass


def JoinTwoLib(infile1, infile2):
    fp1list = pickle.load(open(infile1))
    fp2list = pickle.load(open(infile2))
    infile3 = JoinTwoFPList(fp1list, fp2list)
    return infile3

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Join two libs.')
    parser.add_argument('--infile1')
    parser.add_argument('--infile2')
    arg = parser.parse_args()
    JoinTwoLib(arg.infile1, arg.infile2)
