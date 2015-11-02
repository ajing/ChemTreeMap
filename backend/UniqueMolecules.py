'''
    Get the number of unique molecules
'''
import sys
import pickle

def UniqueMolecules(fp_list):
    print len(fp_list)
    print dir(fp_list[0][1])
    print len(set([x[1].ToBitString() for x in fp_list ]))

def UniqueMolecules2(fp_list1, fp_list2):
    fp_list = [x[1].ToBitString() for x in fp_list1] + [x[1].ToBitString() for x in fp_list2 ]
    #print fp_list
    print len(fp_list)
    print len(set(fp_list))

def Unique1():
    inputfile = sys.argv[1]
    inputobj  = open(inputfile, "rb")
    fp_list   = pickle.load(inputobj)
    UniqueMolecules(fp_list)

def Unique2():
    inputfile1 = sys.argv[1]
    inputobj1  = open(inputfile1, "rb")
    inputfile2 = sys.argv[2]
    inputobj2  = open(inputfile2, "rb")
    fp_list1   = pickle.load(inputobj1)
    fp_list2   = pickle.load(inputobj2)
    UniqueMolecules2(fp_list1, fp_list2)

if __name__ == "__main__":
    Unique2()
