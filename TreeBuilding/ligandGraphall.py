#####################
## Author : ajing
## Date   : 5/14/2013
## Purpose: create a graph from a list of SMILES strings
## Design : Using RDKit to get ECFP6 Taminoto similarity score
## Function: parseLigandFile: read ligand file
##           getSimilarity:   calc similarity between two smiles
##           pairwiseSimilarity: calc similarity for ligand dict pairwise
##           printPairwise: print output from pairwiseSimilarity
######################

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs

# for similarity matrix
import numpy as np

# for clean of smile string
import openbabel

def CleanSmile( smile ):
    # only keep smile string with more than 6 atoms
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smile)
    mol.StripSalts()
    if mol.NumAtoms() < 6:
        return None
    clean_smi = obConversion.WriteString(mol)
    return clean_smi

def get_allinfo(infile):
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

def parseLigandFile(filename, bindingtype = None):
    # parse ligand file to ligand id and smile string as dict obj
    ligand_dict = dict()
    all_info = get_allinfo(filename)
    for index in range(len(all_info["ligandid"])):
        smile = all_info["Canonical_Smiles"][index]
        smile = CleanSmile(smile)
        if smile is None:
            continue
        if bindingtype is None:
            ligand_dict[all_info["ligandid"][index]] = smile
        else:
            if all_info["typeofbinding"][index] == bindingtype:
                ligand_dict[all_info["ligandid"][index]] = smile
    print "ligand_dict length", len(ligand_dict.keys())
    return ligand_dict


def getSimilarity(smile1, smile2):
    # generate similarity score for two smiles strings
    m1 = Chem.MolFromSmiles(smile1)
    m2 = Chem.MolFromSmiles(smile2)
    if (m1 is None or m2 is None):
        return
    # ECFP6
    fp1 = AllChem.GetMorganFingerprint(m1, 3)
    fp2 = AllChem.GetMorganFingerprint(m2, 3)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def getSimilarityAtomPair(smile1, smile2):
    # generate similarity score for two smiles strings
    m1 = Chem.MolFromSmiles(smile1)
    m2 = Chem.MolFromSmiles(smile2)
    if (m1 is None or m2 is None):
        return
    # atom-pair descriptors
    fp1 = Pairs.GetAtomPairFingerprint(m1)
    fp2 = Pairs.GetAtomPairFingerprint(m2)
    return DataStructs.DiceSimilarity(fp1, fp2)

def pairwiseSimilarity(ligand_dict, simfunc):
    #calc similarity for ligand dict pairwise
    pairwise = []
    all_ligand_names = ligand_dict.keys()
    num_ligand = len(all_ligand_names)
    print num_ligand
    for i in range(num_ligand - 1):
        eachid1 = all_ligand_names[i]
        for j in range(i + 1, num_ligand):
            eachid2 = all_ligand_names[j]
            #print "eachid1"
            #print ligand_dict[eachid1]
            #print "eachid2"
            #print ligand_dict[eachid2]
            pairwise.append( [ eachid1, eachid2, simfunc( ligand_dict[eachid1], ligand_dict[eachid2] )] )
        print i
    return pairwise

def printPairwise(pairwise, outfilename):
    outfileobj = open(outfilename, 'w')
    for each in pairwise:
        outfileobj.write( "\t".join(each[:2]) + "\t%.5g" % each[2] + "\n")

def filterOutput(infile, newfile = "newfile.pair"):
    out_obj = open(newfile, "w")
    for line in open(infile):
        content = line.split("\t")
        if float(content[2]) > 0.9:
            out_obj.write(line)
    out_obj.close()

def similarityMatrix(ligand_dict, simfunc):
    # return matrix for similarity between different ligands
    all_ligand_names = ligand_dict.keys()
    num_ligand = len(all_ligand_names)
    print num_ligand
    simmatrix  = np.zeros((num_ligand,num_ligand))
    for i in range(num_ligand):
        print i
        eachid1 = all_ligand_names[i]
        for j in range( i, num_ligand ):
            eachid2 = all_ligand_names[j]
            simmatrix[i,j] = simfunc( ligand_dict[eachid1], ligand_dict[eachid2] )
            simmatrix[j,i] = simfunc( ligand_dict[eachid1], ligand_dict[eachid2] )
    #print simmatrix
    return simmatrix

def NewLigandFile( ligand_dict, filename ):
    # new ligand file with non redundant information
    newfilename = filename + "_new"
    out_obj     = open( newfilename, "w" )
    in_lines    = open( filename, "r" ).readlines()
    for eachLigandID in ligand_dict:
        #print eachLigandID
        for eachline in in_lines:
            content = eachline.split("\t")
            #print "content:", content[1]
            if content[1] == eachLigandID:
                out_obj.write(eachline)
                break
    out_obj.close()
    print "finish writing to new file"

if __name__ == "__main__":
    # test case
    ligandfile = "./Data/ligand_5_7_ppilot.txt"
    liganddict = parseLigandFile(ligandfile)
    NewLigandFile( liganddict, ligandfile )
    #pairwise = pairwiseSimilarity(liganddict, getSimilarity)
    smatrix  = similarityMatrix(liganddict, getSimilarity)
    np.save("similarityMatrix", smatrix)
    #printPairwise(pairwise, "test.pair")
    #filterOutput(ligandfile + "_pair", ligandfile + "_pair_filter")
