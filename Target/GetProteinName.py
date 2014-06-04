'''
    Get ligand protein targets
    I know this can be highly redundant, but at this stage I want to impliment this function
    @Date: 5/18/2014
    @Changelog: add chiral data to tree 6/2/2014
'''

import csv
import os
import json

__FILEDIR__ = "../Data/"

def GetProteinName(prolist, proid):
    for each in prolist:
        if each['proteinid'] == proid:
            return each['description']

def GetLigandTargetName(ligandfile, proteinfile):
    procontent = []
    with open(proteinfile, 'rb') as proobj:
        proreader = csv.DictReader(proobj, delimiter="\t")
        for row in proreader:
            procontent.append(row)
    ligdict = dict()
    with open(ligandfile, 'rb') as ligobj:
        ligreader = csv.DictReader(ligobj, delimiter="\t")
        for row in ligreader:
            ligdict[row['ligandid']] = GetProteinName(procontent, row['proteinid'])
    return ligdict

def TopN(ligdict, N):
    prolist = [ x for x in ligdict.values()]
    proset  = set(prolist)
    prodict = dict()
    for each in proset:
        prodict[each] = prolist.count(each)
    popularpro = sorted(prodict, key = prodict.get, reverse = True)
    return popularpro[:N]

def TopNBlue(ligdict, N, ligandfile):
    ligtype = dict()
    with open(ligandfile, 'rb') as ligobj:
        ligreader = csv.DictReader(ligobj, delimiter="\t")
        for row in ligreader:
            ligtype[row['ligandid']] = row['typeofbinding']
    prolist = [ v for k, v in ligdict.items() if ligtype[k] == "competitive"]
    proset  = set(prolist)
    prodict = dict()
    for each in proset:
        prodict[each] = prolist.count(each)
    popularpro = sorted(prodict, key = prodict.get, reverse = True)
    return popularpro[:N]

def AddMemoization(valname, value):
    def inner(func):
        setattr(func, valname, value)
        return func
    return inner

class NewJSON:
    def __init__(self, jsondict):
        self.json = jsondict

    @AddMemoization("lignames", [])
    def AddField2JSON(self, jsondict, ligdict, fieldname):
        lignames = self.AddField2JSON.lignames
        if not lignames:
            lignames += ligdict.keys()
        if "children" in jsondict:
            for child in jsondict["children"]:
                self.AddField2JSON(child, ligdict, fieldname)
        else:
            if jsondict["name"] in lignames:
                jsondict[fieldname] = ligdict[jsondict["name"]]

    def AddField(self, ligdict, fieldname):
        self.AddField2JSON(self.json, ligdict, fieldname)

    def EmptyMemorization(self):
        self.AddField2JSON.lignames = []

def GetProteinTargetMain():
    filedir = __FILEDIR__
    ligdict = GetLigandTargetName(os.path.join(filedir, "ligand_5_7.txt"), os.path.join(filedir, "proteinseq_5_4.txt"))
    print TopN(ligdict, 12)
    #print TopNBlue(ligdict, 12, os.path.join(filedir, "ligand_5_7.txt"))
    #jdict = json.load(open("all_0.9.json"))
    jdict = json.load(open("test.json"))
    newjson = NewJSON(jdict)
    newjson.AddField(ligdict, "target")
    fileobj  = open("withtarget_small.json", "w")
    fileobj.write(json.dumps(jdict, indent=2))


# The following functions for getting ligand descriptors
def GetDescValue(desclist, smile, descname):
    for each in desclist:
        if each["Canonical_Smiles"] == smile:
            return float(each[descname])

def GetLigandDesc(ligandfile, descfile, descname):
    desccontent = []
    with open(descfile, 'rb') as descobj:
        descreader = csv.DictReader(descobj, delimiter="\t")
        for row in descreader:
            desccontent.append(row)
    ligdict = dict()
    with open(ligandfile, 'rb') as ligobj:
        ligreader = csv.DictReader(ligobj, delimiter="\t")
        for row in ligreader:
            ligdict[row['ligandid']] = GetDescValue(desccontent, row['Canonical_Smiles'], descname)
    return ligdict

def GetLigandChiralMain():
    filedir = __FILEDIR__
    ligdict = GetLigandDesc(os.path.join(filedir, "ligand_5_7.txt"), os.path.join(filedir, "descriptor_5_17.txt"), "chiral")
    #jdict = json.load(open("test.json"))
    jdict = json.load(open("all_0.9.json"))
    newjson = NewJSON(jdict)
    newjson.AddField(ligdict, "chiral")
    ligdict = GetLigandTargetName(os.path.join(filedir, "ligand_5_7.txt"), os.path.join(filedir, "proteinseq_5_4.txt"))
    newjson.AddField(ligdict, "target")
    #fileobj  = open("withchiral_small.json", "w")
    fileobj  = open("withchiral.json", "w")
    fileobj.write(json.dumps(jdict, indent=2))

if __name__ == "__main__":
    #GetProteinTargetMain()
    GetLigandChiralMain()
