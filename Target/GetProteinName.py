'''
    Get ligand protein targets
    I know this can be highly redundant, but at this stage I want to impliment this function
    @Date: 5/18/2014
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

def AddProteinWrapper(valname, value):
    def inner(func):
        setattr(func, valname, value)
        return func
    return inner

@AddProteinWrapper("lignames", [])
def AddProteinTarget2JSON(jsondict, ligdict):
    lignames = AddProteinTarget2JSON.lignames
    if not lignames:
        lignames += ligdict.keys()
    if "children" in jsondict:
        for child in jsondict["children"]:
            AddProteinTarget2JSON(child, ligdict)
    else:
        if jsondict["name"] in lignames:
            jsondict["target"] = ligdict[jsondict["name"]]

if __name__ == "__main__":
    filedir = __FILEDIR__
    ligdict = GetLigandTargetName(os.path.join(filedir, "ligand_5_7.txt"), os.path.join(filedir, "proteinseq_5_4.txt"))
    print TopNBlue(ligdict, 10, os.path.join(filedir, "ligand_5_7.txt"))
    #jdict = json.load(open("all_0.9.json"))
    jdict = json.load(open("test.json"))
    AddProteinTarget2JSON(jdict, ligdict)
    fileobj  = open("withtarget_small.json", "w")
    fileobj.write(json.dumps(jdict, indent=2))
