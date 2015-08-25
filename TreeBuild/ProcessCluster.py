# Linux version for processlingandcluster.py
# Author:ajing
# Date:10/12/2012
# Three input for this program, all are directories.
# ChangeLog:

from datetime import date
import os
import sys
import subprocess
from get_protein_name import get_proteinname_fromsmile
################### be careful for different operating system ###########
_work_dir = '/users/ajing'
path = os.path.join(_work_dir,'pylib')
print _work_dir
print path
sys.path.append(os.path.join(_work_dir,'pylib'))
#from hist_twice import hist_once

class processligandcluster:
    def __init__(self, bindingtype = "allosteric", threshold = 100):
        self.threshold = threshold
        self.bindingtype = bindingtype
        self.ligcluster_file = os.path.join(_work_dir, "Cluster/database/ligandcluster_" + str(date.today().month)+ '_' + str(date.today().day)+".txt")
        print self.ligcluster_file
        self.figfiledir = os.path.join(_work_dir, "Cluster/database")
        self.descriptor_dir = os.path.join(_work_dir, "Cluster/database/descriptor_5_17.txt")
        self.outfile = os.path.join(_work_dir, "Cluster/database/" + bindingtype[:4] + "_" + str(threshold) + "/" + "centerdescriptor_" + bindingtype[:4] + '_' + str(date.today().month)+ '_' + str(date.today().day) + ".txt")
        ##############Check whether certain directories exist###############
        if not os.path.exists(self.ligcluster_file):
            sys.exit("ligand cluster file cannot be found")
        if not os.path.exists(self.descriptor_dir):
            sys.exit("descriptor file cannot be found")
        ############## Partial process of files#############
        ligcluster_info = self.get_allinfo_incluster(self.ligcluster_file, 0)
        #print ligcluster_info
        (desc_info, header) = self.get_allinfo_incluster(self.descriptor_dir, 1)
        self.get_center_descriptor(ligcluster_info, desc_info, header, bindingtype, self.outfile)
        if bindingtype == "competitive":
            self.concatenate_two_files(threshold)

    def get_allinfo_incluster(self, infile, returntype):
        # To understand ligand cluster file and descriptor file
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
            #print rownum
        if returntype == 0:
            # Process ligand cluster file
            return all_info
        else:
            # Process descriptor file
            return (first_col, header)

    def get_center_descriptor(self, ligcluster, desc_info, header, bindingtype = None, outputfile = None):
        num_ligs = []
        num_cluster = []
        idx_start = 0
        all_des = dict()
        smile_list = []
        protein_id_old = ligcluster['proteinid'][0]
        debug_file = open("debug.txt","w")
        n = 0
        for idx in range(1,len(ligcluster['ClusterCenter'])):
            protein_id_now = ligcluster['proteinid'][idx]
            if protein_id_old != protein_id_now:
                idx_end = idx
                sub_cluster = ligcluster['ClusterCenter'][idx_start:idx_end]
                sub_smiles = ligcluster['Column2'][idx_start:idx_end]
                sub_ligand = ligcluster['Column1'][idx_start:idx_end]
                idx_start = idx
                ## append important information in num_cluster & num_ligs
                num_cluster.append(int(sub_cluster[-1]))
                for i in range(1,int(sub_cluster[-1])+1):
                    num_ligs.append(sub_cluster.count(str(i)))
                    smile_index = sub_cluster.index(str(i))
                    smile = sub_smiles[smile_index]
                    smile_list.append(smile)
                    sub_ligand[smile_index]
                    desc_info[smile]
                    all_des[sub_ligand[smile_index] + str(n)] = [protein_id_old] + [sub_ligand[smile_index]] + [smile] + desc_info[smile] + [str(sub_cluster.count(str(i)))]
                    n += 1
                protein_id_old = protein_id_now
            debug_file.write(str(idx_start)+'\n')
        # For the last line
        debug_file.close()
        idx_end = idx + 1
        sub_cluster = ligcluster['ClusterCenter'][idx_start:idx_end]
        sub_smiles = ligcluster['Column2'][idx_start:idx_end]
        sub_ligand = ligcluster['Column1'][idx_start:idx_end]
        ## append important information in num_cluster & num_ligs
        num_cluster.append(int(sub_cluster[-1]))
        for i in range(1,int(sub_cluster[-1])+1):
            num_ligs.append(sub_cluster.count(str(i)))
            smile_index = sub_cluster.index(str(i))
            smile = sub_smiles[smile_index]
            smile_list.append(smile)
            all_des[sub_ligand[smile_index] + str(n)] =[protein_id_old] + [sub_ligand[smile_index]] + [smile] + desc_info[smile] + [str(sub_cluster.count(str(i)))]
        print "num_cluster"
        print num_cluster
        ############# 7/16/2012 find biggest cluster with protein name
        self.find_top_protein(zip(num_ligs, smile_list))
        if outputfile is not None:
            outfile_obj = open(outputfile,'w')
            outfile_obj.write('class' +'\tsmile\t' + '\t'.join(header[1:]) + '\t' + "ligcounts" + '\n')
            for each_ligand in all_des.keys():
                line = bindingtype +'\t' + '\t'.join(all_des[each_ligand]) + '\n'
                outfile_obj.write(line)
            outfile_obj.close()
        new_figfiledir = os.path.join(self.figfiledir, "allo_" + str(self.threshold) + "")
        hist_once(num_cluster, os.path.join(new_figfiledir, "num_cluster_" + bindingtype[:4]), "# of clusters per protein family" , "# protein families")
        hist_once(num_ligs, os.path.join(new_figfiledir, "num_ligs_" + bindingtype[:4]), "# of ligands per cluster", "# clusters")

    def find_top_protein(self, lignum_smile_zipped):
        lignum_smile_sorted = sorted(lignum_smile_zipped, key = lambda lignum_smile:lignum_smile[0], reverse = True)
        smile_list = []
        for i in range(20):
            each = lignum_smile_sorted[i]
            smile = each[1]
            smile_list.append(smile)
        proteinname = get_proteinname_fromsmile(smile_list)
        protein_file  = os.path.join(_work_dir, "Cluster/database/" + self.bindingtype[:4] + "_" + str(self.threshold) + "/" + "proteinname_" + self.bindingtype[:4] + "_" + str(date.today().month)+ '_' + str(date.today().day) + ".txt")
        protein_file_obj = open(protein_file,'w')
        protein_file_obj.write('\n'.join(proteinname))

    def concatenate_two_files(self, threshold):
        file1  = os.path.join(_work_dir, "Cluster/database/allo_" + str(threshold) + "/" + "centerdescriptor_allo_" + str(date.today().month)+ '_' + str(date.today().day) + ".txt")
        file2  = os.path.join(_work_dir, "Cluster/database/comp_" + str(threshold) + "/" + "centerdescriptor_comp_" + str(date.today().month)+ '_' + str(date.today().day) + ".txt")
        final_file = os.path.join(_work_dir, "Cluster/database/allo_" + str(threshold) + "/centerdescriptor_" + str(threshold) + "_" + str(date.today().month)+ "_" + str(date.today().day) + ".txt")
        cmd = ["cp", file1, final_file]
        cmd2 = ["tail", "-n", "+2", file2, ">>", final_file]
        command = subprocess.list2cmdline(cmd)
        command2 = subprocess.list2cmdline(cmd2)
        p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE).communicate()
        p2 = subprocess.Popen(command2, shell = True, stdout = subprocess.PIPE).communicate()

if __name__ == "__main__":
    t = processligandcluster()
