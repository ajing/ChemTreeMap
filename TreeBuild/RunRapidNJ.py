'''
    Run RapidNJ on phylip dist file
'''

import subprocess

RAPIDNJ_COMMAND = "./lib/rapidnj-linux-64"

def RunRapidNJ(distant_file):
    proc   = subprocess.Popen([RAPIDNJ_COMMAND, distant_file, "-i", "pd"], stdout=subprocess.PIPE)
    newick = proc.stdout.read()
    return newick

if __name__ == "__main__":
    #print RAPIDNJ_COMMAND % ("filename")
    RunRapidNJ("./Data/2015-08-24-13h-49m.dist")
