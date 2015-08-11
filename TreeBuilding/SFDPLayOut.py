'''
    Add sfdp layout to dotfile
'''

import subprocess
import datetime
import os

def SFDPonDot(dotfile, size):
    fmt='./Data/%Y-%m-%d-%Hh-%Mm_dot_sfdp.gv'
    newfilename = datetime.datetime.now().strftime(fmt)
    if os.path.isfile(newfilename):
        os.remove(newfilename)
    command = "sfdp -Gsmoothing=triangle -Gsize={size} {infile} > {outfile}".format(size=size, infile=dotfile, outfile=newfilename)
    subprocess.Popen( command, shell = True, stdout = subprocess.PIPE ).communicate()
    RemoveBackSlash(newfilename)
    return newfilename

def RemoveBackSlash(dotfile):
    # remove backslash and replace all " quote sign
    f = open(dotfile, 'r+')
    content = f.readlines()
    newcontent = []
    for line in content:
        line = line.replace("\"", "")
        if line.endswith("\\\n"):
            newcontent.append(line[:-2])
        else:
            newcontent.append(line)
    f.seek(0)
    f.write("".join(newcontent))
    f.truncate()
    f.close()

def test():
    dot = "./Data/test.gv"
    SFDPonDot(dot, 10)
    #dot = "./Data/sfdp.gv"
    #RemoveBackSlash(dot)

if __name__ == "__main__":
    test()
