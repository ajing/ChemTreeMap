# RDKit
sudo apt-get install -y python-rdkit librdkit1 rdkit-data
# NAMS
#easy_install NAMS
# ete
sudo apt-get install -y python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml
sudo easy_install -U ete2
# networkx
sudo apt-get install -y python-pip
sudo pip install networkx
# NAMS
sudo easy_install NAMS
sudo apt-get install -y python-openbabel
sudo apt-get install -y python-munkres
# numpy
sudo apt-get install -y python-numpy
# to draw molecule structure jpg file
# download oasa from http://bkchem.zirael.org/
sudo apt-get install -y bkchem

# graphviz
sudo apt-add-repository ppa:dperry/ppa-graphviz-test
sudo apt-get update
sudo apt-get autoremove graphviz
sudo apt-get install -y graphviz
