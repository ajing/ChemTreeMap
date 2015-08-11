# RDKit
sudo apt-get install -y python-rdkit librdkit1 rdkit-data
# NAMS
#easy_install NAMS
# ete
#sudo apt-get install -y python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml
#sudo easy_install -U ete2
# networkx
sudo apt-get install -y python-pip
sudo pip install networkx
# NAMS
#sudo easy_install NAMS
#sudo apt-get install -y python-openbabel
#sudo apt-get install -y python-munkres
# numpy
sudo apt-get install -y python-numpy
# to draw molecule structure jpg file
# download oasa from http://bkchem.zirael.org/
sudo apt-get install -y bkchem


# Please also install node.js server
# follow the link: https://www.digitalocean.com/community/tutorials/how-to-set-up-a-node-js-application-for-production-on-ubuntu-14-04

# for node.js
curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.25.4/install.sh | bash

# Load nvm and install latest production node
source $HOME/.nvm/nvm.sh
nvm install v0.12.7
nvm use v0.12.7

# Install jshint to allow checking of JS code within emacs
# http://jshint.com/
npm install -g jshint

# Install rlwrap to provide libreadline features with node
# See: http://nodejs.org/api/repl.html#repl_repl
sudo apt-get install -y rlwrap

# install or update the package in package.json
npm update
