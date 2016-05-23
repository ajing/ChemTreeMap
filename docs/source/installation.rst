Installation
============

ChemTreeMap contains two parts: frontend (JavaScript) and backend (Python). Two methods of installation are supported.

- Docker installation: Run ChemTreeMap in a Docker container with a reproducible environment.
- `Source code installation`_: The current source code installation instruction is only tested for Ubuntu (Linux) system.
the dependency rapidnj (http://birc.au.dk/software/rapidnj/) and graphviz (http://www.graphviz.org/)
  do not have command line support for Windows. Please check this
  `wiki page`_ (https://github.com/ajing/ChemTreeMap/wiki/Souce-code-installation) for more details.
  The ChemTreeMap backend can be installed as a python package.


Docker installation (fastest install)
-------------------------------------

Docker provides a reproducible environment for scientific software. All dependencies are pre-installed in the Docker container.

Docker is a system which can be used to build self contained versions of a Linux operating system running on your machine. When you install and run ChemTreeMap via Docker, it completely isolates the installation from pre-existing packages on your machine.

Please follow the following steps to install the ChemTreeMap Docker container. (Superuser privilege may be necessary for the following commands)

See `installing Docker`_ (https://docs.docker.com/engine/installation/) for instructions on installing Docker on your machine.
Note: VirtualBox (https://www.virtualbox.org) must be installed to run Docker.

  Installing Docker on Windows:
https://docs.docker.com/engine/installation/windows/

  Installing Docker on Mac:
https://docs.docker.com/engine/installation/mac/

  Installing Docker on Ubuntu:
https://docs.docker.com/engine/installation/linux/ubuntulinux/



After Docker is installed, pull ChemTreeMap image using the following commands.

.. code-block:: bash

  $ docker pull ajing/chemtreemap

Then, launch a Docker container with the binary image as follows.

.. code-block:: bash

  $ docker run -t -i -p 8000:8000 ajing/chemtreemap /bin/bash

The example code is in ~/examples folder. Please read and run the examples.py file to see how it works with example input files.

.. code-block:: bash

  $ cd examples
  $ python examples.py

Then, open the following URL in Chrome/Firefox web browser. The following URL is for Linux Systems (Mac, Ubuntu, Fedora ...).

http://localhost:8000/dist/#/aff

'aff' may be changed to your input filename without file extension.

In using a Windows system, the docker daemon is running inside a Linux virtual machine. So, you need to use the IP address of your
Linux virtual machine. It is shown in your Docker QuickStart Terminal (Figure 1.2) in the red box:

.. figure:: docker.png

   The IP address (red box) on Docker QuickStart Terminal.

With the IP address, open the visualization with the following URL in your Chrome/Firefox browser. In this example, use

http://192.168.99.100:8000/dist/#/aff

If you have your own molecule files, please prepare your input file following the same tab delimited format as example files (e.g. aff.txt, factorxa.txt, etc.). The following command can be used to import files into Docker container.

.. code-block:: bash

  $ docker cp foo.txt ajing/chemtreemap:/examples/foo.txt

Change examples.py (variables input_file, out_file) accordingly and run the examples.py. (Note: you may also need to change the column header to match.)

Now, you are done with ChemTreeMap installation.


Source Code Installation
-------------------------------------

The following sections (1.3 and 1.4) is for installing from source code (Frontend and Backend). The procedure has been tested on Ubuntu 14.04.

First, pull the source code from the GitHub repository.

.. code-block:: bash

  $ git pull https://github.com/ajing/ChemTreeMap


.. _installing Docker: https://docs.docker.com/engine/installation/
.. _Source code installation: https://github.com/ajing/ChemTreeMap/wiki/Souce-code-installation
.. _wiki page: https://github.com/ajing/ChemTreeMap/wiki/Souce-code-installation