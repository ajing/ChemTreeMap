Installation
============

ChemTreeMap contains two parts: frontend (JavaScript) and backend (Python). Users can install ChemTreeMap backend as python package. We support two ways to install ChemTreeMap.

- Docker installation: Run ChemTreeMap in a Docker container isolated from all other programs on your machine.
- Source code installation: Several packages are needed to fully install ChemTreeMap. Please check this wiki page for more details.


Docker installation (fastest install)
-------------------------------------

Docker is a system which can be used to build self contained versions of a Linux operating system running on your machine. When you install and run ChemTreeMap via Docker, it completely isolates the installation from pre-existing packages on your machine.

Please follow the following steps to install ChemTreeMap Docker container. (Superuser privilege may be necessary for the following commands)

See `installing Docker`_ for instructions on installing Docker on your machine.

After Docker is installed, pull ChemTreeMap image as follows.

.. code-block:: bash

  $ docker pull ajing/chemtreemap

Then, launch a Docker container with the binary image as follows.

.. code-block:: bash

  $ docker run -t -i -p 8000:8000 ajing/chemtreemap /bin/bash

The example code is in /examples folder. Please read and run the examples.py file to see how it works with example input files.

.. code-block:: bash

  $ cd examples
  $ python examples.py

Then, open the following URL in your Chrome/Firefox.

http://localhost:8000/dist/#/aff

'aff' may be changed to your input filename without file extension.

If you have your own molecule files, please prepare your input file following the same tab delimited format as example files (e.g. aff.txt, factorxa.txt, etc.). The following command can be used to import files into Docker container.

.. code-block:: bash

  $ docker cp foo.txt ajing/chemtreemap:/examples/foo.txt

Please change examples.py (variables input_file, out_file) accordingly and run the examples.py.


.. _installing Docker: https://docs.docker.com/engine/installation/