Backend
=======

The backend is written in `Python`, and is where data is processed into the JSON
representation. The API is detailed in :doc:`modules`.

Installation
------------
First, install rdkit (http://www.rdkit.org/), ete (http://etetoolkit.org/), and graphviz (http://www.graphviz.org).

Switch to the backend directory and run the `setup.py` script.

.. code-block:: bash

    $ python setup.py install

Quickstart
----------

A prebuilt script, ``examples.py`` for analysing composite biochemical data
in a standard format is provided in the ``examples`` folder.  This provides
minimum functionality - the software is designed for extensibility,
so consider writing your own scripts like those in the example datasets
included in ``backend/treebuild/examples/``.

The script can be run using:

.. code-block:: bash

    $ python examples.py

The server will serve the output upon running ``examples.py``

Data structure
^^^^^^^^^^^^^^

Compounds
~~~~~~~~~

The example script takes compounds as a smiles file. The smiles should be
provided with column name Canonical_Smiles, and the identifier column is set
as a parameter of TreeBuild object.  The individual activities and properties
are supplied as extra columns.

A sample of a compound dataset would be delimited by tabs:

.. code-block:: bash

    Canonical_Smiles        ligandid        pIC50
    COc1cc2Cc3c(n[nH]c3c4ccc(nc4)c5ccc(O)cc5)c2cc1OC        Chk1N144        9.52432881
    COc1cc2Cc3c(n[nH]c3c4ccc(nc4)C#N)c2cc1OC        Chk1N145        9.09151498
    OC1CCCN(C1)c2ccc(Br)cc2NC(=O)Nc3cnc(cn3)C#N     Chk1N40 8.86646109
    COc1c(OC)cc2c(Cc3c2n[nH]c3c4ccc(c5ccc(O)cc5)cc4)c1      Chk1N115        8.79588002
    ... ... ...

Showing 4 compounds, with pIC50 as activity.


Running the Script and Viewing the Tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The script will use the default fingerprint (ECFP and AtomPair) and
caculated properties (pIC50, SlogP and ligand efficiency) to generate
the ChemTreeMap data, and produce JSON output file, specified
in the out_file parameter of TreeBuild object.


.. code-block:: bash

    $ python examples.py
    ... ... ...

The generated data will be copied to the data directory of the :doc:`frontend` and
viewed in a browser.

Then open the browser Chrome/Firefox with the `link
<http://localhost:8000/dist/index.html#/aff>`_ (http://localhost:8000/dist/index.html#/aff).

`aff` needs to be changed to the filename of your input file without file extention.

Please consider looking at ``treebuild/examples/examples.py`` where there are five examples
using the :doc:`declarative API <modules>` for generating the datasets.