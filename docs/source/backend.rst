Backend
=======

The backend is written in `Python`, and is where data is processed into the JSON
representation. The API is detailed in :doc:`modules`.

Installation
------------

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

The server will server the output in the end of ``examples.py``

Data structure
^^^^^^^^^^^^^^

Compounds
~~~~~~~~~

The example script takes compounds as a smiles file. The smiles should be
provided with column name Canonical_Smiles, and the identifier column is set
as a parameter of TreeBuild object.  The individual activities and properties
are supplied as extra columns.

A sample of a compound dataset would be:

.. code-block:: bash

    $ head aff.txt
    ligandid        Canonical_Smiles        pIC50
    Chk1N144        COc1cc2Cc3c(n[nH]c3c4ccc(nc4)c5ccc(O)cc5)c2cc1OC        9.52432881
    Chk1N145        COc1cc2Cc3c(n[nH]c3c4ccc(nc4)C#N)c2cc1OC        9.09151498
    Chk1N40 OC1CCCN(C1)c2ccc(Br)cc2NC(=O)Nc3cnc(cn3)C#N     8.86646109
    Chk1N115        COc1c(OC)cc2c(Cc3c2n[nH]c3c4ccc(c5ccc(O)cc5)cc4)c1      8.79588002
    ... ... ...

Showing 4 compounds, with pIC50 as activity.


Running the Script and Viewing the Tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The script will use the default fingerprint (ECFP and AtomPair) and
caculated properties (pIC50, SlogP and ligand efficiency) to generate
the ChemTreeMap data, and produce the JSON as an output file, specified
in out_file parameter of TreeBuild object.


.. code-block:: bash

    $ python examples.py
    ... ... ...

The generated data will be copied to the data directory of the :doc:`frontend` and
view the data in browser.

Then open the browser Chrome/Firefox with the `link
<http://localhost:8000/dist/index.html#/aff>`_.

`aff` need to be changed to the filename of your input file without file extention.

Please consider looking at ``treebuild/examples/examples.py`` where there are five examples
using the :doc:`declarative API <modules>` for generating the datasets.