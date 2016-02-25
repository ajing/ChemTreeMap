Documentation
=============

Current readthedoc doesn't include detailed API for modules. More informative
document can be built locally.

The documentation you are reading is included in the docs subdirectory of the
project repository.  It is build using `sphinx`_ and hosted using `readthedocs`_.

Dependencies
------------

To get it running locally, install sphinx and the readthedocs theme:

.. code-block:: bash

  $ pip install sphinx readthedocs

Building
--------

Documentation can then be built using the Makefile (or the make.bat for Windows):

.. code-block:: sh

  $ make html

Live editing
------------

For quick editing, we use sphinx-autodoc.

.. code-block:: bash

  $ pip install sphinx-autodoc
  $ make livehtml


.. _sphinx: http://sphinx-doc.org
.. _readthedocs: http://readthedocs.org
