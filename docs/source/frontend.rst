Frontend
========

The frontend is written in JavaScript, and runs in modern (HTML5 enabled) web browsers.
It is an `angularjs`_ app, and was scaffolded using `yeoman`_, using the `angular`_ generator.
The visualization uses the `D3.js`_ library to render the images.

Dependencies
------------

The project requires the `node`_ runtime to build the app, `Grunt`_ to run tasks,
and the package managers `npm`_ and `bower`_ to install development and production
dependencies respectively.

Building
--------

To build and serve the project (having installed `node`_ and `npm`_), change into the frontend directory and run:

.. code-block:: bash

  $ npm install -g grunt bower      // install grunt and bower
  $ npm install                     // install node modules
  $ bower install                   // install bower modules
  $ grunt build                     // build the server
  $ python -m SimpleHTTPServer 8000 // run the server

Your default browser should open by default to ``localhost:8000``, where the app is being
served.

.. _angularjs: http://angularjs.org
.. _yeoman: http://yeoman.io
.. _angular: https://github.com/yeoman/generator-angular
.. _D3.js: https://d3js.org
.. _node: https://nodejs.org
.. _bower: https://bower.io
.. _npm: http://npmjs.com
.. _grunt: http://gruntjs.com
