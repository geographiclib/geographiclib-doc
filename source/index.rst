.. GeographicLib documentation master file, created by
   sphinx-quickstart on Sun Apr  3 14:53:09 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GeographicLib
=============

**GeographicLib** is collection of libraries descended from an original
C++ library by `Charles Karney <mailto:karney@alum.mit.edu>`_ which
started in 2008.  That `C++ library <C++/doc/index.html>`_ now offers:

* geodesic and rhumb line calculations;

* conversions between geographic, UTM, UPS, MGRS, geocentric, and local
  cartesian coordinates;

* gravity (e.g., EGM2008) and geomagnetic field (e.g., WMM2020)
  calculations.

Portions of the library (particularly the geodesic routines) have also
been implemented in a few :ref:`other languages <languages>`.  Prior to
version 2.0, all the implementations were bundled together in the same
repository.  Now they have been separated into a family of repositories
on `github <https://github.com/geographiclib>`_.

Links
-----

* `C++ library documentation <C++/doc/index.html>`_
* `SourceForge project page <https://sourceforge.net/projects/geographiclib>`_
  for GeographicLib
* `Implementations in other languages <doc/library.html#languages>`_
* `MIT license <LICENSE.txt>`_
* `old home page <oindex.html>`_ for GeographicLib

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   doc/library
   doc/research

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
