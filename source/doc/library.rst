.. _library:

Library
=======

GeographicLib started off as a C++ library.  However the component of
the library solving the geodesic problem was sufficiently novel and
sufficiently important that there was a need to port it to several other
languages.

Infrastructure
--------------

The repositories for all implementations are on `GitHub
<https://github.com/orgs/geographiclib/repositories>`_.  The project is
hosted on `SourceForge
<https://sourceforge.net/projects/geographiclib>`_ and it hosts
GeographicLib's `web pages <../index.html>`_ and
the `file system
<https://sourceforge.net/projects/geographiclib/files>`_ for package
distribution.

.. _languages:

Supported languages
-------------------

The C++ version of GeographicLib offers the broadest range of features.
However, the geodesic routines are available several other languages.
Here is a list of the supported languages.

================= ==============  ================== =========================
language          information     repository\ [#a]_  download\ [#b]_
================= ==============  ================== =========================
C++               `C++ doc`_      `C++ repo`_\ [#c]_ `C++ packages`_\ [#d]_
C                 `C doc`_        `C repo`_          `C packages`_\ [#e]_
Fortran           `Fortran doc`_  `Fortran repo`_    `Fortran packages`_
Python            `Python doc`_   `Python repo`_     `Python packages`_\ [#f]_
Octave\ [#g]_     `Octave doc`_   `Octave repo`_     `Octave packages`_\ [#h]_
Java\ [#i]_       `Java doc`_     `Java repo`_       `Java packages`_\ [#j]_
JavaScript\ [#k]_ `JS doc`_       `JS repo`_         `JS packages`_\ [#l]_
================= ==============  ================== =========================

.. [#a] Active branches are ``main`` and ``devel``.  For the C++ repository,
        releases are available on the ``release`` branch.
.. [#b] The download links offered here are the ones updated as soon as
        there's a release.  The packages are also available, as noted,
        from other downstream providers.  Prior to version 2.0, all the
        languages were bundled into a single source package available
        `here
        <https://sourceforge.net/projects/geographiclib/files/distrib>`_.
.. [#c] This repository also hold the implementations in other languages
        for versions before 2.0.
.. [#d] Also available in several Linux distributions (Debian, Ubuntu,
        Fedora, etc.) and as packages on `conda-forge (C++)
        <https://anaconda.org/conda-forge/geographiclib-cpp>`_,
        `vcpkg <https://vcpkg.info/port/geographiclib>`_, and
        `brew <https://formulae.brew.sh/formula/geographiclib>`_.
.. [#e] The C library is incorporated into `PROJ
        <https://proj.org/geodesic.html>`_ and the R package `geosphere
        <https://cran.r-project.org/package=geosphere>`_.
.. [#f] Also available as a `PyPI package
        <https://pypi.python.org/pypi/geographiclib>`_ and in
        `conda-forge (Python) <https://anaconda.org/conda-forge/geographiclib>`_
.. [#g] This is compatible with both Octave and MATLAB.  This package
        also includes several other components of the C++ library and a
        treatment of great ellipses (not in C++ library).
.. [#h] Also available as a `MATLAB Central package
        <https://www.mathworks.com/matlabcentral/fileexchange/50605>`_.
.. [#i] Also includes the gnomonic projections.
.. [#j] Also available as a `Maven Central package
        <https://search.maven.org/artifact/net.sf.geographiclib/GeographicLib-Java>`_
.. [#k] Consists of 2 packages ``geographiclib-geodesic`` and
        ``geographiclib-dms`` for handling geodesic calculations and DMS
        strings, respectively.
.. [#l] Also available as npm packages `geographiclib-geodesic
        <https://www.npmjs.com/package/geographiclib-geodesic>`_ and
        `geographiclib-dms
        <https://www.npmjs.com/package/geographiclib-dms>`_.

.. _C++ doc:  ../C++/doc/index.html
.. _C++ repo: https://github.com/geographiclib/geographiclib
.. _C++ packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-C++

.. _C doc:  ../C/doc/index.html
.. _C repo: https://github.com/geographiclib/geographiclib-c
.. _C packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-C

.. _Fortran doc:  ../Fortran/doc/index.html
.. _Fortran repo: https://github.com/geographiclib/geographiclib-fortran
.. _Fortran packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-Fortran

.. _Python doc:  ../Python/doc/index.html
.. _Python repo: https://github.com/geographiclib/geographiclib-python
.. _Python packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-Python

.. _Octave doc:
   https://github.com/geographiclib/geographiclib-octave#readme
.. _Octave repo: https://github.com/geographiclib/geographiclib-octave
.. _Octave packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-Octave

.. _Java doc:  ../Java/doc/index.html
.. _Java repo: https://github.com/geographiclib/geographiclib-java
.. _Java packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-Java

.. _JS doc:  ../JavaScript/doc/index.html
.. _JS repo: https://github.com/geographiclib/geographiclib-js
.. _JS packages:
   https://sourceforge.net/projects/geographiclib/files/distrib-JavaScript

Online utilities
----------------

The C++ library provides some online tools using GeographicLib

* `geographic coordinate conversions <../cgi-bin/GeoConvert>`_
  between latitude/longitude, UTM or UPS, and MGRS;

* `direct and inverse geodesic calculations <../cgi-bin/GeodSolve>`_;

* `calculate the perimeter and area of geodesic polygons <../cgi-bin/Planimeter>`_;

* `various geodesic calculations using JavaScript <../scripts/geod-calc.html>`_;

* `a tool for displaying geodesics on Google Maps
  <../scripts/geod-google.html>`_;

* `rhumb line calculator <../cgi-bin/RhumbSolve>`_;

* `evaluate the geoid height <../cgi-bin/GeoidEval>`_ for EGM84, EGM96, and
  EGM2008.


Datasets
--------

Some classes in the C++ library and the geoid functions in the Octave
package require datasets to be available.  Here are the links for the
instructions for how to download and install these.

=============== ================= ========================
Dataset
=============== ================= ========================
geoids          `geoid data`_     `geoid instructions`_
gravity models  `gravity data`_   `gravity instructions`_
magnetic models `magnetic data`_  `magnetic instructions`_
=============== ================= ========================

.. _geoid data:
   https://sourceforge.net/projects/geographiclib/files/geoids-distrib
.. _geoid instructions:
   ../C++/doc/geoid.html#geoidinst
.. _gravity data:
   https://sourceforge.net/projects/geographiclib/files/gravity-distrib
.. _gravity instructions:
   ../C++/doc/gravity.html#gravityinst
.. _magnetic data:
   https://sourceforge.net/projects/geographiclib/files/magnetic-distrib
.. _magnetic instructions:
   ../C++/doc/magnetic.html#magneticinst

Test data
---------

High quality test data is used to test GeographicLib for

* the `geodesic problem <https://doi.org/10.5281/zenodo.32156>`_
* the `transverse Mercator projection <https://doi.org/10.5281/zenodo.32470>`_
* `geoid heights <../C++/doc/geoid.html#testgeoid>`_


Other implementations
---------------------

Others have written independent implementations of the geodesic
algorithms or have built other tools using GeographicLib:

* *Mathematica*: Kei Misawa, `mathematica-geodesic
  <https://github.com/330k/mathematica-geodesic>`_.

* *R*: R. J. Hijmans, `geosphere
  <https://cran.r-project.org/package=geosphere>`_.

* *Rust*: Michael Kirk, `Rust implementation of GeographicLib
  <https://github.com/georust/geographiclib-rs>`_.

* *Go*: Patrick Yukman, `A Go port of GeographicLib
  <https://github.com/pymaxion/geographiclib-go>`_.

* *Python via Cython*: Sergey Serebryakov, `Cython extension module for
  GeographicLib
  <https://github.com/megaserg/geographiclib-cython-bindings>`_.

* *JavaScript via emscripten*: William Wall, `OpenSphere ASM
  <https://github.com/ngageoint/opensphere-asm>`_.

* *JavaScript*: Jaco, `a graphical tool for geodesic calculations
  <http://geo.javawa.nl/coordcalc/index_en.html>`_.
