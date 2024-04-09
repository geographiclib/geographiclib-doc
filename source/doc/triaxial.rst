.. _triaxial:

Triaxial ellipsoids
===================

Jump to:

* :ref:`tricoords`
* :ref:`tricoordconv`
   + :ref:`tricoordconvh`
   + :ref:`tricoordconvu`
* :ref:`trigeod`
   + :ref:`trigeoddirect`
   + :ref:`trigeodinverse`
   + :ref:`trigeodjac`
* :ref:`triutils`
* :ref:`trirefs`

.. _triintro:

Introduction
------------

This page provided some background information for the ``triaxial`` class
in the Octave/MATLAB package ``geographiclib``.  This component of the
package is still in "prerelease", however, you can download the package
from `SourceForge
<https://sourceforge.net/projects/geographiclib/files/distrib-Octave>`_.
Get either ``geographic_toolbox-2.2.zip`` or
``geographiclib-octave-2.2.tar.gz``.  (These two files contain the
same code; the directory layout for the latter file is tailored for
Octave packages.)

Once the package is in your path, type ``triaxial.doc`` to get an overview
of the class.  You can then use, e.g., ``help triaxial.distance``, to get
the documentation on specific routines.  Use ``triaxial.demo`` to list
various demonstations of the geodesic capabilities.

.. _tricoords:

Coordinate systems
------------------

.. math::
   \DeclareMathOperator*{\hypot}{hypot}

A triaxial ellipsoid is the surface defined by

.. math::
  S(\mathbf r) = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1 = 0,

where :math:`a`, :math:`b`, and :math:`c` are the semiaxes of the
ellipsoid.  We will take :math:`a \ge b \ge c > 0`.  Note that the
ellipsoid is centered at the origin and that its principal axes are
aligned with the coordinate axes.

A direction on the surface can be specified by unit cartesian vector
tangential the surface

.. math
   \mathbf v &= (v_x, v_y, v_z), \\
   \lvert v \rvert &= 1, \\
   \mathbf v \cdot \mathbf U &= 0,

where :math:`\mathbf U(\mathbf r) = \nabla\mathbf S(\mathbf r)` is the
normal to the surface.

A point on the surface is specified by a latitude and longitude.  The
*geodetic* latitude and longitude :math:`(\phi, \lambda)` are defined by

.. math::
 \hat{\mathbf U} = (\cos\phi \cos\lambda, \cos\phi \sin\lambda, \sin\phi),

where the "hat" symbol denote a unit vector.  The *parametric* latitude
and longitude :math:`(\phi', \lambda')` are defined by

.. math::
 x &= a \cos\phi' \cos\lambda', \\
 y &= b \cos\phi' \sin\lambda', \\
 z &= c \sin\phi'.

The *geocentric* latitude
and longitude :math:`(\phi'', \lambda'')` are defined by

.. math::
 \hat{\mathbf r} = \begin{bmatrix}
 \cos\phi'' \cos\lambda'' \\ \cos\phi'' \sin\lambda'' \\ \sin\phi''
 \end{bmatrix}.

As with ellipsoids of revolution, the geodetic, parametric, and
geocentric coordinates are closely related to one another.

Finally *ellipsoid* latitude and longitude :math:`(\beta, \omega)` are
defined by

.. math::
  x &= a \cos\omega
      \frac{\sqrt{a^2 - b^2\sin^2\beta - c^2\cos^2\beta}}
           {\sqrt{a^2 - c^2}}, \\
  y &= b \cos\beta \sin\omega, \\
  z &= c \sin\beta
      \frac{\sqrt{a^2\sin^2\omega + b^2\cos^2\omega - c^2}}
           {\sqrt{a^2 - c^2}}.

A notable feature of ellipsoidal coordinates is that they are orthogonal
(unlike either geodetic or parametric coordinates).  The ellipsoidal
azimuth is then well defined.  In the limit of an oblate ellipsoid of
revolution, :math:`(\beta, \omega)` play the roles of the parametric
latitude and the longitude.  For a prolate ellipsoid, these two roles
are switched.

There are two useful representation of arbitrary points in
three-dimensional space.  There first represents positions by

.. math::
   \mathbf r = \mathbf r_0 + h \hat{\mathbf U}(\mathbf r_0),

where :math:`\mathbf r_0` is the closest point on the ellipsoid and
:math:`h` is the height above the ellipsoid.  Since geodetic coordinates
specify the direction of :math:`\mathbf U`, we can also represent this
point be appending the height to the geodetic coordinates to give
:math:`(\phi, \lambda, h)`.

The second uses ellipsoidal coordinates.  For an arbitrary point
:math:`\mathbf r`, we seek the value of :math:`u` such that

.. math::
   \frac{x^2}{u^2 + l_a^2} +
   \frac{y^2}{u^2 + l_b^2} + \frac{z^2}{u^2} = 1,

where

.. math::
   l_a = \sqrt{a^2 - c^2}, \quad
   l_b = \sqrt{b^2 - c^2}

are linear eccentricities.  This is an ellipsoid which is confocal to
the original one (with semiaxes :math:`a, b, c`) and whose minor
semiaxis is :math:`u`.

.. _tricoordconv:

Coordinate system conversions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conversions between thes coordinate on the surface of the ellipsoid are
algebraic exercises.  For example, the conversion from cartesian to
geodetic coordinates proceeds as follows

.. math::
   \xi &= x/a^2, \quad \eta = y/b^2, \quad \zeta = z/c^2, \\
   \phi &= \tan^{-1} \frac\zeta{\hypot(\xi, \eta)}, \\
   \lambda &= \tan^{-1} \frac\eta\xi,

where the quadrant of the angles should be determined by the signs of
the numerators and denominators separately, using, for example, the
library function atan2.

.. _tricoordconvh:

Solving for :math:`h`
^^^^^^^^^^^^^^^^^^^^^

Following [Ligas12]_, we have

.. math::
   \mathbf r_0 &= (x_0, y_0, z_0) = \biggl(
   \frac{a^2x}{p + l_a^2},
   \frac{b^2y}{p + l_b^2},
   \frac{c^2z}{p} \biggr), \\
   h &= \hat{\mathbf U}(\mathbf r_0) \cdot (\mathbf r - \mathbf r_0),

where :math:`p` is the largest root of

.. math::
   f(p) = \biggl(\frac{ax}{p + l_a^2}\biggr)^2 +
   \biggl(\frac{by}{p + l_b^2}\biggr)^2 +
   \biggl(\frac{cz}{p}\biggr)^2 - 1 = 0.

[Ligas12]_ uses Newton's method to find this root; however, with his
choice of starting guess, this sometimes fails to converge.
[Panou+Korakitis22]_ cure this defect by using the bisection method to
find the root.  This is guaranteed to converge but at the high
computation cost of requiring many iterations.

It turns out we can easily fix the problems with Newton's method.  First
of all, we note that :math:`f(p)` has positive double poles at :math:`p
= 0`, :math:`-l_b^2`, and :math:`-l_a^2` and that
:math:`f(p) \rightarrow -1` for :math:`p \rightarrow \pm\infty`.  (For
now, we assume that :math:`x, y, z` are all nonzero.).  Therefore
:math:`f(p)` has a single root for :math:`p \in (0, \infty)`.  In this
region :math:`f'(p) < 0` and :math:`f''(p) > 0`, and, as a consequence,
picking a starting guess for Newton's method between :math:`p = 0` and
the actual root is guaranteed to converge.

To obtain a reasonably tight bound on the root, we use

.. math::
   f(p) &\le \biggl(\frac{cz}{p}\biggr)^2 - 1, \\
   f(p) &\le \biggl(\frac{\hypot(by, cz)}{p + l_b^2}\biggr)^2 - 1, \\
   f(p) &\le \biggl(\frac{\hypot(ax, by, cz)}{p + l_a^2}\biggr)^2 - 1, \\
   f(p) &\ge \biggl(\frac{\hypot(ax, by, cz)}{p}\biggr)^2 - 1.

Because :math:`f'(p) < 0` for :math:`p > 0`, this leads to bounds on
the positive root, :math:`p_{\mathrm{min}} \le p \le p_{\mathrm{max}}`,
where

.. math::

   p_{\mathrm{min}} &= \max(c \lvert z\rvert,
   \hypot(by, cz) - l_b^2,
   \hypot(ax, by, cz) - l_a^2), \\
   p_{\mathrm{max}} &= \hypot(ax, by, cz).

[Panou+Korakitis22]_ substitute :math:`p_{\mathrm{min}} = c \lvert
z\rvert`; they would get better performance using the tighter bound
given here.  [Ligas12]_ uses :math:`p_0 = c\hypot(x, y, z)` for his
initial guess; because :math:`f(p_0)` can then be negative, Newton's
method may fail to converge.

In implementing Newton's method, we neglect any term in the definition
of :math:`f(p)` if its numerator vanishes (even though the denominator
might also vanish).

Provided that :math:`f(p_{\mathrm{min}}) > 0`, we can now start Newton's
method with :math:`p_0 = p_{\mathrm{min}}` and this converges to the
root from below.  If :math:`f(p_{\mathrm{min}}) \le 0` (which can only
happen if :math:`z = 0`), the required solution is :math:`p = 0`.  In
this case, the expression for :math:`\mathbf r_0` is indeterminate, and
we proceed as follows:

* If :math:`x_0` is indeterminate, substitute :math:`x_0 = 0` (this
  can only happen with :math:`x = 0` on a sphere).
* If :math:`y_0` is indeterminate, substitute :math:`y_0 = 0` (this
  can only happen with :math:`y = 0` on an oblate spheroid).
* If :math:`z_0` is indeterminate, substitute :math:`z_0 = c \sqrt{1 -
  x^2/a^2 - y^2/b^2}`.

A few other points to note:

* This prescription obviates the need to enumerate and solve various
  subcases as [Panou+Korakitis22]_ do.
* Newton's method requires about 8 times fewer iterations compared with
  the bisection method.
* The independent variable :math:`f(p)` is shifted with respect to the
  one used by [Ligas12]_ and [Panou+Korakitis22]_.  This gives higher precision
  close to the singularity at :math:`p = 0`.
* We acculate the terms in :math:`f(p)` in a two-word accumulator to
  improved the accuacy near its root.
* To avoid potentially singular behavior, we initially "flush" tiny
  values of the components of :math:`\mathbf r` to zero.
* In the case where :math:`z_0` is indeterminate, the sign of :math:`z`
  should be used to determine the sign of the square root above.
* If need be, this method is easily generalized to ellipsoids in
  higher dimensions.

.. _tricoordconvu:

Solving for :math:`u`
^^^^^^^^^^^^^^^^^^^^^

Writing :math:`u^2 = q`, we need to find the roots of

.. math::
   g(q) = \frac{x^2}{q + l_a^2} + \frac{y^2}{q + l_b^2} + \frac{z^2}{q} - 1
   = 0.

The structure of :math:`g(q)` is very similiar to :math:`f(p)`.  Since
:math:`g(q)` has 3 simple poles with positive coeffients, there are
three real roots and, because the rightmost pole is at :math:`q = 0`,
just one of them is positive.  As before, bounds can be put on this root
:math:`q_{\mathrm{min}} \le q \le q_{\mathrm{max}}`,
where

.. math::
   q_{\mathrm{min}} &= \max(z^2,
   y^2 + z^2 - l_b^2,
   x^2 + y^2 + z^2 - l_a^2), \\
   q_{\mathrm{max}} &= x^2 + y^2 + z^2.

As before, in implementing Newton's method, we neglect any term in the
definition of :math:`g(q)` if its numerator vanishes (even though the
denominator might also vanish).

Provided that :math:`g(q_{\mathrm{min}}) > 0`, we can now start Newton's
method with :math:`q_0 = q_{\mathrm{min}}` and this converges to the
root from below.  If :math:`g(q_{\mathrm{min}}) \le 0` (which can only
happen if :math:`z = 0`), the required solution is :math:`q = u = 0`.

Of course, we can expand out :math:`g(q)` to obtain a cubic polynomial
in :math:`q` which cab be solved analytically, see [DLMF]_,
Secs. 1.11(iii) and 4.43.  This method is advocated by
[Panou+Korakitis21]_.  However, this solution suffers from roundoff
error when the coefficient of :math:`q` is positive; in this case, the
polynomial in :math:`1/q` should be solved instead.  Even so, the
solution may be subject to unacceptable roundoff error; it may be
refined by using as the startig point, :math:`q_0`, for Newton's method.
If :math:`g(q_0) < 0`, then :math:`q_1` should be replaced by
:math:`\max(q_1, q_{\mathrm{min}})`.  Typically only 3--4 iterations
are needed to refine the solution.

Note: tighter bounds can be placed on :math:`q` using

.. math::
   g(q) &\le \frac{y^2}e{q + l_b^2} + \frac{z^2}{q} - 1 \\
   g(q) &\le \frac{x^2+y^2}{q + l_a^2} + \frac{z^2}{q} - 1 \\
   g(q) &\le \frac{x^2}{q + l_a^2} + \frac{y^2+z^2}{q + l_b^2} - 1 \\
   g(q) &\ge \frac{x^2+y^2}{q + l_b^2} + \frac{z^2}{q} - 1 \\
   g(q) &\ge \frac{x^2}{q + l_a^2} + \frac{y^2+z^2}{q} - 1

and solving the resulting quadratic equations.  This yields only a
marginal improvement given that we're starting with the root of the
cubic.

.. _trigeod:

Geodesics
---------

The problem of geodesics on a triaxial ellipsoid was solved by
[Jacobi39]_ who reduced the problem to quadrature.  Even without
evaluating the integrals, this solution allowed the various properties
of geodesics to be found.

Explicit evaluation of Jacobi's integrals was carried out by hand by
[Cayley72]_ and, more recently, by [Baillard13]_.  Accurate evalution of
the integrals involves changing the variable of integration using
elliptic integrals and elliptic functions.  Unfortunately, Octave/MATLAB
has poor support for these special functions, so for this implementation
of the geodesic routines, I instead integrate the geodesic equations in
cartesian coordinates, following [Panou+Korakitis19]_.

.. _trigeoddirect:

The direct problem
^^^^^^^^^^^^^^^^^^

The equation for geodesics on a surface is the same as for the motion of
a particle constrained to move on the surface but subject to no other
forces.  The centrifugal acceleration of the particle is
:math:`-(v^2/R)\hat{\mathbf U}` where :math:`R` is the radius of
curvature in the direction of :math:`\mathbf v`.  We will take the speed
to be unity (and, of course, the speed is a constant in this problem);
thus time can be replaced by :math:`s`, the displacement along the
geodesic, as the independent variable.  The differential equation for
the geodesic is

.. math::
   d\mathbf r / ds &= \mathbf v, \\
   d\mathbf v / ds &= \mathbf A, \\
   d^2 m / ds^2 &= -K m,

where

.. math::
   \mathbf A &= \frac{\mathbf U}{U^2}
   \biggl( \frac{v_x^2}{a^2} + \frac{v_y^2}{b^2} + \frac{v_z^2}{c^2} \biggr),\\
   K &= \frac1{a^2b^2c^2 U^2}.

It is simplest to expression :math:`\mathbf r` and :math:`\mathbf v` is
cartesian coordinates, since then there are no singularities in the
representation.

Here :math:`m` is the reduced length and :math:`K` is the Gaussian
curvature.  It's not necessary to determine :math:`m` to solve the
direct problem; however, it is an important aspect of solving the
inverse problem.

We use the ODE routines provided with Octave and MATLAB to solve these
equations.  To make the control of the error simpler, we first scale the
ellipsoid so that its middle semiaxis :math:`b = 1`; then all the
dependent variables are of order unity.  The ODE solvers take care of
picking the appropriate step size for integration.  In addition, they
allow intermediate points on the path to be found inexpensively by
polynomial interpolation.

The demonstrations ``trixial.demo(n)`` for ``n = 1:5`` and ``n = 11:15`` show
the result of solution of the direct problem of various initial
conditions.  These illustrate the distinctive properties of geodesics,
i.e., that the undulate between either lines of constant :math:`\beta`
or lines of constant :math:`\omega`.  In the limiting case, the geodesic
repeatedly returns to opposite umbilic points.

Note well: Octave is about 40 slower than MATLAB at solving the ODEs.

.. _trigeodinverse:

The inverse problem
^^^^^^^^^^^^^^^^^^^

[Panou13]_ and [Baillard13]_ both attempt to solve the inverse problem,
finding the shortest path between two points.  However, neither offers a
complete solution.  A reliable method of solving the problem is obtained
using the same basic method give by [Karney13]_ for solving the problem
on an oblate ellipoid.  The key observation is that the *cut locus* for
geodesics emanating from a given point is a segment of the line of the
opposite ellipsoidal latitude; see ``triaxial.demo(6)``.

The solution in the general case, involves starting with the point with
the large absolute latitude, varying the azimuth at this point and find
the longitude where this geodesic intersects the line of latitude for
the other point.  This makes use of the ability for the ODE solvers in
Octave/MATLAB to stop at the occurrence of certain "events".  The
azimuth can be corrected using Newton's method (this is where the
reduced length :math:`m` is needed) to find the azimuth where the
longitude matches that of the other point.

About 6 iterations are required for random pairs of points on a
terrestrial ellipsoid and the overall accuracy is probably about 1 μm.
The method is somewhat fragile in that it expects geodesics to behave in
the way dictated by Jacobi's solution; however, the ODE solver cannot
guarantee that this is so.  However by setting reasonably tight error
tolerances are set on the ODE solver and deploying some other defensive
tricks, I believe that the method work as long as the ellipsoid is not
too eccentric.  (To be safe, the ellipsoid should satisfy :math:`a/b \le
2` and :math:`b/c \le 2`.  Also avoid ellipsoids which are nearly but
not quite ellipsoids of revolution; triaxial models of the earth are
fine, but expect problems if the difference in the equatorial semiaxes
is 1 μm.)

.. _trigeodjac:

Jacobi's solution
^^^^^^^^^^^^^^^^^

I have coded up the solution to the direct problem in MATLAB using the
[Chebfun]_ package.  This allows the indefinite integrals in Jacobi's
solution to be evaluated accurately.  I do not include this
functionality in the ``triaxial`` class because

* Chebfun is not compatible with Octave;
* MATLAB's support for elliptic integrals and elliptic functions with
  modulus close to 1 is deficient --- this leads to inaccuracies for
  geodesics which graze the umbilic points.

I will reimplement the solution in the C++ version of GeographicLib.
This will make more consistent use of Fourier series (Chebfun switches
to a Chebyshev series when asked to integrate a Fourier series) and use
GeographicLib's implementation of elliptic integrals and elliptic
functions.

With this in place, the solution of the inverse problem should be
straightforward.  Jacobi does not include an expression for the reduced
length :math:`m`, so I will use some method other than Newton's for
finding the azimuth.

.. _triutils:

Utilities
---------

You can sample points (and directions) uniformly on the ellipsoid with
``cart2rand``, see [Marples+Williams23]_

The function ``horizon`` returns points on the horizon of the ellipsoid
when viewed from view poitn :math:`\mathbf V`.  These points satisfy

.. math::
   \mathbf U \cdot \mathbf V &= 0\\
   \biggl(\frac x{a^2}, \frac y{b^2}, \frac z{c^2}\biggl) \cdot \mathbf V&= 0\\
   \biggl(\frac xa, \frac yb, \frac zc\biggl) \cdot
   \biggl(\frac{V_x}a, \frac{V_y}b, \frac{V_z}c\biggl) &= 0

The first vector in the last equation gives points on a unit sphere, and
these are on the horizon of the sphere when viewed from the direction
given by the second vector.  So the ellipsoidal horizon is obtained by
computing this spherical horizon (a circle) and scaling the cartesian
components by :math:`(a, b, c)`.

.. _trirefs:

References
----------

.. [Baillard13] Baillard. `Geodesics on a triaxial ellipsoid for the
   HP-41 <https://hp41programs.yolasite.com/geod3axial.php>`__ (2013).

.. [Cayley72] Cayley, `On the geodesic lines on an ellipsoid
   <https://books.google.com/books?id=S4znAAAAMAAJ&pg=PA31>`__ (1872).

.. [Chebfun] Chebfun, `Numerical computing with functions
   <https://www.chebfun.org>`__ (2014).

.. [DLMF] Olver et al., `NIST Handbook of Mathematical Functions
   <https://dlmf.nist.gov>`__ (2010).

.. [Jacobi39] Jacobi, `Note von der geodätischen Linie auf einem
   Ellipsoid und den verschiedenen Anwendungen einer merkwürdigen
   analytischen Substitution
   <https://doi.org/10.1515/crll.1839.19.309>`__ (1839).

.. [Karney13] Karney, `Algorithms for geodesics
   <http://dx.doi.org/10.1007/s00190-012-0578-z>`__ (2013).

.. [Ligas12] Ligas, `Two modified algorithms to transform Cartesian to
   geodetic coordinates on a triaxial ellipsoid
   <http://dx.doi.org/10.1007/s11200-011-9017-5>`__ (2012).

.. [Marples+Williams23] Marples & Williams, `Patch area and uniform
   sampling on the surface of any ellipsoid
   <https://doi.org/10.1007/s11075-023-01628-4>`_ (2023).

.. [Panou13] Panou, `The geodesic boundary value problem and its
   solution on a triaxial ellipsoid
   <10.https://doi.org/2478/jogs-2013-0028>`__ (2013).

.. [Panou+Korakitis19] Panou & Korakitis, `Geodesic equations and their
   numerical solution in Cartesian coordinates on a triaxial ellipsoid
   <http://dx.doi.org/10.1515/jogs-2019-0001>`__ (2019).

.. [Panou+Korakitis21] Panou & Korakitis, `Analytical and numerical
   methods of converting Cartesian to ellipsoidal coordinates
   <http://dx.doi.org/10.1515/jogs-2020-0126>`__ (2021).

.. [Panou+Korakitis22] Panou & Korakitis, `Cartesian to geodetic
   coordinates conversion on a triaxial ellipsoid using the bisection
   method <http://dx.doi.org/10.1007/s00190-022-01650-9>`__ (2022).
