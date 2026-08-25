#ifndef OPT_MAX_FN_H
#define OPT_MAX_FN_H

#include <cfloat>
#include <cmath>

/*
below is an adapted version of `Brent_fmin` taken from `src/library/stats/optimize.c` of the R source repo @ trunk 89985
it is only very lightly adapted to:
- remove unneeded interop code w/ R
- adapt to C++ conventions
- maximize f instead of minimizing
- minor reformatting (by clang-tidy)

copyright declaration below:
  R : A Computer Language for Statistical Data Analysis
  Copyright (C) 1998--2025  The R Core Team
  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
  Copyright (C) 2003-2004  The R Foundation

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, a copy is available at
  https://www.R-project.org/Licenses/
*/

const double DEFAULT_TOL = std::pow(DBL_EPSILON, 0.25);
template <class Fn1> inline double optimize_fmax(const Fn1 &f, const double x_min = -30., const double x_max = 30., const double tol = DEFAULT_TOL) {
  /*  c is the squared inverse of the golden ratio */
  const double c = (3. - std::sqrt(5.)) * .5;

  /* Local variables */
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
  eps = DBL_EPSILON;
  tol1 = eps + 1.; /* the smallest 1.000... > 1 */
  eps = std::sqrt(eps);

  a = x_min;
  b = x_max;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.; /* -Wall */
  e = 0.;
  fx = -f(x);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  /*  main loop starts here ----------------------------------- */

  for (;;) {
    xm = (a + b) * .5;
    tol1 = eps * std::abs(x) + tol3;
    t2 = tol1 * 2.;

    /* check stopping criterion */

    if (std::abs(x - xm) <= t2 - (b - a) * .5)
      break;
    p = 0.;
    q = 0.;
    r = 0.;
    if (std::abs(e) > tol1) { /* fit parabola */

      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.)
        p = -p;
      else
        q = -q;
      r = e;
      e = d;
    }

    if (std::abs(p) >= std::abs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

      if (x < xm)
        e = b - x;
      else
        e = a - x;
      d = c * e;
    } else { /* a parabolic-interpolation step */

      d = p / q;
      u = x + d;

      /* f must not be evaluated too close to x_min or x_max */

      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm)
          d = -d;
      }
    }

    /* f must not be evaluated too close to x */

    if (std::abs(d) >= tol1)
      u = x + d;
    else if (d > 0.)
      u = x + tol1;
    else
      u = x - tol1;

    fu = -f(u);

    /*  update  a, b, v, w, and x */

    if (fu <= fx) {
      if (u < x)
        b = x;
      else
        a = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  /* end of main loop */

  return x;
}

#endif