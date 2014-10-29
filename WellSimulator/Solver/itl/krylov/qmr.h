// -*- c++ -*-
//
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//                                                                           

#ifndef ITL_QMR_H
#define ITL_QMR_H


#include "itl/itl.h"

namespace itl {

  //: Quasi-Minimal Residual
  //
  //  This routine solves the unsymmetric linear system Ax = b using the
  //  Quasi-Minimal Residual method.
  //<table align=center border=1>
  // <tr><td> return value </td>   <td>   meaning </td> </tr>
  // <tr><td>      0   </td><td>   convergence within maximum iterations </td> </tr>
  // <tr><td>      1   </td><td>     no convergence after maximum iterations</td> </tr>
  //  <tr><td>     2  </td><td>      breakdown in       rho </td> </tr>
  //  <tr><td>     3  </td><td>      breakdown in       beta </td> </tr>
  //  <tr><td>     4  </td><td>      breakdown in       gamma</td> </tr>
  //  <tr><td>     5  </td><td>      breakdown in       delta</td> </tr>
  //  <tr><td>     6  </td><td>      breakdown in       ep</td> </tr>
  //  <tr><td>     7  </td><td>      breakdown in       xi</td> </tr>
  // </table>
  //See: R. W. Freund and N. M. Nachtigal, A quasi-minimal residual method for 
  //non-Hermitian linear systems, Numerical Math., 60(1991), pp. 315-339
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: qmr.h
  //!example: qmr.cc
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot,add,scaled,two_norm,trans_mult,scale */

  template <class Matrix, class Vector, class VectorB, class Precond1, class Precond2, class Iteration>
  int 
  qmr(const Matrix &A, Vector &x, const VectorB &b, const Precond1 &M1,
      const Precond2 &M2, Iteration& iter)
  {
    typedef typename itl_traits<Vector>::value_type value_type;
    value_type delta(0), ep(0), beta(0), rho_1(0), gamma_1(0), theta_1(0);
  
    typedef Vector TmpVec;
    TmpVec r(size(x)), v_tld(size(x)), y(size(x)), w_tld(size(x)), z(size(x));
    TmpVec v(size(x)), w(size(x)), y_tld(size(x)), z_tld(size(x));
    TmpVec p(size(x)), q(size(x)), p_tld(size(x)), d(size(x)), s(size(x));
  
    itl::mult(A, itl::scaled(x, -1.0), b, r);
    itl::copy(r, v_tld);

    itl::solve(M1, v_tld, y);
    value_type rho = itl::two_norm(y);

    itl::copy(r, w_tld);
    itl::trans_solve(M2, w_tld, z);
    value_type xi = itl::two_norm(z);
  
    value_type gamma = 1.0, eta = -1.0, theta = 0.0;
  
    while (! iter.finished(r)) {
    
      if (rho == 0.0) {
	iter.fail(2, "qmr breakdown");
	break;
      }
      if (xi == 0.0) {
	iter.fail(7, "qmr breakdown");
	break;
      }

      itl::copy(itl::scaled(v_tld, 1./rho), v);
      itl::scale(y, 1./rho);

      itl::copy(itl::scaled(w_tld, 1./xi), w);
      itl::scale(z, 1./xi);

      delta = itl::dot(z, y);
      if (delta == 0.0) {
	iter.fail(5, "qmr breakdown");
	break;
      }

      itl::solve(M2, y, y_tld);		
      itl::trans_solve(M1, z, z_tld);

      if (iter.first()) {
	itl::copy(y_tld, p);
	itl::copy(z_tld, q);
      } else {
	itl::add(y_tld, itl::scaled(p, -(xi  * delta / ep)), p);
	itl::add(z_tld, itl::scaled(q, -(rho * delta / ep)), q);
      }
    
      itl::mult(A, p, p_tld);

      ep = itl::dot(q, p_tld);
      if (ep == 0.0) {
	iter.fail(6, "qmr breakdown");
	break;
      }

      beta = ep / delta;
      if (beta == 0.0) {
	iter.fail(6, "qmr breakdown");
	break;
      }

      itl::add(p_tld, itl::scaled(v, -beta), v_tld);
      itl::solve(M1, v_tld, y);

      rho_1 = rho;
      rho = itl::two_norm(y);

      itl::trans_mult(A, q, w_tld);
      itl::add(w_tld, itl::scaled(w, -beta), w_tld);
      itl::trans_solve(M2, w_tld, z);

      xi = itl::two_norm(z);

      gamma_1 = gamma;
      theta_1 = theta;

      theta = rho / (gamma_1 * beta);
      gamma = 1.0 / sqrt(1.0 + theta * theta);

      if (gamma == 0.0){
	iter.fail(4, "qmr breakdown");
	break;
      }

      eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);

      if (iter.first()) {
	itl::copy(itl::scaled(p, eta), d);
	itl::copy(itl::scaled(p_tld, eta), s);
      } else {
	value_type tmp = (theta_1 * theta_1 * gamma * gamma);
	itl::add(itl::scaled(p, eta), itl::scaled(d, tmp), d);
	itl::add(itl::scaled(p_tld, eta), itl::scaled(s, tmp), s);
      }
      itl::add(d, x);
      itl::add(itl::scaled(s, -1.), r);

      ++iter;
    }

    return iter.error_code();
  }


}

#endif 

