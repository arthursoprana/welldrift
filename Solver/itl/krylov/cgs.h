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

#ifndef ITL_CGS_H
#define ITL_CGS_H

#include "itl/itl.h"

namespace itl {

  //: Conjugate Gradient Squared
  //
  //  Solves the unsymmetric linear system Ax = b 
  //  using the Conjugate Gradient Squared method.
  //  <p>
  //  A return value of 0 indicates convergence within the
  //  maximum number of iterations (determined by the iter object).
  //  A return value of 1 indicates a failure to converge.
  //  <p>
  //  <p>
  //See: P. Sonneveld, CGS, a fast Lanczos-type solver for nonsymmetric linear,
  //systems, SIAM, J.Sci. Statist. Comput., 10(1989), pp. 36-52
  //!category: itl,algorithms
  //!component: function
  //!definition: cgs.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: VectorX - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot,add,itl::scaled */
template <class Matrix, class Vector, class VectorB, 
          class Preconditioner, class Iteration>
int
cgs(const Matrix& A, Vector& x, const VectorB& b,
    const Preconditioner& M, Iteration& iter)
{
  typename itl_traits<Vector>::value_type rho_1, rho_2, alpha, beta;
  typedef Vector TmpVec;
  TmpVec p(size(x)), phat(size(x)), q(size(x));
  TmpVec qhat(size(x)), vhat(size(x)), u(size(x));
  TmpVec uhat(size(x)), r(size(x)), rtilde(size(x));
  
  itl::mult(A, itl::scaled(x, -1.0), b, r);	  
  itl::copy(r, rtilde);

  while (! iter.finished(r)) {
    rho_1 = itl::dot(rtilde, r);

    if (rho_1 == 0.) {
      iter.fail(2, "cgs breakdown");
      break;
    }

    if (iter.first()) {
      itl::copy(r, u);
      itl::copy(u, p);
    } else {
      beta = rho_1 / rho_2;

      itl::add(r, itl::scaled(q, beta), u); 
      itl::add(q, itl::scaled(p, beta), p);   
      itl::add(u, itl::scaled(p, beta), p);
    }
    itl::solve(M, p, phat);
    itl::mult(A, phat, vhat);
    alpha = rho_1 / itl::dot(rtilde, vhat);
    itl::add(u, itl::scaled(vhat, -alpha), q);
    
    //add(u, q, u);
    itl::add(q, u);
    itl::solve(M, u, uhat);

    itl::add(x, itl::scaled(uhat, alpha), x);
    itl::mult(A, uhat, qhat);
    itl::add(r, itl::scaled(qhat, -alpha), r);

    rho_2 = rho_1;

    ++iter;
  }
  
  return iter.error_code();
}

}

#endif
