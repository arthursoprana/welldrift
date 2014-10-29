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
#ifndef ITL_GCR_H
#define ITL_GCR_H

#include "itl/itl.h"
#include "itl/number_traits.h"
#include <vector>

namespace itl {

  //: Generalized Conjugate Residual
  //
  //   This solve the linear system Ax = b with restarted Preconditioned 
  //   Generalized Conjugate Residual Algorithm.
  //   <p>
  //   A return value of 0 indicates convergence within the
  //   maximum number of iterations (determined by the iter object).
  //   A return value of 1 indicates a failure to converge.
  //   <p>
  //<p>
  //See: Y. Saad, Iterative Methods for Sparse Linear System, PWS Publishing
  //   Company, 1996
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: gcr.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot_conj,add,itl::scaled, two_norm */
template < class Matrix, class Vector, class VectorB, class Preconditioner, class Iteration >
int 
gcr(const Matrix& A, Vector& x, const VectorB& b,
    const Preconditioner& M, int m, Iteration& outer)
{
  typedef typename itl_traits<Vector>::value_type value_type;
  typedef typename itl::number_traits<value_type>::magnitude_type Real;
  typedef typename internal_matrix_traits<value_type>::Matrix gmx;
  typedef Vector TmpVec;

  gmx  p(size(x), m+1), w(size(x), m+1);

  std::vector<value_type> beta(m+1);

  TmpVec r(size(x)), q(size(x)), u(size(x));

  itl::mult(A, itl::scaled(x, -1.), b, u);

  itl::solve(M, u, r);

  value_type alpha;

  Real normr = itl::two_norm(r);
  
  while (! outer.finished(normr)) {
    
    Iteration inner(outer.normb(), m, outer.tol(), outer.atol());
    
    itl::copy(itl::scaled(r, 1./normr), p[0]);
    
    int j = 0;
    
    while (! inner.finished(r) ) {
      itl::mult(A, p[j], u);
      itl::solve(M, u, w[j]);
      
      beta[j] = itl::dot_conj(w[j], w[j]);
      alpha = itl::dot_conj(w[j], r) / beta[j];
      
      itl::add(x, itl::scaled(p[j], alpha), x);
      itl::add(r, itl::scaled(w[j], -alpha), r);
      
      itl::mult(A, r, u);
      itl::solve(M, u, q);

      for (int i=0; i<=j; i++)
	itl::add(r, itl::scaled(p[i], 
				-itl::dot_conj(w[i], q)/beta[i]), 
	    p[j+1]);

      ++inner; ++outer; ++j;
    } 

     normr = itl::two_norm(r);
  }

  return outer.error_code();
}

}

#endif
