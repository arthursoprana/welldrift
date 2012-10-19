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

#ifndef ITL_BICGSTAB_H
#define ITL_BICGSTAB_H

#include "itl/itl.h"

namespace itl {

  //:  BiConjugate Gradient Stabilized
  //  A return value of 0 indicates convergence within the
  //  maximum number of iterations (determined by the iter object).
  //  A return value of 1 indicates a failure to converge.
  //  <p>
  //  <p>
  //  See: H. Van der Vorst, Bi-CGSTAB: A fast and smoothly converging variant
  //  of BiCG for the solution of nonsysmmetric linear systems, SIAM J. Sci. 
  //  Statist. Comput. 13(1992), pp. 631-644
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: bicgstab.h
  //!example: bicgstab.cc
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot,add,itl::scaled */
  template <class Matrix, class Vector, class VectorB, 
  class Preconditioner, class Iteration>
  int
  bicgstab(const Matrix& A, Vector& x, const VectorB& b,
	   const Preconditioner& M, Iteration& iter)
  {
    typedef typename itl_traits<Vector>::value_type T;
    T rho_1, rho_2, alpha, beta, omega;
    typedef Vector TmpVec;
    TmpVec p(size(x)), phat(size(x)), s(size(x)), shat(size(x)), 
           t(size(x)), v(size(x)), r(size(x)), rtilde(size(x));

    itl::mult(A, itl::scaled(x, -1.0), b, r);	  
    itl::copy(r, rtilde);	          

    while (! iter.finished(r)) {
      
      rho_1 = itl::dot(rtilde, r);
      if (rho_1 == T(0.)) {
	iter.fail(2, "bicg breakdown #1");
	break;
      }
    
      if (iter.first())
	itl::copy(r, p);
      else {
	if (omega == T(0.)) {
	  iter.fail(3, "bicg breakdown #2");
	  break;
	}
      
	beta = (rho_1 / rho_2) * (alpha / omega);
      
	itl::add(itl::scaled(v, -omega), p); 
	itl::add(r, itl::scaled(p, beta), p);      
      }
      itl::solve(M, p, phat);
      itl::mult(A, phat, v);	
      alpha = rho_1 / itl::dot(v, rtilde);
      itl::add(r, itl::scaled(v, -alpha), s);
    
      if (iter.finished(s)) {
	itl::add(itl::scaled(phat, alpha), x); 
	break;
      }
    
      itl::solve(M, s, shat);	
      itl::mult(A, shat, t);	
      omega = itl::dot(t, s) / itl::dot(t, t);
    
      itl::add(itl::scaled(phat, alpha), x); 
      itl::add(itl::scaled(shat, omega), x);
      itl::add(s, itl::scaled(t, -omega), r); 
    
      rho_2 = rho_1;
    
      ++iter;
    }

    return iter.error_code();
  }

}

#endif
