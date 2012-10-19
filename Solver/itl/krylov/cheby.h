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

#ifndef ITL_CHEBY_H
#define ITL_CHEBY_H

#include "itl/itl.h"

namespace itl {

  //: Chebyshev Iteration
  //This solves the unsymmetric  linear system Ax = b using 
  //the Preconditioned Chebyshev Method.
  //<p>
  //A return value of 0 indicates convergence within the
  //maximum number of iterations (determined by the iter object).
  //A return value of 1 indicates a failure to converge.
  //<p>
  //See: T. Manteuffel, The Chebyshev iteration for nonsysmmetric linear systems
  //Numer. Math. 28(1977), pp. 307-327
  //G. H. Golub and C. F. Van Loan, Matrix Computations, The Johns Hopkins 
  //University Press, Baltimore, Maryland, 1996 
  //
  //!category: itl,algorithms 
  //!component: function 
  //!definition: cheby.h 
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods 
  //!tparam: Vector  - Vector  
  //!tparam: VectorB - Vector 
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner. 
  //!tparam: Iteration - Controls the stopping criteria 

  /* required operations: mult,copy,add,scaled */

template < class Matrix, class Vector, class VectorB, 
           class Preconditioner, class Iteration >
int 
cheby(const Matrix &A, Vector &x, const VectorB &b,
      const Preconditioner &M, Iteration& iter,
      typename Vector::value_type eigmin, 
      typename Vector::value_type eigmax)
{

  typedef typename itl_traits<Vector>::value_type Real;
  Real alpha, beta, c, d;
  typedef Vector TmpVec;
  TmpVec p(size(x)), q(size(x)), z(size(x)), r(size(x));

  itl::mult(A, itl::scaled(x, -1.0), b, r);

  c = (eigmax - eigmin) / 2.0;
  d = (eigmax + eigmin) / 2.0;

  while ( ! iter.finished(r) ) {
    itl::solve(M, r, z);         

    if ( iter.first() ) {
      itl::copy(z, p);          
      alpha = 2.0 / d;
    } else {
      beta = c * alpha / 2.0;    
      beta = beta * beta;
      alpha = 1.0 / (d - beta);  
      itl::add(z, itl::scaled(p, beta), p);
    }

    itl::mult(A, p, q);
    itl::add(x, itl::scaled(p, alpha), x);
    itl::add(r, itl::scaled(q, -alpha), r);

    ++iter;
  }

  return iter.error_code();
}

}

#endif
