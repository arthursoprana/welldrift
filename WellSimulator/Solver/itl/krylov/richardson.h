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

#ifndef ITL_IR_H
#define ITL_IR_H

#include "itl/itl.h"

namespace itl {

  //: Preconditioned Richardson
  //
  //  This solves the unsymmetric linear system Ax = b using 
  //  Iterative Refinement (preconditioned Richardson iteration).
  //  <p>
  //  A return value of 0 indicates convergence within the
  //  maximum number of iterations (determined by the iter object).
  //  A return value of 1 indicates a failure to converge.
  //  <p>
  //  <p>
  //  See: R.S. Varga, Matrix Iterative Analysis, Automatic Computation Series,
  //  Pentice Hall Inc, Englewood Cliffs, New Jersey, 1962
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: richardson.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,add,scaled */

template < class Matrix, class Vector, class VectorB, class Preconditioner, class Iteration >
int 
richardson(const Matrix& A, Vector& x, const VectorB& b,
   const Preconditioner& M, Iteration& iter)
{
  typedef Vector TmpVec;
  TmpVec z(size(x)), r(size(x));

  itl::mult(A, itl::scaled(x, -1.), b, r);   

  while ( ! iter.finished(r) ) {
    itl::solve(M, r, z);
    itl::add(z, x);		   
    itl::mult(A, itl::scaled(x, -1.), b, r); 
    ++iter;
  }
  return iter.error_code();
}

}

#endif

