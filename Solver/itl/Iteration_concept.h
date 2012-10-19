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

//: Iteration
//
//  The Iteration object calculates whether the solution has reached the
//  desired accuracy, or whether the maximum number of iterations has
//  been reached. The method finished() checks both convergence and
//  number of iterations. The method converged() only checks
//  convergence. The error code() method is used to determine the return
//  value for the this iterative solver function. The first() method is
//  used to determine the first iteration of the loop.
//  <p>
//  <p>
//  For all algorithms, if the error_code() is 0, it suggests the algorithm 
//  converges. Otherwise, if the error_code() returns 1, it means the maximum 
//  number of iteration has been reached but the desired accuacy is not reached.
//  For other return codes, see the respective document.
//!category: itl,utilities
//!component: concept
concept Iteration {

  //: Testing if stop criteria is satisfied
  //!exp: finish( Vector );
  bool finish(const Vector& r);
  //: Testing if stop criteria is satisfied for the case of qusi residual case.
  //!exp: finish( TrivialType );
  bool finish(const value_type& r); /*used in TFQMF*/
  //: Testing if it is converged
  //!exp: converged( Vector );
  bool converged(const Vector& r);
  //: Testing if it is converged for the case of qusi residual case.
  //!exp: converged( TrivialType );
  bool converged(const value_type& r); /*used in TFQMF*/
  //: to increment number of iteration by one
  //!exp: ++Iter;
  void operator++();
  //: to Check if this is the first iteration
  //!exp: Iter.first();
  bool first();
  //: to return error_code. Zero means success.
  //!exp: Iter.error_code();
  int error_code();
  //: to return number of iteration to be performed.
  //!exp: Iter.iteration();
  int iteration();
  //: to return residual or qusi-residual
  //!exp: Iter.resid();
  magnitude_type resid();
  //: to return tolerance.
  //!exp: Iter.tol();
  magnitude_type tol();
  //: set fail reason once fail happens during iteration
  //!exp: fail( int );
  void fail(int);
  //: set fail reason once fail happens during iteration
  //!exp: fail( int, string );
  void fali(int, const std::string&);
};
