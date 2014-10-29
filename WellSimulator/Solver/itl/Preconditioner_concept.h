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

//: Preconitioner
//  The Preconditioner object performs a preconditioning operation based
//  on vector x and stores the result in vector z.  The trans solve()
//  method only need be defined when the preconditioner is used with an
//  iterative solver that requires it.
// 
//!category: itl,utilities
//!component: concept
//!models: preconditioner
//!notation: M - A Preconditioner object
//!notation: x, z - A Vector
concept Preconditioner {
  //: solve M z = x
  //!exp: solve(M, x, z);
  solve(const Preconditioner& M, const Vecter& x,  Vector& z);
};


concept HermitianPreconditioner : concept Preconditioner {
  //: solve M^T z = x
  //!exp: trans_solve(M, x, z); 
  trans_solve(const Preconditioner& M, const Vector& x, Vector& z);
};
