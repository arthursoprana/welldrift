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

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include "mtl/matrix.h"
#include "mtl/mtl.h"

/* Rich, the documentation for this needs help */

namespace itl {
  
//:  Preconditioner application class
//
// For left or right preconditioner.
//
//!tparam: Matrix1 - An MTL Matrix type
//!tparam: Matrix2 - Another MTL Matrix type
//!tparam: Lower - Uplo for the lower triangle (can be unit)
//!tparam: Upper - Uplo for the upper triangle (can be unit)
template <class Matrix1, class Matrix2,  int Lower, int Upper >
class preconditioner {
public:
  inline preconditioner() {}
  inline preconditioner(const Matrix1& L, const Matrix2& U)
    :LA(L), 
     UA(U) {}
  ///  //LUz = x ==>  Ly=x, Uz=y 
  template <class VecX, class VecZ>
  inline void solve(const VecX& x, const VecZ& z) const {
    mtl::copy(x, const_cast<VecZ&>(z));      
    mtl::tri_solve(LA, const_cast<VecZ&>(z));
    mtl::tri_solve(UA, const_cast<VecZ&>(z));
  }
  /////(LU)^T z = x ==>U^T y = x //L^T z = y
  template <class VecX, class VecZ>
  inline void trans_solve(const VecX& x, const VecZ& z) const {
    mtl::copy(x, const_cast<VecZ&>(z));
    mtl::tri_solve(mtl::trans(UA), const_cast<VecZ&>(z)); 
    mtl::tri_solve(mtl::trans(LA), const_cast<VecZ&>(z));
  }

private:
  //it is unit_lower triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix1, Lower>::type LA;

  //it is upper triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix2, Upper>::type UA;

};


///for split preconditioner, the left part
template <class Matrix1, class Matrix2, int Lower, int Upper>
class preconditioner1 {

public:
  inline preconditioner1() {}
  inline preconditioner1(const Matrix1& L, const Matrix2& U)
    :LA(L), 
     UA(U) {}
  ///
  template <class VecX, class VecZ>
  inline void solve(const VecX& x, const VecZ& z) const {  //LUz = x ==>  Ly=x, Uz=y 
    mtl::copy(x, const_cast<VecZ&>(z));                
    mtl::tri_solve(LA, const_cast<VecZ&>(z));
  }
  ///
  template <class VecX, class VecZ>
  inline void trans_solve(const VecX& x, const VecZ& z) const {//(LU)^T z = x ==>U^T y = x 
    mtl::copy(x, const_cast<VecZ&>(z));                                           //L^T z = y
    mtl::tri_solve(mtl::trans(LA), const_cast<VecZ&>(z)); 
  }

private:

  //it is unit_lower triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix1, Lower>::type LA;

  //it is upper triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix2, Upper>::type UA;

};


///for split preconditioner, the right part
template <class Matrix1, class Matrix2, int Lower, int Upper>
class preconditioner2 {

public:
  inline preconditioner2() {}
  inline preconditioner2(const Matrix1& L, const Matrix2& U)
    :LA(L), 
     UA(U) {}
  ///
  template <class VecX, class VecZ>
  inline void solve(const VecX& x, const VecZ& z) const {  //LUz = x ==>  Ly=x, Uz=y 
    mtl::copy(x, const_cast<VecZ&>(z));                
    mtl::tri_solve(UA, const_cast<VecZ&>(z));
  }
  ///
  template <class VecX, class VecZ>
  inline void trans_solve(const VecX& x, const VecZ& z) const {//(LU)^T z = x ==>U^T y = x 
    mtl::copy(x, const_cast<VecZ&>(z));                                           //L^T z = y
    mtl::tri_solve(mtl::trans(UA), const_cast<VecZ&>(z)); 
  }

private:

  //it is unit_lower triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix1, Lower>::type LA;

  //it is upper triangular matrix from real preconditioner
  typename mtl::triangle_view<Matrix2, Upper>::type UA;

};


//@}

}

#endif
