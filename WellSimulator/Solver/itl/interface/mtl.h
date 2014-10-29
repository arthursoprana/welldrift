// -*- c++ -*-
//===========================================================================
//  CVS Information:                                                         
//                                                                           
//     $RCSfile: mtl.h,v $  $Revision: 1.5 $  $State: Exp $ 
//     $Author: llee $  $Date: 2001/10/26 15:02:15 $ 
//     $Locker:  $ 
//---------------------------------------------------------------------------
//                                                                           
// DESCRIPTION                                                               
//   This is the file contain the implementations of basic operations.
//   It provides mtl implementations of those operations. Therefore,
//   It requires mtl package. The tested mtl release is 2.1.2
//                                                                           
//---------------------------------------------------------------------------
//                                                                           
// LICENSE AGREEMENT                                                         
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
//---------------------------------------------------------------------------
//                                                                           
// REVISION HISTORY:                                                         
//                                                                           
// $Log: mtl.h,v $
// Revision 1.5  2001/10/26 15:02:15  llee
// *** empty log message ***
//
// Revision 1.4  2001/10/23 16:15:01  llee
// *** empty log message ***
//
// Revision 1.3  2001/10/18 21:37:06  llee
// *** empty log message ***
//
// Revision 1.2  2001/10/18 15:44:32  llee
// compile and test for examples
//
// Revision 1.1  2001/10/18 14:08:32  llee
// re-organize the directory structures
//
// Revision 1.5  2000/07/26 17:14:45  llee1
// *** empty log message ***
//
//                                                                           
//===========================================================================

#ifndef ITL_MTL_INTERFACE_H
#define ITL_MTL_INTERFACE_H

/*
  This is to provide an interface for ITL to use MTL. 
  
  Matrix: mtl matrices
  Vector: mtl dense1D
  Preconditioners: all preconditioners in itl/preconditioners
 */

#include <itl/itl_config.h>
#include <mtl/dense1D.h>
#include "mtl/matrix.h"
#include "mtl/mtl.h"
#include <complex>
#include <itl/itl_tags.h>
#include <itl/number_traits.h>

namespace itl {
  
  //: The vector type used inside of the ITL routines for work space
  template <class Vec>
  struct itl_traits {
    typedef referencing_object_tag vector_category;
    typedef typename Vec::value_type value_type;
    typedef typename Vec::size_type size_type;
  };

  template <class Matrix>
  inline typename Matrix::size_type nrows(const Matrix& A)
  { return A.nrows(); }

  template <class Vec>
  inline 
  typename itl::number_traits< typename Vec::value_type >::magnitude_type
  two_norm(const Vec& v) {
    return mtl::two_norm(v);
  }

  //deal with the case when b is a handle, 
  template <class VecA, class VecB>
  inline void copy(const VecA& a, const VecB& b) {
    mtl::copy(a, const_cast<VecB&>(b));
  }
  
  template <class Matrix, class VecX, class VecY>
  inline void mult(const Matrix& A, const VecX& x, const VecY& y) {
    mtl::mult(A, x, const_cast<VecY&>(y));
  }

  template <class Matrix, class VecX, class VecY, class VecZ>
  inline void mult(const Matrix& A, const VecX& x, const VecY& y, 
		   const VecZ& z) {
    mtl::mult(A, x, y, const_cast<VecZ&>(z));
  }

  template <class VecA, class VecB>
  inline typename VecA::value_type dot(const VecA& a, const VecB& b) {
    return mtl::dot(a, b);
  }
  
  template <class VecA, class VecB>
  inline typename VecA::value_type  dot_conj(const VecA& a, const VecB& b) {
    return mtl::dot_conj(a, b);
  }
  

  template <class VecX, class VecY>
  inline void add(const VecX& x, const VecY& y) {
    mtl::add(x, const_cast<VecY&>(y));
  }

  template <class VecX, class VecY, class VecZ>
  inline void add(const VecX& x, const VecY& y, const VecZ& z) {
    mtl::add(x, y, const_cast<VecZ&>(z));
  }

  template <class VecX, class VecY, class VecZ, class VecR>
  inline void add(const VecX& x, const VecY& y, const VecZ& z, const VecR& r) {
    mtl::add(x, y, z, const_cast<VecR&>(r));
  }

  template <class Scalable, class T>
  inline typename Scalable::scaled_type
  scaled(const Scalable& v, T t) {
    return mtl::scaled(v, t);
  }

  template <class Vec, class T>
    inline void scale(const Vec& v, T t) {
    mtl::scale(const_cast<Vec&>(v), t);
  }

  template <class Vector>
  inline typename Vector::size_type
  size(const Vector& x) 
  {
    return x.size();
  }


  template <class Vector, class Size>
  inline void
  resize(Vector& x, const Size& sz) 
  {
    x.resize(sz);
  }


  template <class Matrix, class VecX, class VecY>
  inline void trans_mult(const Matrix& A, const VecX& x, const VecY& y) {
    mtl::mult(mtl::trans(A), x, const_cast<VecY&>(y));
  }

  //used inside of GCR and GMRES algorithm

  template <class T>
  struct internal_matrix_traits {
    typedef typename mtl::matrix<T, mtl::rectangle<>,
				 mtl::dense<>,
				 mtl::column_major>::type Matrix;
  };

  
  template <class Hessenberg, class Vec>
  inline void upper_tri_solve(const Hessenberg& hh, Vec& rs, int i) {
    mtl::tri_solve(mtl::tri_view<mtl::upper>()(hh.sub_matrix(0, i, 0, i)), rs);
  }

}

//for classical gram schmidt 
#include "itl/interface/detail/mtl_classical_gram_schmidt.h"

#endif /*ITL_MTL_INTERFACE_H*/
