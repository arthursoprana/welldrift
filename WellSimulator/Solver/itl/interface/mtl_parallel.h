// -*- c++ -*-
//===========================================================================
//  CVS Information:                                                         
//                                                                           
//     $RCSfile: mtl_parallel.h,v $  $Revision: 1.3 $  $State: Exp $ 
//     $Author: llee $  $Date: 2001/10/26 15:02:15 $ 
//     $Locker:  $ 
//---------------------------------------------------------------------------
//                                                                           
// DESCRIPTION                                                               
//   This is the file contain the implementations of basic operations.
//   It provides mtl implementations of those operations. Therefore,
//   It requires mtl package. The tested mtl release is 2.1.2
//
//   It provides the parallel itl interface using mtl package.
//   Matrix is mtl row-wise matrices
//   Vector is mtl dense1D
//   Preconditioners: block_ilu, identity_preconditioner
//   Krylov subspace methode: gmres, tfqmr, bicgstab, cgs
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
// $Log: mtl_parallel.h,v $
// Revision 1.3  2001/10/26 15:02:15  llee
// *** empty log message ***
//
// Revision 1.2  2001/10/19 05:23:06  llee
// *** empty log message ***
//
//                                                                           
//===========================================================================

#ifndef PARALLEL_INTERFACE_H
#define PARALLEL_INTERFACE_H

#include <itl/itl_config.h>
#include "mtl/dense1D.h"
#include "mtl/mtl.h"
#include "mtl/matrix.h"

#include <utility>
#include "itl/interface/detail/Manager.h"
#include "itl/interface/detail/parallel_dot.h"
#include "itl/interface/detail/parallel_matrix_v2.h"
#include <itl/itl_tags.h>

namespace itl {
  
  template<class Vector>
  struct itl_traits {
    //This type will be used for temp vectors inside of itl algorithm 
    typedef referencing_object_tag vector_category;
    typedef typename Vector::value_type value_type;
    typedef typename Vector::size_type  size_type;
  };


  template <class Matrix, class VectorX, class VectorB>
  inline 
  void mult(const parallel_matrix<Matrix>& A, const VectorX& x,
	    const VectorB& b) {
    A.process_vector(x);
    A.mult(const_cast<VectorB&>(b));
  }

  template <class Matrix, class VectorX, class VectorY, class VectorB>
  inline 
  void mult(const parallel_matrix<Matrix>& A, const VectorX& x,
	    const VectorY& y, const VectorB& b) {
    A.process_vector(x);
    A.mult(const_cast<VectorB&>(b));
    mtl::add(y, b);
  }

  template <class VecA, class VecB>
  inline void copy(const VecA& a, const VecB& b) {
    mtl::copy(a, const_cast<VecB&>(b));
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

  template <class Vec, class T>
  inline void scale(const Vec& v, T t) {
    mtl::scale(const_cast<Vec&>(v), t);
  }

  template <class Scalable, class T>
  inline typename Scalable::scaled_type
  scaled(const Scalable& v, T t) {
    return mtl::scaled(v, t);
  }

  template <class VectorX, class VectorY>
  inline double  dot_conj(const VectorX& a, const VectorY& b) {
    double ret = parallel_dot_conj(a, b);    
    return ret;
  }

  template <class VectorX, class VectorY>
  inline double  dot(const VectorX& a, const VectorY& b) {
    double ret = parallel_dot(a, b);
    return ret;
  }

  template <class Vector>
  inline double two_norm(const Vector& v) {
    double ret = sqrt(itl::dot_conj(v, v));
    return ret;
  }
  
  template <class Vector>
  inline int size(const Vector& x) {
    return x.size();
  }

  template <class Vector>
  inline void resize(Vector& x, int sz) {
    x.resize(sz);
  }

  template <class Hessenberg, class Vec>
  inline void upper_tri_solve(const Hessenberg& hh, Vec& rs, int i) {
    mtl::tri_solve(mtl::tri_view<mtl::upper>()(hh.sub_matrix(0, i, 0, i)), rs);
  }

  template <class T>
  struct internal_matrix_traits {
    typedef typename mtl::matrix<T, mtl::rectangle<>,
				 mtl::dense<>,
				 mtl::column_major>::type Matrix;
  };

}

#endif
