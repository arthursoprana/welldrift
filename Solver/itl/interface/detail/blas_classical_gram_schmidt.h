// -*- c++ -*-
//===========================================================================
//  CVS Information:                                                         
//                                                                           
//     $RCSfile: blas_classical_gram_schmidt.h,v $  $Revision: 1.1 $  $State: Exp $ 
//     $Author: llee $  $Date: 2001/10/25 16:40:45 $ 
//     $Locker:  $ 
//---------------------------------------------------------------------------
//                                                                           
// DESCRIPTION                                                               
//
// This is a module for the last template in gmres function. It
// provides the basis stored in MTL dense column-wise matrix
// format. The class Gram-Schmidt orthogonalization with iterative
// refinement is performed.  Matrix-vector operations consist of the
// lowel level linear algebra operations.
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
// $Log: blas_classical_gram_schmidt.h,v $
// Revision 1.1  2001/10/25 16:40:45  llee
// *** empty log message ***
//
// Revision 1.3  2000/07/26 16:55:08  llee1
// *** empty log message ***
//
// Revision 1.2  2000/07/26 14:52:22  llee1
// *** empty log message ***
//
// Revision 1.1  2000/07/26 04:19:35  llee1
// *** empty log message ***
//
// Revision 1.2  2000/07/25 22:57:54  llee1
// *** empty log message ***
//
//                                                                           
//===========================================================================


#ifndef ITL_CLASSICAL_GRAM_SCHMIDT_ORTHOGONALIZATION_H
#define ITL_CLASSICAL_GRAM_SCHMIDT_ORTHOGONALIZATION_H

#include <itl/itl.h>
#include <itl/itl_tags.h>
#include <itl/number_traits.h>

namespace itl {
  //in parallel envirnment, this is an all-reduce operation
  template<class VecX, class VecY>
  inline void all_reduce(const VecX& x, const VecY& yy, int k) {
    VecY& y = const_cast<VecY&>(yy);
    for (int i=0; i<k; i++)
      y[i] += x[i];
  }


  //This is to use mtl dense column-wise matrix for orthogonal vectors
  template <class Vec>
  class classical_gram_schmidt {
    typedef typename itl_traits<Vec>::value_type value_type;
    typedef typename internal_matrix_traits<value_type>::Matrix ColumnMatrix;
  public:
    typedef size_t size_type;

    template <class Size>
    classical_gram_schmidt(int restart, const Size& s)
      : V(s, restart+1), ss(restart+1), w(s) {}
    
    
    typedef typename ColumnMatrix::OneD OneD;
    OneD operator[](size_type i) const {
      return V[i];
    }

    ColumnMatrix V;
    Vec ss;
    Vec w;
  };


    template <class Vec, class VecHi>
    void orthogonalize(classical_gram_schmidt<Vec>& KS,
		       const VecHi& _Hi, int i) {
      /*classical Gram-Schmidt with 1 iterative refinement.
	This is faster than modified Gram-Schmidt since here 
	we use BLAS2 operations instead of BLAS1 operations.
      */
      VecHi& Hi = const_cast<VecHi&>(_Hi);
      char trans = 'T' ;
      char not_trans = 'N';
      int m = KS.V.nrows();
      double alpha = 1.0;
      double beta = 0.0;
      int one = 1;
      int i1 = i + 1; 

      for (int k=0; k<=i; k++) Hi[k] = 0.0;
      
      for (int ii=0; ii<2; ii++) {
	//V is a N x m+1  column-wise matrix
	
	SP_NAME(dgemv)(&trans, &m, &i1, &alpha, KS.V.get_val(),
		       &m, KS[i+1].data(), &one, &beta, get_data(KS.ss), &one);
	
	//aggregate ss
	itl::all_reduce(KS.ss, Hi, i+1);
	
	//itl::mult(V.V.sub_matrix(0, size(V[i+1]), 0, i+1), V.ss, V.w);

	SP_NAME(dgemv)(&not_trans, &m, &i1, &alpha, KS.V.get_val(),
		       &m, get_data(KS.ss), &one, &beta, get_data(KS.w), &one);

	itl::add(itl::scaled(KS.w, -1.0), KS[i+1]);
      }

    }
      
    template <class Vec, class VecS, class VecX>
    void combine(classical_gram_schmidt<Vec>& KS,
		 const VecS& s, VecX& x, int i) {
      char trans = 'N';
      int one = 1;
      double alpha = 1.0;
      double beta  = 1.0;
      int m = KS.V.nrows();
      //itl::mult(KS.V.sub_matrix(0, size(x), 0, i), s, x, x);
      SP_NAME(dgemv)(&trans, &m, &i, &alpha, KS.V.get_val(),
		     &m, get_data(s), &one, &beta, get_data(x), &one);
    }
  
}

#endif
