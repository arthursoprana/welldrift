//
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
//

#ifndef ITL_SSOR_H
#define ITL_SSOR_H

#include <vector>
#include <functional>
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/meta_if.h"
#include "mtl/meta_equal.h"
#include "mtl/matrix.h"

namespace itl {

  using mtl::ROW_MAJOR; using mtl::COL_MAJOR;
  using mtl::generators_error;

  //:   SSOR preconditioner.
  //<codeblock>
  // Usage:
  //     Matrix A;
  //     SSOR<Matrix> precond(A);
  //     qmr(A, x, b, precond.left(), precond.right(), iter);
  //     cg(A, x, b, precond(), iter);
  //</codeblock>
  // Matrix has to be unsymmetric. 
  // For symmetric one, use incomlete cholesky.
  // Notes: The idea under a concrete Preconditioner such 
  // as Incomplete Cholesky is to create a Preconditioner
  // object to use in iterative methods. 
  //
  //!definition: ssor.h
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix - Matrix 
  //!example: ssor.cc
  //
  template <class Matrix>
  class SSOR {
    typedef typename Matrix::value_type T;
    typedef typename Matrix::orientation Orien;
    typedef typename mtl::matrix< T, mtl::rectangle<>, 
				  mtl::compressed<int, mtl::external>,
				  Orien >::type LUMatrix;
        enum { Orien_id = Orien::id };
  public:

    typedef typename mtl::IF< EQUAL < Orien_id, ROW_MAJOR >::RET, 
                         preconditioner<LUMatrix, LUMatrix, mtl::lower, 
                                        mtl::unit_upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                         preconditioner<LUMatrix, LUMatrix, mtl::unit_lower,
                                        mtl::upper>,
                         generators_error
    >::RET
    >::RET  Precond;

    typedef typename mtl::IF< EQUAL < Orien_id, ROW_MAJOR >::RET, 
                         preconditioner1<LUMatrix, LUMatrix, mtl::lower, 
                                         mtl::unit_upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                         preconditioner1<LUMatrix, LUMatrix, mtl::unit_lower, 
                                         mtl::upper>,
                         generators_error
    >::RET
    >::RET  Left;

    typedef typename mtl::IF< EQUAL < Orien_id, ROW_MAJOR >::RET, 
                         preconditioner2<LUMatrix, LUMatrix, mtl::lower, 
                                         mtl::unit_upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                         preconditioner2<LUMatrix, LUMatrix, mtl::unit_lower, 
                                         mtl::upper>,
                         generators_error
    >::RET
    >::RET  Right;

  
    SSOR(const Matrix& A) 
      :L_val(A.nnz()),     U_val(A.nnz()),
       L_ind(A.nnz()),     U_ind(A.nnz()),
       L_ptr(A.nrows()+1), U_ptr(A.nrows()+1)
    {
      do_ssor(A, Orien());
    }

  private:
    void
    do_ssor(const Matrix& A, mtl::row_tag)
    {

      using std::upper_bound;
      using std::bind2nd;
      using std::multiplies;
      using std::transform;
      int L_loc=0, U_loc=0;
      L_ptr[0] = 0; 
      U_ptr[0] = 0;
#if 1
      typename Matrix::const_iterator A_i = A.begin();
      if (mtl::not_at(A_i,A.end())) do {
        typename Matrix::OneD A_row = *A_i;
	typedef typename Matrix::OneD::IndexArray::value_type IA_T;
        IA_T i=A_i.index(); 
        typename Matrix::OneD::IndexArray::const_iterator 
          diag_location = 
          upper_bound(A_row.nz_struct().begin(), A_row.nz_struct().end(), i);

        int d = diag_location - A_row.nz_struct().begin();

        std::copy(A_row.begin(), A_row.begin()+d, L_val.begin()+L_loc);
        std::copy(A_row.nz_struct().begin(), A_row.nz_struct().begin()+d, 
             L_ind.begin()+L_loc);

        L_loc += d;
        L_ptr[i+1] = L_loc; 

        T A_ii = *(A_row.begin()+d-1);
        A_ii = T(1)/A_ii;

        transform(A_row.begin()+d, A_row.end(), U_val.begin()+U_loc, 
                  bind2nd(multiplies<T>(), A_ii));

        std::copy(A_row.nz_struct().begin()+d, A_row.nz_struct().end(), 
             U_ind.begin()+U_loc);

        U_loc += A_row.nnz() - d;
        U_ptr[i+1] = U_loc;

        ++A_i;
      } while (mtl::not_at(A_i,A.end()));
#else
      //This is better way in term of general, however, It need to 
      // measure and compare the performance
      typename Matrix::const_iterator A_i = A.begin();
      if (mtl::not_at(A_i,A.end())) do {
	typename Matrix::size_type i = A_i.index();

	T A_ii;
	
        typename Matrix::OneD::const_iterator A_ij = (*A_i).begin();
	if (mtl::not_at(A_ij,(*A_i).end())) do {
	  typename Matrix::size_type j = A_ij.index();
	  if ( j > i ) {
	    U_val[U_loc] = *A_ij;
	    U_ind[U_loc] = j;
	    ++U_loc;
	  } else {
	    L_val[L_loc] = *A_ij;
	    L_ind[L_loc] = j;
	    ++L_loc;
	    if ( i == j ) A_ii = T(1) / *A_ij;
	  }
	  ++A_ij;
	} while (mtl::not_at(A_ij, (*A_i).end()));
	
	U_ptr[i+1] = U_loc;
	L_ptr[i+1] = L_loc;

	//scale upper parts
        transform(U_val.begin()+U_ptr[i], U_val.begin()+U_loc, 
		  U_val.begin()+U_ptr[i], bind2nd(multiplies<T>(), A_ii));

        ++A_i;
      } while (mtl::not_at(A_i,A.end()));
#endif
      L_val.resize(L_loc);
      U_val.resize(U_loc);
      L_ind.resize(L_loc);
      U_ind.resize(U_loc);

      L = LUMatrix(A.nrows(), A.ncols(), 
                   L_loc, &(*L_val.begin()), &(*L_ptr.begin()), 
                   &(*L_ind.begin()));

      U = LUMatrix(A.nrows(), A.ncols(), 
                   U_loc, &(*U_val.begin()), &(*U_ptr.begin()), 
                   &(*U_ind.begin()));
    }

    void
    do_ssor(const Matrix& A, mtl::column_tag)
    {
      using std::upper_bound;
      using std::bind2nd;
      using std::multiplies;
      using std::transform;
      int L_loc=0, U_loc=0;
      L_ptr[0] = 0; 
      U_ptr[0] = 0;

      typename Matrix::const_iterator A_i = A.begin();
      if (mtl::not_at(A_i,A.end())) do {
        const typename Matrix::OneD& A_col = *A_i;

        int i=A_i.index();
        
        typedef typename Matrix::OneD::IndexArray::const_iterator index_iter;
        
        index_iter diag_location = 
          upper_bound(A_col.nz_struct().begin(), A_col.nz_struct().end(), i);

        int d = diag_location - A_col.nz_struct().begin();

        std::copy(A_col.begin(), A_col.begin()+d, U_val.begin()+U_loc);
        std::copy(A_col.nz_struct().begin(), A_col.nz_struct().begin()+d, 
             U_ind.begin()+U_loc);

        U_loc += d;
        U_ptr[i+1] = U_loc; 
        T A_ii = *(A_col.begin()+d-1);
        A_ii = 1./A_ii;

        transform(A_col.begin()+d, A_col.end(), L_val.begin()+L_loc, 
                  bind2nd(multiplies<T>(), A_ii));

        std::copy(A_col.nz_struct().begin()+d, A_col.nz_struct().end(), 
             L_ind.begin()+L_loc);

        L_loc += A_col.nnz() - d;
        L_ptr[i+1] = L_loc;

        ++A_i;
      } while (mtl::not_at(A_i,A.end()));

      L_val.resize(L_loc);
      U_val.resize(U_loc);
      L_ind.resize(L_loc);
      U_ind.resize(U_loc);


      L = LUMatrix(A.nrows(), A.ncols(), 
                   L_loc, &(*L_val.begin()), &(*L_ptr.begin()), 
                   &(*L_ind.begin()));

      U = LUMatrix(A.nrows(), A.ncols(), 
                   U_loc, &(*U_val.begin()), &(*U_ptr.begin()), 
                   &(*U_ind.begin()));
    }
  
  public:
    //: return  a right or Left Preconditioner object.
    Precond operator()() { return Precond(L, U); }
    //: return the Left part of a Split Preconditioner
    Left left() { return Left(L, U); }
    //: return the Right part of a Split Preconditioner
    Right right() { return Right(L, U); }

  private:

    LUMatrix L, U;

    std::vector<T> L_val;
    std::vector<T> U_val;
    std::vector<int> L_ind;
    std::vector<int> U_ind;
    std::vector<int> L_ptr;
    std::vector<int> U_ptr;

  };


}

#endif
