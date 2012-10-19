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

#ifndef ITL_ILU_H
#define ITL_ILU_H


#include <algorithm>
#include <vector>
#include "mtl/entry.h"
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/norm.h"


namespace itl {

  //:  Incomplete LU without fill-in Preconditioner.
  //<codeblock>
  //  Usage:
  //    Matrix A;
  //    ILU<Matrix> precond(A);
  //    qmr(A, x, b, precond.left(), precond.right(), iter);
  //</codeblock>
  // Matrix has to be unsymmetric. 
  //For symmetric one, use incomlete cholesky.
  //Notes: The idea under a concrete Preconditioner such 
  //as Incomplete Cholesky is to create a Preconditioner
  //object to use in iterative methods.
  //
  //!definition: ilu.h
  //!example: ilu.cc
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix -  An MTL Matrix 
  //
  template <class Matrix>
  class ILU  {
    typedef typename Matrix::value_type T;
    typedef typename Matrix::size_type sizeT;
    typedef typename Matrix::orientation Orien;
    typedef typename mtl::matrix< T, mtl::rectangle<>, 
			 mtl::compressed<sizeT, mtl::external>, Orien >::type LUMatrix; 
  public:
    //: The preconditioner type
    typedef  preconditioner<LUMatrix, LUMatrix,
			    mtl::unit_lower, mtl::upper> Precond;
    //: The left preconditioner type
    typedef preconditioner1<LUMatrix, LUMatrix, 
			    mtl::unit_lower, mtl::upper> Left;
    //: The right preconditioner type
    typedef preconditioner2<LUMatrix, LUMatrix,
			    mtl::unit_lower, mtl::upper> Right;

    //: Default Constructor  
    ILU() {}
  
    //: Construct from Matrix A
    ILU(const Matrix& A_) 
      :L_val(A_.nnz()),     U_val(A_.nnz()),
       L_ind(A_.nnz()),     U_ind(A_.nnz()),
       L_ptr(A_.nrows()+1), U_ptr(A_.nrows()+1)
    {
      do_ilu(A_, Orien());
    }
  private:
    void 
    do_ilu(const Matrix& A_, mtl::column_tag)
    {
      using std::upper_bound;
      sizeT L_loc=0, U_loc=0;
      L_ptr[0] = 0; 
      U_ptr[0] = 0;

      typename Matrix::const_iterator A_i = A_.begin();
      if (mtl::not_at(A_i,A_.end())) do {
	typename Matrix::OneD A_row = *A_i;
	sizeT i = A_i.index();
	typename Matrix::OneD::IndexArray::const_iterator diag_location = 
	  std::upper_bound(A_row.nz_struct().begin(), 
			   A_row.nz_struct().end(), i);

	int d = diag_location - A_row.nz_struct().begin();

	std::copy(A_row.begin(), A_row.begin()+d, U_val.begin()+U_loc);
	std::copy(A_row.nz_struct().begin(), A_row.nz_struct().begin()+d, 
	     U_ind.begin()+U_loc);

	U_loc += d;
	U_ptr[i+1] = U_loc; 

	std::copy(A_row.begin()+d, A_row.end(), L_val.begin()+L_loc);
	std::copy(A_row.nz_struct().begin()+d, A_row.nz_struct().end(), 
	     L_ind.begin()+L_loc);

	L_loc += A_row.nnz() - d;
	L_ptr[i+1] = L_loc;

	++A_i;
      } while (mtl::not_at(A_i,A_.end()));

      L_val.resize(L_loc);
      U_val.resize(U_loc);
      L_ind.resize(L_loc);
      U_ind.resize(U_loc);

      sizeT i, j, qn, pn, rn; 
      for (i = 0; i < A_.nrows() - 1; i++) {
	T multiplier = U_val[U_ptr[i+1]-1];
    
	for (j = L_ptr[i]; j < L_ptr[i+1]; j++)
	  L_val[j] /= multiplier;
    
	for (j = U_ptr[i+1]; j < U_ptr[i+2]-1; j++) {
	  multiplier = U_val[j];
	  qn = j + 1;
	  rn = L_ptr[i+1];
	  for (pn = L_ptr[U_ind[j]]; 
	       L_ind[pn] <= i + 1 && pn < L_ptr[U_ind[j]+1]; 
	       pn++) {
	    while (U_ind[qn] < L_ind[pn] && qn < U_ptr[i+2])
	      qn++;
	    if (L_ind[pn] == U_ind[qn] && qn < U_ptr[i+2])
	      U_val[qn] -= multiplier * L_val[pn];
	  }
	  for (; pn < L_ptr[U_ind[j]+1]; pn++) {
	    while (L_ind[rn] < L_ind[pn] && rn < L_ptr[i+2])
	      rn++;
	    if (L_ind[pn] == L_ind[rn] && rn < L_ptr[i+2])
	      L_val[rn] -= multiplier * L_val[pn];
	  }
	}
      }

      L = LUMatrix(A_.nrows(), A_.ncols(), 
		   L_loc, &L_val[0], &L_ptr[0], 
		   &L_ind[0]);
      
      U = LUMatrix(A_.nrows(), A_.ncols(), 
		   U_loc, &U_val[0], &U_ptr[0], 
		   &U_ind[0]);
    }

    void
    do_ilu(const Matrix A_, mtl::row_tag)
    {
      sizeT L_loc=0, U_loc=0;
      L_ptr[0] = 0; 
      U_ptr[0] = 0;

      typename Matrix::const_iterator A_i = A_.begin();
      if (mtl::not_at(A_i,A_.end())) do {
	typename Matrix::OneD A_row = *A_i;
	typedef typename Matrix::OneD::IndexArray::value_type IA_T;
	IA_T i = A_i.index();
	typename Matrix::OneD::IndexArray::const_iterator diag_location = 
	  std::lower_bound(A_row.nz_struct().begin(), 
			   A_row.nz_struct().end(), i);

	int d = diag_location - A_row.nz_struct().begin();

	std::copy(A_row.begin(), A_row.begin()+d, L_val.begin()+L_loc);
	std::copy(A_row.nz_struct().begin(), A_row.nz_struct().begin()+d, 
	     L_ind.begin()+L_loc);

	L_loc += d;
	L_ptr[i+1] = L_loc; 

	std::copy(A_row.begin()+d, A_row.end(), U_val.begin()+U_loc);
	std::copy(A_row.nz_struct().begin()+d, A_row.nz_struct().end(), 
	     U_ind.begin()+U_loc);

	U_loc += A_row.nnz() - d;
	U_ptr[i+1] = U_loc;

	++A_i;
      } while (mtl::not_at(A_i,A_.end()));

      L_val.resize(L_loc);
      U_val.resize(U_loc);
      L_ind.resize(L_loc);
      U_ind.resize(U_loc);

      sizeT i, j, qn, pn, rn; 
      for (i = 1; i < A_.nrows(); i++) {
	for (j = L_ptr[i]; j < L_ptr[i+1]; j++) {
	  pn = U_ptr[L_ind[j]];

	  T multiplier = (L_val[j] /= U_val[pn]);

	  qn = j + 1;
	  rn = U_ptr[i];

	  for (pn++; U_ind[pn] < i && pn < U_ptr[L_ind[j]+1]; pn++) {
	    while (L_ind[qn] < U_ind[pn] && qn < L_ptr[i+1])
	      qn++;
	    if (U_ind[pn] == L_ind[qn] && qn < L_ptr[i+1])
	      L_val[qn] -= multiplier * U_val[pn];
	  }
	  for (; pn < U_ptr[L_ind[j]+1]; pn++) {
	    while (U_ind[rn] < U_ind[pn] && rn < U_ptr[i+1])
	      rn++;
	    if (U_ind[pn] == U_ind[rn] && rn < U_ptr[i+1])
	      U_val[rn] -= multiplier * U_val[pn];
	  }
	}
      }


      L = LUMatrix(A_.nrows(), A_.ncols(), 
		   L_loc, &(*L_val.begin()), &(*L_ptr.begin()), 
		   &(*L_ind.begin()));

      U = LUMatrix(A_.nrows(), A_.ncols(), 
		   U_loc, &(*U_val.begin()), &(*U_ptr.begin()), 
		   &(*U_ind.begin()));

    }

  public:
  

    //:return a left or right Preconditioner object.
    inline Precond operator()()  { return Precond(L, U); }
    //: return the left part of a Split Preconditioner object
    inline Left left() { return Left(L, U); }
    //: return the right part of a Split Preconditioner object
    inline Right right() { return Right(L, U); }
    void print() {
      print_all_matrix(L);
      print_all_matrix(U);
    }

  private:
    LUMatrix L; 
    LUMatrix U;
    std::vector<T> L_val;
    std::vector<T> U_val;
    std::vector<sizeT> L_ind;
    std::vector<sizeT> U_ind;
    std::vector<sizeT> L_ptr;
    std::vector<sizeT> U_ptr;
  };


}

#endif
