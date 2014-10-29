#ifndef BLOCK_ILU_H
#define BLOCK_ILU_H

#include "itl/preconditioner/detail/preconditioner.h"
#include <vector>

template <class Matrix>
class block_ilu {
  typedef double T;
  typedef typename mtl::matrix< double, mtl::rectangle<>, 
				mtl::compressed<int, mtl::external>, 
				mtl::row_major >::type LUMatrix; 
public:
  
  //: The preconditioner type
  typedef  itl::preconditioner<LUMatrix, LUMatrix,
			       mtl::unit_lower, mtl::upper> Precond;
  //: The left preconditioner type
  typedef itl::preconditioner1<LUMatrix, LUMatrix, 
			       mtl::unit_lower, mtl::upper> Left;
  //: The right preconditioner type
  typedef itl::preconditioner2<LUMatrix, LUMatrix,
			       mtl::unit_lower, mtl::upper> Right;

  block_ilu() {}

  block_ilu(const Matrix& A_, int disp) 
    : L_ptr(A_.nrows()+1),         U_ptr(A_.nrows()+1)
  {
    L_val.reserve(A_.nnz()),     U_val.reserve(A_.nnz());
    L_ind.reserve(A_.nnz()),     U_ind.reserve(A_.nnz());

    _do_ilu(A_, disp, typename Matrix::orientation());
  }

  void operator()(const Matrix& A_, int disp) {
    L_ptr.resize(A_.nrows()+1);
    U_ptr.resize(A_.nrows()+1);

    L_val.reserve(A_.nnz()),     U_val.reserve(A_.nnz());
    L_ind.reserve(A_.nnz()),     U_ind.reserve(A_.nnz());
    _do_ilu(A_, disp, typename Matrix::orientation());
  }

protected:
  void _do_ilu(const Matrix& A_, int disp, mtl::row_tag)
  {
    int block_size = A_.nrows();
    L_ptr[0] = 0;
    U_ptr[0] = 0;
    typename Matrix::const_iterator Ai=A_.begin(), iend = A_.end();
    typename Matrix::OneD::const_iterator Aj, jend;
    for (; Ai != iend; ++Ai) {
      Aj = (*Ai).begin(); jend = (*Ai).end();
      int rowind = Aj.row();
      for (; Aj != jend; ++Aj) {
	int colind = Aj.column() - disp;
	if ( colind < 0 ) continue;
	if ( colind >= block_size ) break;
	
	if ( colind < rowind ) {
	  //L
	  L_ind.push_back(colind);
	  L_val.push_back(*Aj);
	} else {
	  //U
	  U_ind.push_back(colind);
	  U_val.push_back(*Aj);
	}
      }
      U_ptr[rowind+1] = U_ind.size();
      L_ptr[rowind+1] = L_ind.size();
    }
    
    //factorize, assuming column indices in one row are in order
    int i, j, qn, pn, rn; 
    for (i = 1; i < block_size; i++) {
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
    
    L = LUMatrix(block_size, block_size, L_val.size(),
		 &(L_val[0]), &(L_ptr[0]), &(L_ind[0]));
    
    U = LUMatrix(block_size, block_size, U_val.size(),
		 &(U_val[0]), &(U_ptr[0]), &(U_ind[0]));
    
  }
public:
    inline Precond operator()()  { return Precond(L, U); }
    //: return the left part of a Split Preconditioner object
    inline Left left() { return Left(L, U); }
    //: return the right part of a Split Preconditioner object
    inline Right right() { return Right(L, U); }
    void print() {
      print_all_matrix(L);
      print_all_matrix(U);
    }

protected:
  LUMatrix L; 
  LUMatrix U;
  std::vector<T> L_val;
  std::vector<T> U_val;
  std::vector<int> L_ind;
  std::vector<int> U_ind;
  std::vector<int> L_ptr;
  std::vector<int> U_ptr;
};


#endif
