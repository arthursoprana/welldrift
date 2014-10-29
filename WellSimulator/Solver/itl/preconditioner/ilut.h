
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

#ifndef ITL_ILUT_H
#define ITL_ILUT_H

#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>
#include "mtl/entry.h"
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/norm.h"
#include "mtl/meta_if.h"
#include "mtl/meta_equal.h"

namespace itl {

  using mtl::IF; using mtl::ROW_MAJOR; using mtl::COL_MAJOR;
  using mtl::lower; using mtl::upper; 
  using mtl::unit_lower; using mtl::unit_upper;
  using mtl::generators_error;

  struct entry1_value_less {
    template <class real>
    inline bool operator()(const mtl::entry1<real>& a, 
			   const mtl::entry1<real>& b) const 
    {
      return ( MTL_ABS(a.value) < MTL_ABS(b.value));
    }
  };
  

  /*
    Performane comparing for SSOR, ILU and ILUT based on sherman 5 matrix 
    in Harwell-Boeing collection on Sun Ultra 30 UPA/PCI (UltraSPARC-II 296MHz)
    Preconditioner & Factorization time  &  Number of Iteration \\ \hline
    SSOR        &   0.010577  & 41 \\
    ILU         &   0.019336  & 32 \\
    ILUT with 0 fill-in and threshold of 1.0e-6 & 0.343612 &  23 \\
    ILUT with 5 fill-in and threshold of 1.0e-6 & 0.343612 &  18 \\ \hline
  */

  //: ILUT:  Incomplete LU with threshold and K fill-in Preconditioner.
  //  The algorithm of ILUT(A, 0, 1.0e-6) is slower than ILU(A). If No fill-in 
  //  is arrowed , please use ILU instead of ILUT.
  //
  // <codeblock>
  //  Usage:
  //    Matrix A;
  //    int fill_in=5;
  //    double threshold = 1.0e-5;
  //    ILUT<Matrix> precond(A, fill_in, threshold);
  //    qmr(A, x, b, precond.left(), precond.right(), iter);
  // </codeblock>
  //
  // Matrix has to be unsymmetric. 
  // For symmetric one, use incomlete cholesky.
  // Notes: The idea under a concrete Preconditioner such 
  // as Incomplete LU is to create a Preconditioner
  // object to use in iterative methods. 
  //
  //!definition: ilut.h
  //!example: ilut.cc
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix - Matrix 
  //
  template <class Matrix>
  class ILUT  {
    typedef typename Matrix::orientation Orien;
    typedef typename Matrix::value_type T;
        enum { Orien_id = Orien::id };
        typedef std::vector< mtl::entry1<T> > entry_vec;
        typedef std::vector<T> T_vec;
  public:
    //: The preconditioner type    
    typedef typename mtl::IF < EQUAL < Orien_id, ROW_MAJOR >::RET, 
                          preconditioner<Matrix, Matrix, unit_lower, upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                          preconditioner<Matrix, Matrix, lower, unit_upper>,
                          generators_error
    >::RET
    >::RET Precond;
    //: The left preconditioner type
    typedef typename mtl::IF < EQUAL < Orien_id, ROW_MAJOR >::RET,
                          preconditioner1<Matrix, Matrix, unit_lower, upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                          preconditioner1<Matrix, Matrix, lower, unit_upper>,
                          generators_error
    >::RET
    >::RET Left;

    //: The right preconditioner type
    typedef typename mtl::IF < EQUAL < Orien_id, ROW_MAJOR >::RET,
                          preconditioner2<Matrix, Matrix, unit_lower, upper>,
            typename mtl::IF< EQUAL < Orien_id, COL_MAJOR >::RET, 
                          preconditioner2<Matrix, Matrix, lower, unit_upper>,
                          generators_error
    >::RET
    >::RET Right;

    //: Default Constructor
    ILUT() : K(0), eps(0.), dropped(0) {}

    //: Construct from Matrix A, max fill-ins (k), and threshold (eps)
    ILUT(const Matrix& A, int k_, double eps_) 
      : K(k_), eps(eps_), dropped(0), 
      L(A.nrows(), A.ncols()), U(A.nrows(), A.ncols())
    {
      do_ilut(A, Orien());
    }

    void print() {
      print_all_matrix(L);
      print_all_matrix(U);
    }
  private:
    void 
    do_ilut(const Matrix& A, mtl::row_tag)
    {
      typedef Matrix RowMatrix;

      T_vec indiag(A.nrows());
      
      entry_vec* w = new entry_vec;
      entry_vec* wswap = new entry_vec;

      typename RowMatrix::const_iterator A_i = A.begin();

      if (mtl::not_at(A_i,A.end())) do {
        const typename RowMatrix::OneD& A_row = *A_i;

        typename RowMatrix::OneD::const_iterator A_ij = A_row.begin();
        typename RowMatrix::OneD::const_iterator row_end = A_row.end();

        int i=A_i.index();

        int ninrow = A_row.nnz();
        w->resize(ninrow);
        double norm_row = 0.;
        int nL = 0, nU = 0;
        int inrow = 0;
        if (mtl::not_at(A_ij, row_end)) do {

          (*w)[inrow].value = *A_ij;
          (*w)[inrow].index = A_ij.index();
          inrow++;
          norm_row += std::norm(T(*A_ij));
          if ( i > A_ij.index() ) nL++;
          
          ++A_ij;
        } while (mtl::not_at(A_ij, row_end));
 
        norm_row = sqrt(norm_row);

        nU = ninrow - nL - 1;

        norm_row = 1./norm_row;

        typename entry_vec::iterator wk= w->begin();
        typename entry_vec::iterator wkend = w->end();
        entry_vec tmp_w;
        int krow = 0;
        while( wk != wkend ) {
          int k = (*wk).index;
          if ( k >= i ) break;
          T tmp = (*wk).value;
          tmp = tmp * indiag[k];  
          if ( std::fabs(tmp) * norm_row < eps ) {
            w->erase(wk); 
            dropped++;
          } else {
            (*wk).value = tmp;
            
            typename RowMatrix::RowRef U_k = U[k];
            typename RowMatrix::RowRef::iterator U_kj = U_k.begin();
            ++U_kj;
            typename RowMatrix::RowRef::iterator U_kend = U_k.end();
            tmp_w.resize( U_k.nnz() - 1);

            int irow = 0;
            if (mtl::not_at(U_kj, U_kend)) do {
              tmp_w[irow].value = -tmp*( *U_kj );
              tmp_w[irow].index = U_kj.index();
              irow++;
              ++U_kj;
            } while (mtl::not_at(U_kj, U_kend));
            
            wswap->resize(w->size() + tmp_w.size());
            typename entry_vec::iterator wj = w->begin();
            typename entry_vec::iterator tmp_wj = tmp_w.begin();
            int j = 0;
            while ( 1 ) {
              if ( wj == w->end() ) {

                if ( mtl::not_at(tmp_wj, tmp_w.end()) ) do {
                  (*wswap)[j] = *tmp_wj; ++tmp_wj; ++j;
                } while ( mtl::not_at(tmp_wj, tmp_w.end()) );
                
                break;
              }
              
              if ( tmp_wj == tmp_w.end() ) {
                if ( mtl::not_at(wj, w->end()) ) do {
                  (*wswap)[j] = *wj; ++wj; ++j;
                } while (mtl::not_at(wj, w->end()));
                break;
              }
              
              if ( *wj == *tmp_wj ) {
                (*wswap)[j] = *wj; 
                (*wswap)[j].value += (*tmp_wj).value;
                ++tmp_wj; ++wj; ++j;
                continue;
              }

              if ( *wj < *tmp_wj ) {
                (*wswap)[j] = *wj; 
                ++wj; ++j;
                continue;
              }

              (*wswap)[j] = *tmp_wj; 
              ++tmp_wj;  ++j;

            }

            wswap->resize(j);
            entry_vec* tswap = wswap;
            wswap = w;
            w = tswap;

            krow++; 
          }
          wk = w->begin()+krow;
          wkend = w->end();
        }

        if (i) {
          typename entry_vec::iterator wi = w->begin();
          int jj = 0;
          while( wi != w->end() ) {
            if ( (*wi).index != i)
              if ( std::fabs((*wi).value)*norm_row < eps ) {
                w->erase(wi);
                dropped++;
                jj--;
              }
            jj++;
            wi = w->begin() + jj;
          }

          typename entry_vec::iterator diag = 
            std::find(w->begin(), w->end(), mtl::entry1<T>(i));
          int m = diag-w->begin();
          std::make_heap(w->begin(), diag, entry1_value_less());
          int jmax = MTL_MIN(nL+K, m);
          if (jmax != m) dropped += m-jmax;
          for ( int j=m; j>m-jmax; j-- ) {   
            typename entry_vec::iterator first = w->begin();
          
            L(i, (*first).index) = (*first).value;

            std::pop_heap(first, first+j, entry1_value_less());
          }

          U(i, i) = (*diag).value;

          m = w->end() - diag - 1;
          jmax = MTL_MIN(nU+K, m);
          if ( m > 0 ) {
            if (jmax != m) dropped += m-jmax;
            std::make_heap(diag+1, w->end(), entry1_value_less());

            for ( int j=m; j>m-jmax; j-- ) {   
              typename entry_vec::iterator first = diag+1;
            
              U(i, (*first).index) = (*first).value;

              std::pop_heap(first, first+j, entry1_value_less());
            }
          }
        } else { 
          for (typename entry_vec::iterator wj= w->begin();
               wj != w->end(); ++wj)
            U(i, (*wj).index) = (*wj).value;
        }

        indiag[i] = 1.0 / T( *(U[i].begin()) ); 

        ++A_i;
      } while (mtl::not_at(A_i,A.end()));

      delete w;
      delete wswap;
    }


    void 
    do_ilut(const Matrix A_, mtl::column_tag)
    {
      typedef typename Matrix::transpose_type  RowMatrix;

      const RowMatrix A=trans(A_);
      
      RowMatrix LA = trans(U);
      RowMatrix UA = trans(L);

      T_vec indiag(A.nrows()); 
      
      entry_vec* w = new entry_vec;
      entry_vec* wswap = new entry_vec;

      typename RowMatrix::const_iterator A_i = A.begin();

      if (mtl::not_at(A_i,A.end())) do {
        const typename RowMatrix::OneD& A_row = *A_i;

        typename RowMatrix::OneD::const_iterator A_ij = A_row.begin();
        typename RowMatrix::OneD::const_iterator row_end = A_row.end();

        int i=A_i.index();

        int ninrow = A_row.nnz();
        w->resize(ninrow);
        double norm_row = 0.;
        int nL = 0, nU = 0;
        int inrow = 0;
        if (mtl::not_at(A_ij, row_end)) do {

          (*w)[inrow].value = *A_ij;
          (*w)[inrow].index = A_ij.index();
          inrow++;
          norm_row += std::norm(T(*A_ij));
          if ( i > A_ij.index() ) nL++;
          
          ++A_ij;
        } while (mtl::not_at(A_ij, row_end));
 
        norm_row = sqrt(norm_row);

        nU = ninrow - nL - 1;

        norm_row = 1./norm_row;

        typename entry_vec::iterator wk= w->begin();
        typename entry_vec::iterator wkend = w->end();
        entry_vec tmp_w;
        int krow = 0;
        while( wk != wkend ) {
          int k = (*wk).index;
          if ( k >= i ) break;
          T tmp = (*wk).value;
          tmp = tmp * indiag[k];  
          if ( std::fabs(tmp) * norm_row < eps ) {
            w->erase(wk); 
            dropped++;
          } else {
            (*wk).value = tmp;
            
            typename RowMatrix::Row U_k = UA[k];
            typename RowMatrix::Row::iterator U_kj = U_k.begin();
            ++U_kj;
            typename RowMatrix::Row::iterator U_kend = U_k.end();
            tmp_w.resize( U_k.nnz() - 1);

            int irow = 0;
            if (mtl::not_at(U_kj, U_kend)) do {
              tmp_w[irow].value = -tmp*( *U_kj );
              tmp_w[irow].index = U_kj.index();
              irow++;
              ++U_kj;
            } while (mtl::not_at(U_kj, U_kend));
            
            wswap->resize(w->size() + tmp_w.size());
            typename entry_vec::iterator wj = w->begin();
            typename entry_vec::iterator tmp_wj = tmp_w.begin();
            int j = 0;
            while ( 1 ) {
              if ( wj == w->end() ) {

                if ( mtl::not_at(tmp_wj, tmp_w.end()) ) do {
                  (*wswap)[j] = *tmp_wj; ++tmp_wj; ++j;
                } while ( mtl::not_at(tmp_wj, tmp_w.end()) );
                
                break;
              }
              
              if ( tmp_wj == tmp_w.end() ) {
                if ( mtl::not_at(wj, w->end()) ) do {
                  (*wswap)[j] = *wj; ++wj; ++j;
                } while (mtl::not_at(wj, w->end()));
                break;
              }
              
              if ( *wj == *tmp_wj ) {
                (*wswap)[j] = *wj; 
                (*wswap)[j].value += (*tmp_wj).value;
                ++tmp_wj; ++wj; ++j;
                continue;
              }

              if ( *wj < *tmp_wj ) {
                (*wswap)[j] = *wj; 
                ++wj; ++j;
                continue;
              }

              (*wswap)[j] = *tmp_wj; 
              ++tmp_wj;  ++j;

            }

            wswap->resize(j);
            entry_vec* tswap = wswap;
            wswap = w;
            w = tswap;

            krow++; 
          }
          wk = w->begin()+krow;
          wkend = w->end();
        }

        if (i) {
          typename entry_vec::iterator wi = w->begin();
          int jj = 0;
          while (wi != w->end()) {
            if ( (*wi).index != i)
              if ( MTL_ABS((*wi).value)*norm_row < eps ) {
            w->erase(wi);
            dropped++;
            jj--;
              }
        jj++;
        wi = w->begin() + jj;
          }

          typename entry_vec::iterator diag = 
            std::find(w->begin(), w->end(), mtl::entry1<T>(i));
          int m = diag - w->begin();
          std::make_heap(w->begin(), diag, entry1_value_less());
          int jmax = MTL_MIN(nL+K, m);
          if (jmax != m) dropped += m-jmax;
          for (int j = m; j > m - jmax; j-- ) {  
            typename entry_vec::iterator first = w->begin();
          
            LA(i, (*first).index) = (*first).value;

            std::pop_heap(first, first+j, entry1_value_less());
          }

          UA(i, i) = (*diag).value;

          m = w->end() - diag - 1;
          jmax = MTL_MIN(nU+K, m);
          if ( m > 0 ) {
            if (jmax != m) dropped += m-jmax;
            std::make_heap(diag+1, w->end(), entry1_value_less());

            for ( int j=m; j>m-jmax; j-- ) {
              typename entry_vec::iterator first = diag+1;
            
              UA(i, (*first).index) = (*first).value;

              std::pop_heap(first, first+j, entry1_value_less());
            }
          }
        } else {
          for (typename entry_vec::iterator wj= w->begin();
               wj != w->end(); ++wj)
            UA(i, (*wj).index) = (*wj).value;
        }

        indiag[i] = 1.0 / T( *(UA[i].begin()) ); 

        ++A_i;
      } while (mtl::not_at(A_i,A.end()));

      delete w;
      delete wswap;
    }

  public:
    //: return  a right or Left Preconditioner object
    inline Precond operator()() { return Precond(L, U); }
    //: return the Left part of a Split Preconditioner object
    inline Left left() { return Left(L, U); }
    //: return the Right part of a Split Precondtioner objet
    inline Right right() { return Right(L, U); }
  private:
    int K;
    double eps;
    int dropped;
    Matrix L; 
    Matrix U;

  };



}

#endif
