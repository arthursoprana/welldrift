#ifndef ITL_APPROXIMATE_INVERSE_PRECOND_H
#define ITL_APPROXIMATE_INVERSE_PRECOND_H

/* Approximate inverse via MR iteration
 * see P301 of Saad book
*/
#include <mtl/matrix.h>

namespace itl {

  /*
    sparse-sparse mode

    row_matrix A * sparse_vector x:
     for ( iter = x.begin(); 

   
  */
  template <class Matrix, class SparseVectorX, class SparseVectorY>
  inline void sparse_mult(const Matrix& A, const SparseVectorX& x,
			  const SparseVectorY& yy) {
    SparseVectorY& y = const_cast<SparseVectorY&>(yy);
    typedef typename mtl::matrix_traits<Matrix>::orientation Orien;
    sparse_mult_matrix_orien(A, x, y, Orien());
  }

  template <class Matrix, class SparseVectorX, class SparseVectorY>
  inline void sparse_mult_matrix_orien(const Matrix& A, const SparseVectorX& x,
				       SparseVectorY& y, mtl::row_tag) {
    assert(false);
  }

  template <class Matrix, class SparseVectorX, class SparseVectorY>
  inline void sparse_mult_matrix_orien(const Matrix& A, const SparseVectorX& x,
				       SparseVectorY& y, mtl::column_tag) {
    typedef typename SparseVectorX::const_iterator Iter;
    for (Iter i = x.begin(); i != x.end(); ++i) {
      typename SparseVectorX::size_type ind = i.index();

      // need a sparse vector add:
      // current mtl does not have a good implementation !!!
      // if A is symmtric shape, only part of column is there: woo!
      SparseVectorY yy(y.size());
      mtl::add(mtl::scaled(A[ind], *i), y, yy);
      y = yy;
    }
  }


  template <class Matrix>
  struct AI_precond {
    AI_precond(const Matrix& _M) : M(_M) {}
    Matrix M;
  };

  template <class Matrix, class VectorX, class VectorY>
  inline void solve(const AI_precond<Matrix>& M, const VectorX& x, 
		    const VectorY& yy) {
    VectorY& y = const_cast<VectorY&> (yy);
    mtl::mult(M.M, x, y);
  }


  template <class T>
  struct dropping_predicate {
    inline dropping_predicate(T _eps) : eps (_eps) {}
    template <class X>
    inline bool operator()(const X& x) {
      if ( x < eps ) 
	return true;
      else
	return false;
    }
    T eps;
  };


  template <class Mat>
  class approximate_inverse {
    typedef typename mtl::matrix_traits<Mat>::value_type T;
    typedef typename mtl::matrix<T, mtl::rectangle<>, 
				 mtl::array< mtl::compressed<> >,
				 mtl::column_major >::type MMatrix;
  public:
    approximate_inverse(const Mat& A, int iter_n, 
			typename mtl::matrix_traits<Mat>::value_type threshold)
      : M(A.nrows(), A.ncols()) {

      typedef typename mtl::matrix_traits<Mat>::size_type size_type;

      T alpha = 0;

      {
	T a = 0, b = 0;
	typename Mat::const_iterator i;  
	typename Mat::OneD::const_iterator j, jend;
	
	for (i = A.begin(); i != A.end(); ++i) {
	  size_type diag = i.index();
	  j = (*i).begin(); jend = (*i).end();
	  for (; j != jend; ++j) {
	    if ( diag == size_type(j.index()) )
	      a += *j;
	    b += *j * *j;
	  }
	}

	alpha = a / b;
      }      

      for (size_type i = 0; i < A.nrows(); ++i) {
	mtl::compressed1D<T> m(A.ncols());
	m.push_back(i, alpha);
	mtl::compressed1D<T> ei(A.ncols());      
	ei.push_back(i, 1);	


	mtl::compressed1D<T> r(A.ncols());

	for (int iter_i = 0; iter_i < iter_n; ++iter_i) {
	  sparse_mult(A, mtl::scaled(m, -1.0), r);
	  
	  mtl::compressed1D<T> rr(A.ncols());
	  mtl::add(ei, r, rr);
	  r = rr;

	  mtl::compressed1D<T> Ar( A.ncols() );
	  sparse_mult(A, r, Ar);

	  T beta = mtl::dot(r, Ar) / mtl::dot(Ar, Ar);

	  mtl::compressed1D<T> mm(A.ncols());
	  mtl::add(mtl::scaled(r, beta), m, mm);
	 
	  //apply numerical dropping to m
#if 0
	  m = mm;
	  typename mtl::compressed1D<T>::iterator pos = std::remove_copy_if(m.begin(), m.end(), m.begin(), dropping_predicate<T>(threshold*mtl::two_norm(m)));
	  m.erase(pos, m.end());
#else
	  m.clear();
	  T act_thres = threshold*mtl::two_norm(mm);
	  for ( typename  mtl::compressed1D<T>::iterator iter = mm.begin();
		iter != mm.end(); ++iter) {
	    if ( std::fabs(*iter) > act_thres )
	      m.push_back(iter.index(), *iter);
	  }
#endif
	}

	mtl::copy(m, M[i]);
      }

    }

    typedef AI_precond<MMatrix> Precond;
    
    inline Precond operator()() 
      { return AI_precond<MMatrix>(M); }

  protected:
    MMatrix M;
  };


}



#endif // ITL_APPROXIMATE_INVERSE_PRECOND_H
