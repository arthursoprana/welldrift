#ifndef ITL_DIAGONAL_PRECOND_H
#define ITL_DIAGONAL_PRECOND_H

namespace itl {

  template <class V>
  struct D_precond {
    D_precond(const V& v) : diag(v) {}
    V diag;
  };


 template <class V, class VectorX, class VectorY>
  inline void solve(const D_precond<V>& M, const VectorX& x, 
		    const VectorY& yy) {
    VectorY& y = const_cast<VectorY&> (yy);
    mtl::ele_mult(M.diag, x, y);
  }

  template <class V, class VectorX, class VectorY>
  inline void trans_solve(const D_precond<V>& M, const VectorX& x, 
		    const VectorY& yy) {
    VectorY& y = const_cast<VectorY&> (yy);
    mtl::ele_mult(M.diag, x, y);
  }

  template <class Matrix>
  struct diagonal_precond {
    typedef typename mtl::matrix_traits<Matrix>::value_type T;
    diagonal_precond(const Matrix& A): diag(A.nrows()) {
      typedef typename mtl::matrix_traits<Matrix>::size_type size_type;
      for (size_type i=0; i<A.nrows(); ++i)
	diag[i] = 1./A(i, i);
    }

    typedef D_precond<mtl::dense1D<T> > Precond;
    inline Precond operator()() { return Precond(diag); }

    mtl::dense1D<T> diag;
  };


}


#endif //ITL_DIAGONAL_PRECOND_H
