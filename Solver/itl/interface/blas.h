//-*-c++-*-
//	=======================================================	
//	CVS Information:                                     	
//	                                                       	
//	$RCSfile: blas.h,v $  $Revision: 1.7 $  $State: Exp $
//	     $Author: llee $  $Date: 2001/10/26 16:10:15 $ 	
//	     $Locker:  $ 	
//	-------------------------------------------------------	
//	                                                       	
//	 DESCRIPTION                                           	
//
//       This is to provide an itl interface for BLAS with 
//       Matrix as CSR format and Vector as std::vector.
//       In this interface, it require matrix has three 
//       member functions:
//       get_ptr()  //row pointers   in CSR
//       get_ind()  //column indices in CSR
//       get_val()  //values  in CSR
//	                                                       	
//	-------------------------------------------------------	
//	                                                       	
//	 LICENSE AGREEMENT                                     	
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
//	-------------------------------------------------------	
//	                                                       	
//	 REVISION HISTORY:                                     	
//	                                                       	
//	 $Log: blas.h,v $
//	 Revision 1.7  2001/10/26 16:10:15  llee
//	 *** empty log message ***
//	
//	 Revision 1.6  2001/10/26 15:42:54  llee
//	 *** empty log message ***
//	
//	 Revision 1.5  2001/10/26 15:02:15  llee
//	 *** empty log message ***
//	
//	 Revision 1.4  2001/10/19 15:27:40  llee
//	 parallel bratu examples: done
//	
//	 Revision 1.3  2001/10/18 18:55:03  llee
//	 *** empty log message ***
//	
//	 Revision 1.2  2001/10/18 15:44:32  llee
//	 compile and test for examples
//	
//	 Revision 1.1  2001/10/18 14:08:32  llee
//	 re-organize the directory structures
//	
//	 Revision 1.7  2001/07/05 22:28:57  llee1
//	 gcc 3.0 fix
//	
//	 Revision 1.6  2000/07/26 16:28:07  llee1
//	 *** empty log message ***
//	
//	 Revision 1.5  2000/07/26 15:17:05  llee1
//	 *** empty log message ***
//		
//	                                                       	
//	=======================================================	
#ifndef ITL_INTERFACE_BLAS_H
#define ITL_INTERFACE_BLAS_H

#include <algorithm>
#include <vector>

#include <itl/itl_config.h>
#include <itl/itl_tags.h>
#include <itl/interface/detail/pointer_vector.h>
#include <itl/interface/detail/sparse_matrix.h>

#if (defined ITL_NO_STD_ABS) && !(defined MTL_CONFIG_H)
//mtl has included abs if it is not in the compiler header 
namespace std {
  inline double abs(double a) {
    return a > 0 ? a : -a;
  }
  inline float abs(float a) {
    return a > 0 ? a : -a;
  }
  inline int abs(int a) {
    return a > 0 ? a : -a;
  }
}
#endif

//blas library individual function name
#define SP_NAME(X) X##_

#ifdef __cplusplus
extern "C" {
#endif
  
  double SP_NAME(ddot)(int*, double*, int*, double*, int*);
  void   SP_NAME(daxpy)(int*, double*, double*, int*, double*, int*);
  void   SP_NAME(dcopy)(int*, double*, int*, double*, int*);
  double SP_NAME(dnrm2)(int *, double*, int*);
  void   SP_NAME(dscal)(int*, double*, double*, int*);
  
  //triangular solver
  void   SP_NAME(dtrsv)(char* uplo, char* trans, char* diag, int* n, 
			double *da, int* lda, double *dx, int* incx);
  //dense matrix-vector product
  void   SP_NAME(dgemv)(char* trans, int* m, int* n,  double*  alpha,
			double *da,  int* lda, double *dx, int* incx, 
			double* dbeta, double *dy, int* incy);

  void amux(int* n, double* x, double* y, double* a, int* ja, int* ia) {
    for(int i=0; i<*n; ++i) {
      double t = 0;
      int s = ia[i], e = ia[i+1];
      for(int j=s; j<e; ++j)
	t += a[j] * x[ ja[j] ];

      y[i] = t;
    }
  }

  void atmux(int* n, double* x, double* y, double* a, int* ja, int* ia) {
    for(int i=0; i<*n; ++i) y[i] = 0;

    for(int i=0; i<*n; ++i) {
      double t = x[i];
      int s = ia[i], e = ia[i+1];
      for(int j=s; j<e; ++j)
	y[ja[j]] += a[j] * t;
    }
  }
  
#ifdef __cplusplus
} // extern "C"
#endif


namespace itl {
  
  template<class Vector>
  struct itl_traits {
    typedef non_referencing_object_tag vector_category;
    typedef typename Vector::value_type value_type;
    typedef typename Vector::size_type size_type;
  };

  //in VC++ 6 does not support partial specialization yet.
  template<class T>
  struct itl_traits <std::vector<T> > {
    typedef std::vector<T> Vector;
    typedef non_referencing_object_tag vector_category;
    typedef typename Vector::value_type value_type;
    typedef typename Vector::size_type size_type;
  };
  
  template <class Vec, class T>
  struct Scaled { 
    
    inline Scaled(const Vec& v, const T& alpha) : _alpha(alpha), _v(v) { }
    inline T alpha() const { return _alpha; }
    inline const Vec& vec() const { return _v; }
    inline int size() const { return _v.size(); }

  protected:
    T _alpha;  
    const Vec _v;
  };
  
  template <class Vec, class T>
  inline Scaled<Vec, T>
  scaled(const Vec& v, const T& alpha)
  {
    return Scaled<Vec, T> (v, alpha);
  }

  template <class Vec, class T>
  inline Scaled<Vec, T>
  scaled(const Scaled<Vec, T>& v, const T& beta)
  {
    return Scaled<Vec, T> (v.vec(), v.alpha()*beta);
  }

  template <class Matrix, class VecX, class VecY>
  inline void mult(const Matrix& A, const VecX& x, VecY& y) {
    int n = nrows(A);
    amux(&n, const_cast<double*>(get_data(x)), get_data(y), 
	 const_cast<double*>(A.get_val()), 
	 const_cast<int*>(A.get_ind()), 
	 const_cast<int*>(A.get_ptr()));
  }

  template <class Matrix, class VecX, class T, class VecY>
  inline void mult(const Matrix& A, const Scaled<VecX, T>& x, VecY& y) {
    int n = nrows(A);
    //should have Matrix::value_type and Matrix::size_type
    //but amux and blas function only take int*
    amux(&n, const_cast<double*>(get_data(x.vec())), get_data(y), 
	 const_cast<double*>(A.get_val()), 
	 const_cast<int*>(A.get_ind()), 
	 const_cast<int*>(A.get_ptr()));
    itl::scale(y, x.alpha());
  }


  template <class Matrix, class VecX, class VecY, class VecZ>
  inline void mult(const Matrix& A, const VecX& x, const VecY& y, VecZ& z) {
    itl::mult(A, x, z);
    itl::add(y, z);// y has to be not the alias of z
  }

  template <class Matrix, class VecX, class VecY>
  inline void trans_mult(const Matrix& A, const VecX& x, VecY& y) {
    int n = nrows(A);
    //should have Matrix::value_type Matrix::size_type
    atmux(&n, const_cast<double*>(get_data(x)), get_data(y), 
	  const_cast<double*>(A.get_val()), 
	  const_cast<int*>(A.get_ind()), 
	  const_cast<int*>(A.get_ptr()));
  }


  template <class VecA, class VecB>
  inline 
  void add(const VecA& v1, const VecB& v2)
  {
    int n = v2.size();
    int one = 1;
    double alpha = 1.0;
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(v1)),
		   &one, get_data(v2), &one);
  }

  template <class VecA, class T, class VecB>
  inline 
  void add(const Scaled<VecA,T>& sv, const VecB& v2)
  {
    int n = v2.size();
    int one = 1;
    double alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(sv.vec())),
		   &one, get_data(v2), &one);
  }

  template <class Vec, class T>
  inline 
  void add(const Scaled<Vec,T>& sv, const Vec& v2, Vec& v3)
  {
    int n = v3.size();
    int one = 1;
    double alpha = sv.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(get_data(v2)), 
		   &one, get_data(v3), &one);
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(sv.vec())), 
		   &one, get_data(v3), &one);
  }

  template <class Vec, class T>
  inline 
  void add(const Vec& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    itl::add(sv, v2, v3);
  }

  template <class Vec, class T>
  inline
  void add(const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = v2.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(get_data(v2.vec())), 
		   &one, get_data(v3), &one);
    SP_NAME(dscal)(&n, &alpha, get_data(v3), &one);
    alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(sv.vec())),
		   &one, get_data(v3), &one);
  }

  template <class Vec, class T>
  inline void 
  add(const Vec& g, const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = v2.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(get_data(g)), 
		   &one, get_data(v3), &one);
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(v2.vec())),
		   &one, get_data(v3), &one);
    alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(get_data(sv.vec())),
		   &one, get_data(v3), &one);
  }

  //this is not for complex yet
  template <class VecA, class VecB>
  inline double dot_conj(const VecA& a, const VecB& b) {
    double ret = itl::dot(a, b);
    return ret;
  }

  template <class VecA, class VecB>
  inline double dot(const VecA& a, const VecB& b) {
    double local = 0;
    int n = a.size();
    int one = 1;
    local = SP_NAME(ddot)(&n, const_cast<double*>(get_data(a)), &one, 
			  const_cast<double*>(get_data(b)), &one);
    return local;
  }


  template <class Vec>
  inline double 
  two_norm(const Vec& v) {
    double ret = std::sqrt(itl::dot_conj(v, v));
    return ret;
  }
  

  template <class VecA, class VecB>
  inline void copy(const VecA& a, const VecB& b) {
    int n = a.size();
    int one = 1;
    SP_NAME(dcopy)(&n, const_cast<double*>(get_data(a)), &one, 
		   get_data(b), &one);
  }


  template <class VecA, class T, class VecB>
  inline void copy(const Scaled<VecA, T>& a, const VecB& b) {
    int n = b.size();
    int one = 1;
    double alpha = a.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(get_data(a.vec())),
		   &one, get_data(b), &one);
    SP_NAME(dscal)(&n, &alpha, get_data(b), &one);
  }
  
  
  template <class Vec, class T>
  inline void scale(const Vec& v, T t) {
    double alpha = t;
    int n = v.size();
    int one = 1;
    SP_NAME(dscal)(&n, &alpha, get_data(v), &one);
  }

  template <class Vector>
  inline int size(const Vector& x)
  {
    return x.size();
  }

  template <class Vector>
  inline void resize(Vector& x, int sz)
  {
    x.resize(sz);
  }

  template <class Matrix, class VecX>
  inline void upper_tri_solve(const Matrix& A, VecX& x, int i) {
    int n = A.nrows();
    int one = 1;
    char uplo = 'U', trans='N', diag = 'N';
    SP_NAME(dtrsv)(&uplo, &trans, &diag, &i,
		   const_cast<double*>(A.get_val()), &n, get_data(x), &one);
  }

  template <class Vec>
  typename Vec::value_type* get_data(const Vec& v) {
    return const_cast<typename Vec::value_type*>(v.data());
  }


  template <class T>
  T* get_data(const std::vector<T>& v) {
    return const_cast<T*>(&v[0]);
  }

  template <class T>
  class internal_matrix {
  public:
    typedef T value_type;
    typedef pointer_vector<value_type> OneD;

    internal_matrix(int _m, int _n) : m(_m), n(_n) {
      A = new T[m*n];
    }

    internal_matrix(const internal_matrix<T>& M) : m(M.m), n(M.n) {
      A = new T[m*n];
      std::copy(M.A, M.A+m*n, A);
    }
    
    OneD operator[](int i) { //0 <= i < n
      return OneD(A+i*m, m);
    }

    OneD operator[](int i) const { //0 <= i < n
      return OneD(A+i*m, m);
    }

    value_type& operator()(int i, int j) {
      return A[i+j*m];
    }

    value_type operator()(int i, int j) const {
      return A[i+j*m];
    }
    
    value_type* get_val() {
      return A;
    }

    const value_type* get_val() const {
      return A;
    }

    int nrows() const { return m; }
    int ncols() const { return n; }

    ~internal_matrix() {
      if ( A != NULL ) 
	delete [] A;
    }

  protected:
    T* A;
    int m; //number of rows
    int n; //number of columns
  };

  template <class T>
  struct internal_matrix_traits {
    typedef internal_matrix<T> Matrix;
  };

  template<class Matrix>
  int nrows(const Matrix& A) { return A.nrows(); }

}

#include "itl/interface/detail/blas_classical_gram_schmidt.h"

#endif//ITL_INTERFACE_BLAS_H
