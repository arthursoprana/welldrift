//-*-c++-*-
//	=======================================================	
//	CVS Information:                                     	
//	                                                       	
//	$RCSfile: blas_parallel.h,v $  $Revision: 1.4 $  $State: Exp $
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
//	 $Log: blas_parallel.h,v $
//	 Revision 1.4  2001/10/26 16:10:15  llee
//	 *** empty log message ***
//	
//	 Revision 1.3  2001/10/26 15:02:15  llee
//	 *** empty log message ***
//	
//	 Revision 1.2  2001/10/19 15:27:40  llee
//	 parallel bratu examples: done
//	
//	=======================================================	

#ifndef NTL_ITL_BLAS_INTERFACE_H
#define NTL_ITL_BLAS_INTERFACE_H

#include <algorithm>

#include <itl/itl_config.h>
#include <itl/itl_tags.h>
#include <itl/interface/detail/pointer_vector.h>

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

//sun performance library individual function name
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


#ifdef __cplusplus
} // extern "C"
#endif

#define ITL_USE_BLAS
#include <itl/interface/detail/parallel_matrix_v2.h>
#undef  ITL_USE_BLAS

namespace itl {

  template<class Vector>
  struct itl_traits {
    //This type will be used for temp vectors inside of itl algorithm 
    typedef referencing_object_tag vector_category;
    typedef typename Vector::value_type value_type;
    typedef typename Vector::size_type  size_type;
  };

  template <class Vec, class T>
  struct Scaled { 
    
    inline Scaled(const Vec& v, const T& alpha) : _alpha(alpha), _v(v) { }
    
    inline T alpha() const { return _alpha; }
    inline const Vec& vec() const { return _v; }
    inline int size() const { return _v.size(); }

  protected:
    T _alpha;  
    const Vec& _v;
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


  template <class Matrix, class VectorX, class VectorB>
  inline 
  void mult(const parallel_matrix<Matrix>& A, const VectorX& x, VectorB& b) {
    A.process_vector(x);
    A.mult(b);
  }

  template <class Matrix, class VectorX, class VectorY, class VectorB>
  inline 
  void mult(const parallel_matrix<Matrix>& A, const VectorX& x,
	    const VectorY& y, VectorB& b) {
    A.process_vector(x);
    A.mult(b);

    itl::add(y, b);
  }

  template <class VecA, class VecB>
  inline 
  void add(const VecA& v1, VecB& v2)
  {
    int n = v2.size();
    int one = 1;
    double alpha = 1.0;
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(v1.data()),
		   &one, v2.data(), &one);
  }

  template <class VecA, class T, class VecB>
  inline 
  void add(const Scaled<VecA,T>& sv, VecB& v2)
  {
    int n = v2.size();
    int one = 1;
    double alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(sv.vec().data()),
		   &one, v2.data(), &one);
  }

  template <class Vec, class T>
  inline 
  void add(const Scaled<Vec,T>& sv, const Vec& v2, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = sv.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(v2.data()), &one, v3.data(), &one);
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(sv.vec().data()), 
		   &one, v3.data(), &one);
    //v3 = sv.alpha() * sv.vec() + v2;
  }

  template <class Vec, class T>
  inline 
  void add(const Vec& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = sv.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(v2.data()), &one, v3.data(), &one);
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(sv.vec().data()),
		   &one, v3.data(), &one);
    //v3 = sv.alpha() * sv.vec() + v2;
  }

  template <class Vec, class T>
  inline
  void add(const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = v2.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(v2.vec().data()), 
		   &one, v3.data(), &one);
    SP_NAME(dscal)(&n, &alpha, v3.data(), 1);
    alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(sv.vec().data()),
		   &one, v3.data(), &one);
    //v3 = sv.alpha() * sv.vec() + v2.alpha() * v2.vec();
  }
  template <class Vec, class T>
  inline void 
  add(const Vec& g, const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
  {
    int n = v2.size();
    int one = 1;
    double alpha = v2.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(g.data()), &one, v3.data(), &one);
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(v2.vec().data()),
		   &one, v3.data(), &one);
    alpha = sv.alpha();
    SP_NAME(daxpy)(&n, &alpha, const_cast<double*>(sv.vec().data()),
		   &one, v3.data(), &one);
    //v3 = g + sv.alpha() * sv.vec() + v2.alpha() * v2.vec();
  }

  template <class VecA, class VecB>
  inline double dot_conj(const VecA& a, const VecB& b) {
    double ret = itl::dot(a, b);
    return ret;
  }

  template <class VecA, class VecB>
  inline double dot(const VecA& a, const VecB& b) {
    double local, g = 0;
    int n = a.size();
    int one = 1;
    local = SP_NAME(ddot)(&n, const_cast<double*>(a.data()), &one, 
			  const_cast<double*>(b.data()), &one);
    MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, Manager::comm());
    return g;
  }


  template <class Vec>
  inline double 
  two_norm(const Vec& v) {
    double ret = std::sqrt(itl::dot_conj(v, v));
    return ret;
  }
  

  template <class VecA, class VecB>
  inline void copy(const VecA& a, VecB& b) {
    int n = a.size();
    int one = 1;
    SP_NAME(dcopy)(&n, const_cast<double*>(a.data()), &one, b.data(), &one);
  }


  template <class VecA, class T, class VecB>
  inline void copy(const Scaled<VecA, T>& a, VecB& b) {
    int n = b.size();
    int one = 1;
    double alpha = a.alpha();
    SP_NAME(dcopy)(&n, const_cast<double*>(a.vec().data()),
		   &one, b.data(), &one);
    SP_NAME(dscal)(&n, &alpha, b.data(), &one);
   //b = a.vec() * a.alpha();
  }
  
  
  template <class Vec, class T>
  inline void scale(Vec& v, T t) {
    double alpha = t;
    int n = v.size();
    int one = 1;
    SP_NAME(dscal)(&n, &alpha, v.data(), &one);
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
  inline typename Vec::value_type* get_data(const Vec& v) {
    return const_cast<typename Vec::value_type*>(v.data());
  }


  template <class T>
  inline T* get_data(const std::vector<T>& v) {
    return const_cast<T*>(&v[0]);
  }

  template <class T>
  class internal_matrix {
  public:
    typedef T value_type;
    typedef pointer_vector<value_type> OneD;

    inline internal_matrix(int _m, int _n) : m(_m), n(_n) {
      A = new T[m*n];
    }

    inline internal_matrix(const internal_matrix<T>& M) : m(M.m), n(M.n) {
      A = new T[m*n];
      std::copy(M.A, M.A+m*n, A);
    }
    
    inline OneD operator[](int i) { //0 <= i < n
      return OneD(A+i*m, m);
    }

    inline OneD operator[](int i) const { //0 <= i < n
      return OneD(A+i*m, m);
    }

    inline value_type& operator()(int i, int j) {
      return A[i+j*m];
    }

    inline value_type operator()(int i, int j) const {
      return A[i+j*m];
    }
    
    inline value_type* get_val() {
      return A;
    }

    inline const value_type* get_val() const {
      return A;
    }

    inline int nrows() const { return m; }
    inline int ncols() const { return n; }

    inline ~internal_matrix() {
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

}

#endif
