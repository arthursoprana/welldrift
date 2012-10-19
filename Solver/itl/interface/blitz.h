//-*-c++-*-
//	=======================================================	
//	CVS Information:                                     	
//	                                                       	
//	$RCSfile: blitz.h,v $  $Revision: 1.5 $  $State: Exp $
//	     $Author: llee $  $Date: 2001/10/26 16:33:55 $ 	
//	     $Locker:  $ 	
//	-------------------------------------------------------	
//	                                                       	
//	 DESCRIPTION                                           	
//
//       This is to provide an itl interface for Blitz++
//       Matrix is matrix-free operator
//       Vector is Blitz++ Array 
//       Only no transpose matrix-vector operation is allowed.
//       Tested Krylov subspace methods are:
//             GMRES
//             BiCGSTAB
//             CGS
//             TFQMR
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
//	 $Log: blitz.h,v $
//	 Revision 1.5  2001/10/26 16:33:55  llee
//	 *** empty log message ***
//	
//	 Revision 1.4  2001/10/26 15:02:15  llee
//	 *** empty log message ***
//	
//	 Revision 1.3  2001/10/23 16:54:53  llee
//	 blitz does not have conj(v)??
//	
//	 Revision 1.2  2001/10/18 21:37:06  llee
//	 *** empty log message ***
//	
//	=======================================================	
#ifndef ITL_INTERFACE_BLITZ_H
#define ITL_INTERFACE_BLITZ_H

#include <itl/itl_config.h>
#include <algorithm>
#include <itl/itl_tags.h>
#include <itl/number_traits.h>
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

namespace std {

  double conj( double& x) 
  { return x; }

  float conj( float& x) 
  { return x; }

}

namespace itl {

  template<class Vector>
  struct itl_traits {
    //This type will be used for temp vectors inside of itl algorithm 
    typedef non_referencing_object_tag vector_category;
    typedef typename Vector::size_type  size_type;
    typedef typename Vector::value_type value_type;
  };


  template <class prec, int dim>
  struct itl_traits<Array<prec, dim> > {

    typedef non_referencing_object_tag vector_category;
    typedef size_t size_type;
    typedef prec value_type;
  };


template <class Vec, class T>
struct Scaled { 

  Scaled(const Vec& v, const T& alpha) : _alpha(alpha), _v(v) { }
  
  T alpha() const { return _alpha; }
  const Vec& vec() const { return _v; }

protected:
  T _alpha;  
  const Vec& _v;
};


template <class prec, int dim, class T>
Scaled<Array<prec, dim>, T>
scaled(const Array<prec, dim>& v, const T& alpha)
{
  return Scaled<Array<prec, dim>, T> (v, alpha);
}

template <class prec, int dim, class T>
Scaled<Array<prec, dim>, T>
scaled(const Scaled<Array<prec, dim>, T>& v, const T& beta)
{
  return Scaled<Array<prec, dim>, T> (v.vec(), v.alpha()*beta);
}

template <class Vec, class T>
void add(const Scaled<Vec,T>& sv, Vec& v2)
{
  v2 += sv.alpha() * sv.vec();
}

template <class Vec>
void add(const Vec& sv, Vec& v2)
{
  v2 += sv;
}

template <class Vec, class T>
void add(const Scaled<Vec,T>& sv, const Vec& v2, Vec& v3)
{
  v3 = sv.alpha() * sv.vec() + v2;
}

template <class Vec, class T>
void add(const Vec& v2, const Scaled<Vec,T>& sv, Vec& v3)
{
  v3 = sv.alpha() * sv.vec() + v2;
}

template <class Vec, class T>
void add(const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
{
  v3 = sv.alpha() * sv.vec() + v2.alpha() * v2.vec();
}
template <class Vec, class T>
void add(const Vec& g, const Scaled<Vec,T>& v2, const Scaled<Vec,T>& sv, Vec& v3)
{
  v3 = g + sv.alpha() * sv.vec() + v2.alpha() * v2.vec();
}

template <class VecA, class VecB>
inline typename VecA::T_numtype dot_conj(const VecA& a, const VecB& b) {
  //return sum(conj(a) * b);
  //blitz does not have conj(a)??
  return sum(a*b);
}

template <class VecA, class VecB>
inline typename VecA::T_numtype dot(const VecA& a, const VecB& b) {
  return sum(a * b);
}


template <class Vec>
inline typename itl::number_traits< typename Vec::T_numtype >::magnitude_type
two_norm(const Vec& v) {
  return std::sqrt(itl::dot_conj(v, v));
}
  

template <class VecA>
inline void copy(const VecA& a, VecA& b) {
  b = a;
  //std::copy(a.begin(), a.end(), b.begin());
}


template <class VecA, class T>
inline void copy(const Scaled<VecA, T>& a, VecA& b) {
  b = a.vec() * a.alpha();
}
  
  
template <class Vec, class T>
inline void scale(Vec& v, T t) {
  v *= t;
}


inline 
Array<double, 2>::T_index
size(const Array<double, 2>& A) 
{
  return A.shape();
}


template <class Size>
inline void
resize(Array<double, 2>& A, const Size& sz) 
{
  A.resize(sz);
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
  typename Matrix::size_type
  nrows(const Matrix& A) { return A.nrows(); }


  template <class Hessenberg, class Vec, class Size>
  void upper_tri_solve(const Hessenberg& hh, Vec& rs, Size i) {
    i--;

    rs[i] /= hh(i, i);
    for (Size ii = 1; ii <= i; ii++) {
      int k  = i - ii;
      int k1 = k + 1;
      double t  = rs[k];
      for (Size j = k1; j <= i; j++) 
	t -= hh(k, j) * rs[j];
      rs[k] = t / hh(k, k);
    }
  }


} //namespace itl

#endif //ITL_INTERFACE_BLITZ_H
