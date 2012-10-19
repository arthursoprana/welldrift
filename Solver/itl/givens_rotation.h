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
#ifndef ITL_GIVENS_ROTATION_H
#define ITL_GIVENS_ROTATION_H
#include <math.h>

#define ITL_BLAS_GROT 1
/*
  If I use lapack style in givens_rotation,
  For the following matrix, gmres with ssor does not converge well (It
  should converge at fifth step since it has complete space at that
  step. But the outer iteration of gmres reported that the residual is 
  relatively large.)

  With ITL_BLAS_GROT defined (old blas version), gmres with ssor
  converges very well. That is the reason why I defined ITL_BLAS_GROT above.

  Results of gmres with ssor iterations are:

     1      2      0      0      3   x          0.205847  =       1
     4      5      0      6      0   x         0.0419363  =       1
     0      7      8      0      9   x         -0.178049  =       1
     0      0     10     11     12   x       -0.00551162  =       1
     0      0     13      0     14   x           0.23676  =       1


iteration 0: resid 4.27219
iteration 1: resid 0.510264
iteration 2: resid 0.0372096
iteration 3: resid 0.00883377
iteration 4: resid 1.78958e-15
iteration 4: resid 2.15058e-15
finished! error code = 0
4 iterations
2.15058e-15 is actual final residual.
5.0339e-16 is actual relative tolerance achieved.
Relative tol: 1e-06  Absolute tol: 0
*/


#include <complex>

namespace itl {

template <class T>
class givens_rotation {
public:

  //: Default constructor
  inline givens_rotation() 
    : 
#ifdef ITL_BLAS_GROT
    a_(0), b_(0),
#endif
    c_(0), s_(0)
#ifndef ITL_BLAS_GROT
    , r_(0) 
#endif
  { }

  //: Givens Plane Rotation Constructor
  inline givens_rotation(const T& a_in, const T& b_in) {
#ifdef ITL_BLAS_GROT // old BLAS version
    T roe;
    if (std::fabs(a_in) > std::fabs(b_in))
      roe = a_in;
    else
      roe = b_in;
    
    T scal = std::fabs(a_in) + std::fabs(b_in);
    T r, z;
    if (scal != T(0)) {
      T a_scl = a_in / scal;
      T b_scl = b_in / scal;
      r = scal * std::sqrt(a_scl * a_scl + b_scl * b_scl);
      if (roe < T(0)) r *= -1;
      c_ = a_in / r;
      s_ = b_in / r;
      z = 1;
      if (std::fabs(a_in) > std::fabs(b_in))
        z = s_;
      else if (std::fabs(b_in) >= std::fabs(a_in) && c_ != T(0))
        z = T(1) / c_;
    } else {
      c_ = 1; s_ = 0; r = 0; z = 0;      
    }
    a_ = r;
    b_ = z;
#else // similar LAPACK slartg version, modified to the NEW BLAS proposal
#if 1
    /*
     *     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
     *     [ -SN  CS  ]     [ G ]     [ 0 ]
     *
     *  This is a slower, more accurate version of the BLAS1 routine DROTG,
     *  with the following other differences:
     *     F and G are unchanged on return.
     *     If G=0, then CS=1 and SN=0.
     *     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
     *        floating point operations (saves work in DBDSQR when
     *        there are zeros on the diagonal).
     *
     *  If F exceeds G in magnitude, CS will be positive.
     */
    if (b_in == T(0)) {
      c_ = T(1);
      s_ = T(0);
      r_ = a_in;
    } else if (a_in == T(0)) {
      c_ = T(0);
      s_ = T(1);
      r_ = b_in;
    } else {
      T a = a_in, b = b_in;
      T eps = std::numeric_limits<T>::epsilon();
      T safemin = std::numeric_limits<T>::min();
      T safemn2 = std::pow(2., std::log(safemin/eps)/std::log(2.0)/2.);
      T safemx2 = 1./safemn2;

      T scale = std::max(std::fabs(a), std::fabs(b));

      int count = 0;
      if ( scale >= safemx2 ) {
	while ( scale >= safemx2 ) {
	  ++count;
	  a *= safemn2;
	  b *= safemn2;
	  scale = std::max(std::fabs(a), std::fabs(b));
	}
	
	r_ = std::sqrt( a * a + b * b);
	c_ = a / r_;
	s_ = b / r_;
	
	for (int i=0; i<count; ++i)
	  r_ *= safemx2;
      } else if ( scale <= safemn2 ) {
	while ( scale <= safemn2 ) {
	  ++count;
	  a *= safemx2;
	  b *= safemx2;
	  scale = std::max(std::fabs(a), std::fabs(b));
	}
	
	r_ = std::sqrt( a * a + b * b);
	c_ = a / r_;
	s_ = b / r_;

	for (int i=0; i<count; ++i)
	  r_ *= safemn2;
      } else {
	r_ = std::sqrt( a * a + b * b);
	c_ = a / r_;
	s_ = b / r_;
      }

      if ( (std::fabs(a_in) > std::fabs(b_in)) && c_ < T(0) ) {
	c_ = -c_;
	s_ = -s_;
	r_ = -r_;
      }
    }
#else
    //who knows where the following code came from!
    T a = a_in, b = b_in;
    if (b == T(0)) {
      c_ = T(1);
      s_ = T(0);
      r_ = a;
    } else if (a == T(0)) {
      c_ = T(0);
      s_ = sign(b);
      r_ = b;
    } else {
      // cs = |a| / std::sqrt(|a|^2 + |b|^2)
      // sn = sign(a) * b / std::sqrt(|a|^2 + |b|^2)
      T abs_a = std::fabs(a);
      T abs_b = std::fabs(b);
      if (abs_a > abs_b) {
        // 1/cs = std::sqrt( 1 + |b|^2 / |a|^2 )
        T t = abs_b / abs_a;
        T tt = std::sqrt(T(1) + t * t);
        c_ = T(1) / tt;
        s_ = t * c_;
        r_ = a * tt;
      } else {
        // 1/sn = sign(a) * std::sqrt( 1 + |a|^2/|b|^2 )
        T t = abs_a / abs_b;
        T tt = std::sqrt(T(1) + t * t);
        s_ = sign(a) / tt;
        c_ = t * s_;
        r_ = b * tt;
      }
    }
#endif 

#endif
  }

  inline void set_cs(const T& cin, const T& sin) { c_ = cin; s_ = sin; }

  //: Apply plane rotation to two real scalars. (name change a VC++ workaround)
  inline void scalar_apply(T& x, T& y) {
    T tmp = c_ * x + s_ * y;
    y = c_ * y - s_ * x;
    x = tmp;
  }

  //: Apply plane rotation to two vectors.
  template <class VecX, class VecY>
  inline void apply(const VecX& x_, const VecY& y_) {
    VecX& x = const_cast<VecX&>(x_);
    VecY& y = const_cast<VecY&>(y_);

    typename VecX::iterator xi = x.begin();
    typename VecX::iterator xend = x.end();
    typename VecY::iterator yi = y.begin();

    //while (mtl::not_at(xi, xend)) {
    while ( xi != xend ) {
      scalar_apply(*xi, *yi);
      ++xi; ++yi;
    }
  }

#ifdef ITL_BLAS_GROT
  inline T a() { return a_; }
  inline T b() { return b_; }
#endif
  inline T c() { return c_; }
  inline T s() { return s_; }
#ifndef ITL_BLAS_GROT
  inline T r() { return r_; }
#endif
protected:
  T sign(const T& t) { T ret = 1; if ( t < T(0) ) ret = -1; return ret; }
#ifdef ITL_BLAS_GROT
  T a_, b_;
#endif
  T c_, s_;
#ifndef ITL_BLAS_GROT
  T r_;
#endif
};



//:  The specialization for complex numbers.

  /*
  Using a and b to represent elements of an input complex vector, the CROTG
  and ZROTG functions calculate the elements real c and complex s of an
  orthogonal matrix such that:

              c*a + s*b = r
  -conjugate(s)*a + c*b = 0
  */


//!category: functors
//!component: type
template <class T>
class givens_rotation < std::complex<T> > {
  typedef std::complex<T> C;
public:
  //:
  inline givens_rotation() : cs(0), sn(0)
#ifndef ITL_BLAS_GROT
    , r_(0)
#endif
  { }
  
  inline T abs_sq(C t) 
  { return std::real(t) * std::real(t) + std::imag(t) * std::imag(t); }
#ifndef ITL_BLAS_GROT
  inline T abs1(C t) 
  { return std::fabs(std::real(t)) + std::fabs(std::imag(t)); }
  inline T abs2(C t) 
  { return std::max(std::fabs(std::real(t)), std::fabs(std::imag(t))); }
  
  inline T scaled_abs(C t) {
    T ret;
    T x = std::fabs(std::real(t));
    T y = stdLLabs(std::imag(t));
    T w = std::max(x, y);
    T z = std::min(x, y);
    if ( z == T(0) )
      ret = w;
    else 
      ret = w * std::sqrt( T(1.0) + (z/w) * (z/w) );
    return ret;
  }
#endif 
    
  //:
  inline givens_rotation(const C& a_in, const C& b_in) {
#ifdef ITL_BLAS_GROT

    T a = std::fabs(a_in), b = std::fabs(b_in);
    if ( a == T(0) ) {
      cs = T(0);
      sn = C(1.);
    } else {
      T scale = a + b;
      T norm = std::sqrt(abs_sq(a_in/scale)+abs_sq(b_in/scale)) * scale;
    
      cs = a / norm;
      sn = a_in/a * std::conj(b_in)/norm;
      
      //in zrotg there is an assignment for ca, what is that for? 
    }
#else // LAPACK version, clartg
    /* lapack comments in clartg
     *  CLARTG generates a plane rotation so that
     *
     *     [  CS  SN  ]     [ F ]     [ R ]
     *     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
     *     [ -SN  CS  ]     [ G ]     [ 0 ]
     *
     *  This is a faster version of the BLAS1 routine CROTG, except for
     *  the following differences:
     *     F and G are unchanged on return.
     *     If G=0, then CS=1 and SN=0.
     *     If F=0, then CS=0 and SN is chosen so that R is real.
     *
     *  Arguments
     *  =========
     *
     *  F       (input) COMPLEX
     *          The first component of vector to be rotated.
     *
     *  G       (input) COMPLEX
     *          The second component of vector to be rotated.
     *
     *  CS      (output) REAL
     *          The cosine of the rotation.
     *
     *  SN      (output) COMPLEX
     *          The sine of the rotation.
     *
     *  R       (output) COMPLEX
     *          The nonzero component of the rotated vector.
     *
     *  Further Details
     *  ======= =======
     *
     *  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
     *
     *  =====================================================================
     */
#  if 1
    T eps = std::numeric_limits<T>::epsilon();
    T safemin = std::numeric_limits<T>::min();
    T safemn2 = std::pow(2., std::log(safemin/eps)/std::log(2.0)/2.);
    T safemx2 = 1./safemn2;
    
    T scale  = std::max(abs2(a_in), abs2(b_in));
    C f = a_in, g = b_in;
    int count = 0;
    if ( scale >= safemx2 ) {
      while ( scale >= safemx2 ) {
	++count;
	f *= safemn2;
	g *= safemn2;
	scale *= safemn2;
      }
    } else if (scale <= safmn2 ) {
      if ( g == C(0) ) {
	cs = T(1);
	sn = C(0);
	r_ = f;
	return ; ///??????
      }

      while (scale <= safmn2 ) {
	--count;
	f *= safemx2;
	g *= safemx2;
	scale *= safemx2;
      }
    }
    
    T f2 = abs_sq(f), g2 = abs_sq(g);

    if ( f2 < std::max(g2, 1.0)*safemin ) {
      //rare case a_in is very small
      if ( f == C(0) ) {
	cs = C(0);
	r_ = scaled_abs(b_in);
	sn = std::conj(g)/scaled_abs( g );
	return ; /////
      }

      T f2s = scaled_abs(f);
      T g2s = std::sqrt(g2);
      cs = f2s/g2s;

      T ff;
      if ( abs2(a_in) > T(1.) )
	ff = a_in/scaled_abs(a_in);
      else {
	C d = a_in * safemx2;
	ff = d / scaled(d);
      }
      sn = ff * std::conj(g)/g2s;
      r_ = cs * a_in + sn * b_in;
    } else {
      //common case
      T f2s = std::sqrt(T(1.0)+g2/f2);
      r_ = f2s * f;
      cs = T(1.)/f2s;
      sn = r/(f2+g2);
      sn *= std::conj(g);

      if ( count > 0 ) {
	for (int i=0; i<count; ++i)
	  r *= safmx2;
      } else if ( count < 0 ) {
	for (int i=0; i<count; ++i)
	  r *= safmn2;
      }
    }

#  else
    C f(a_in), g(b_in);
    if (g == C(0)) {
      cs = T(1);
      sn = C(0);
      r_ = f;
    } else if (f == C(0)) {
      cs = T(0);
      sn = std::conj(g) / std::fabs(g);
      r_ = std::fabs(g);
    } else {
      C fs, gs, ss;
      T d, di, f1, f2, fa, g1, g2, ga;
      f1 = abs1(f);
      g1 = abs1(g);
      if (f1 >= g1) {
        gs = g / f1;
        g2 = abs_sq(gs);
        fs = f / f1;
        f2 = abs_sq(fs);
        d = std::sqrt(T(1) + g2 / f2);
        cs = T(1) / d;
        sn = std::conj(gs) * fs * (cs / f2);
        r_ = f * d;
      } else {
        fs = f / g1;
        f2 = abs_sq(fs);
        fa = std::sqrt(f2);
        gs = g / g1;
        g2 = abs_sq(gs);
        ga = std::sqrt(g2);

        d = std::sqrt(T(1) + f2 / g2);
        di = T(1) / d;
        cs = (fa / ga ) * di;
        ss = (std::conj(gs) * fs) / (fa * ga);
        sn = ss * di;
        r_ = g * ss * d;
      }
    }
#  endif

#endif
  }
  //:  Apply plane rotation to two vectors.
  template <class VecX, class VecY>
  inline void apply(const VecX& x_, const VecY& y_) {
    VecX& x = const_cast<VecX&>(x_);
    VecY& y = const_cast<VecY&>(y_);
    
    typename VecX::iterator xi = x.begin();
    typename VecX::iterator xend = x.end();
    typename VecY::iterator yi = y.begin();
    
    //while (mtl::not_at(xi, xend)) {
    while (xi != xend ) {
      scalar_apply(*xi, *yi);
      ++xi; ++yi;
    }
  }
  //: Apply plane rotation to two complex scalars.
  inline void scalar_apply(C& x, C& y) {
    //complex<T> temp  =  std::conj(cs) * x + std::conj(sn) * y;
    //y = cs * y - sn * x; 
    std::complex<T> temp  =  cs * x + sn * y;
    y = cs * y - std::conj(sn) * x;     
    x = temp;
  }
  inline void set_cs(const T& cs_, const C& sn_) {
    cs = cs_; sn = sn_;
  }

  inline T c() { return cs; }
  inline C s() { return sn; }
#ifndef ITL_BLAS_GROT
  inline C r() { return r_; }
#endif

protected:
  T cs;
  C sn;
#ifndef ITL_BLAS_GROT
  C r_;
#endif
};
}


#undef ITL_BLAS_GROT //do not escape it out of the file scope

#endif
