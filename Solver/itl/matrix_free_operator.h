/* -*- c++ -*- */
// 	$Id: matrix_free_operator.h,v 1.5 2001/10/26 14:28:25 llee Exp $	
#ifndef MATRIX_FREE_OPERATOR_H
#define MATRIX_FREE_OPERATOR_H

namespace itl {


  template <class Vector, class NonlinearFunction>
    class matrix_free_operator {
    
    public:
    
    typedef typename itl_traits<Vector>::value_type value_type;
    typedef typename itl_traits<Vector>::size_type  size_type;
      
    inline matrix_free_operator(NonlinearFunction f,
				const Vector& x,
				const Vector& z) : 
      f0(f), x0(size(x)), z0(size(z)), tmp0(size(x)) {
      itl::copy(x, x0);
      itl::copy(z, z0);
      
      sigma = 1.e-6 * itl::two_norm(x) + 1.e-8;
    }
    
    template<class VectorX, class VectorY>
      inline void
      apply(const VectorX& x, VectorY& y) const {
      
      // tmp0 <- x0 + sigma*x
      itl::add(x0, itl::scaled(x, sigma), tmp0);
      
      // y <- f(tmp0)
      f0(tmp0, y);
      
      // y <- y - f(x0)
      itl::add(itl::scaled(z0, -1.0), y);
      
      itl::scale(y, 1.0/sigma);
    }
    
    template<class VectorX, class VectorY, class VectorZ>
      inline void
      apply(const VectorX& x, const VectorY& y, VectorZ& z) const {
      
      
      // tmp0 <- x0 + sigma*x
      itl::add(x0, itl::scaled(x, sigma), tmp0);
      
      // z <- f(tmp0)
      f0(tmp0, z);
      
      itl::add(itl::scaled(z0, -1.0), z);
      
      itl::scale(z, 1.0/sigma);

      itl::add(y, z);
    }
    
    private:
    NonlinearFunction f0;
    Vector x0;
    Vector z0;
    mutable Vector tmp0;
    double sigma;
  };
  

  // Overload mult for the matrix_free_operator
  template <class Vector, class NonlinearFunction, class VectorX, 
            class VectorY>
  inline void
  mult(const matrix_free_operator<Vector, NonlinearFunction>& A,
       const VectorX& x, VectorY& y)
  {
    A.apply(x, y);
  }
  template <class Vector, class NonlinearFunction,
            class VectorW, class VectorX, class VectorY>
  inline void
  mult(const matrix_free_operator<Vector, NonlinearFunction>& A,
       const VectorW& w, const VectorX& x, VectorY& y)
  {
    A.apply(w, x, y);
  }

} //namespace itl

#endif //MATRIX_FREE_OPERATOR_H
