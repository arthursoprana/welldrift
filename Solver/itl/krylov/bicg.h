// -*- c++ -*-
//
//     $Copyright$
//

#ifndef ITL_BICG_H
#define ITL_BICG_H

#include "itl/itl.h"

namespace itl {

  //: BiConjugate Gradient
  //
  //This solves the unsymmetric linear system Ax = b 
  //using the Preconditioned BiConjugate Gradient method.
  //<p>
  //A return value of 0 indicates convergence within the
  //maximum number of iterations (determined by the iter object).
  //A return value of 1 indicates a failure to converge.
  //<p>
  //<p>
  //See: R. Fletcher, Conjugate gradient methods for indefinite systems, 
  //In Numerical Analysis Dundee 1975, G. Watson, ed., Springer Verlag, 
  //Berlin, New York, 1976 pp. 73-89
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: bicg.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot,add,scaled,trans_mult */
template <class Matrix, class Vector, class VectorB, 
          class Preconditioner, class Iteration>
int
bicg(const Matrix& A, Vector& x, const VectorB& b,
     const Preconditioner& M, Iteration& iter)
{

  typename itl_traits<Vector>::value_type rho_1, rho_2, alpha, beta;
  typedef Vector TmpVec;
  TmpVec r(size(x)), z(size(x)), p(size(x)), q(size(x));
  TmpVec r_tilde(size(x)), z_tilde(size(x));
  TmpVec p_tilde(size(x)), q_tilde(size(x));
  
  itl::mult(A, itl::scaled(x, -1), b, r);	  
  itl::copy(r, r_tilde);	          

  while (! iter.finished(r)) {
    itl::solve(M, r, z);		
    itl::trans_solve(M, r_tilde, z_tilde);

    rho_1 = itl::dot(z, r_tilde);

    if (rho_1 == 0.) {
      iter.fail(2, "bicg breakdown");
      break;
    }

    if (iter.first()) {
      itl::copy(z, p);
      itl::copy(z_tilde, p_tilde);	  
    } else {
      beta = rho_1 / rho_2;
      itl::add(z, itl::scaled(p, beta), p); 
      itl::add(z_tilde, itl::scaled(p_tilde, beta), p_tilde);
                                  
    }
    
    itl::mult(A, p, q);		  
    itl::trans_mult(A, p_tilde, q_tilde);
    
    alpha = rho_1 / itl::dot(p_tilde, q);
    
    itl::add(x, itl::scaled(p, alpha), x);   
    itl::add(r, itl::scaled(q, -alpha), r);  
    itl::add(r_tilde, itl::scaled(q_tilde, -alpha), r_tilde);

    rho_2 = rho_1;

    ++iter;
  }

  return iter.error_code();
}

}

#endif
