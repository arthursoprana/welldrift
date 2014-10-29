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

#ifndef ITL_TFQMR_H
#define ITL_TFQMR_H


#include "itl/itl.h"

namespace itl {

  //: Transpose Free Quasi-Minimal Residual
  //
  //  Transpose free QMR. First solve Q_1 A Q_2 x = Q_1 b. Then,
  //  return z which is Q_2 x. Here Q1 and Q2 are precondtioners.
  //  Suppose M is about equal to A and M = M_1 * M_2, then
  //  Q_1 = M_1^{-1} and Q_2 = M_2^{-1}
  //  <p>
  //  The residual holds |b - A * x_m| < sqrt{m+1} * tau_m. 
  //  The algorithm check the latter to see if convergence arrives instead of 
  //  checking real residual.
  //  <p>
  //<table align=center border=1>
  // <tr><td> return value </td>   <td>   meaning </td> </tr>
  // <tr><td>      0   </td><td>   convergence within maximum iterations </td> </tr>
  // <tr><td>      1   </td><td>     no convergence after maximum iterations</td> </tr>
  //  <tr><td>     2  </td><td>      breakdown in       tau </td> </tr>
  //  <tr><td>     3  </td><td>      breakdown in       alpha </td> </tr>
  //  <tr><td>     4  </td><td>      breakdown in       gamma</td> </tr>
  //  <tr><td>     5  </td><td>      breakdown in       rho</td> </tr>
  // </table>
  //
  //  <p>
  //  See: R. W. Freund, A Transpose-Free Quasi-Minimal Residual algorithm for 
  //  non-Hermitian linear system. SIAM J. on Sci. Comp. 14(1993), pp. 470-482
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: tfqmr.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: Vector - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold, SSOR or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot,add,scaled,two_norm */
template <class Matrix, class Vector, class VectorB, class Precond1, class Precond2, class Iteration>
int tfqmr(const Matrix& A, Vector& x, const VectorB& b,
	  const Precond1& M1, const Precond2& M2, Iteration& iter)
{
  typedef typename itl_traits<Vector>::value_type Real;
  typedef Vector TmpVec;
  TmpVec tmp(size(x)), r0(size(x)), v(size(x));

  Real sigma, alpha, c, kappa, beta;

  TmpVec  h(size(x));

  //x is initial value

  // 1. r0 = Q1 (b - A Q2 x)
  itl::solve(M2, x, r0);
  itl::mult(A, r0, tmp);
  itl::add(b, itl::scaled(tmp, -1.), tmp);
  itl::solve(M1, tmp, r0);

  // 2. w=y=r
  TmpVec w(size(x));
  itl::copy(r0, w);
  TmpVec y1(size(x));
  itl::copy(r0, y1);

  // 3. g=v=Q1AQ2y
  itl::solve(M2, y1, v);
  itl::mult(A, v, tmp);
  itl::solve(M1, tmp, v);

  TmpVec g(size(x));
  itl::copy(v, g);

  // 4. d=0
  TmpVec d(size(x));

  // 5. tau=||r||2
  Real tau = itl::two_norm(r0);

  // 6. theta=eta=0
  Real theta = 0.0;
  Real eta = 0.0;

  // 7. rtilde=r
  TmpVec rtilde(size(x));
  itl::copy(r0, rtilde);

  // 8. rho=dot(rtilde,r)
  Real rho = itl::dot(rtilde, r0);
  Real rho0 = rho;
  TmpVec y0(size(x));
  for (;;) {
    // 9. 10. 11.
    // sigma=dot(rtilde,v)
    // alpha=rho/sigma 
    // y2k=y(2k-1)-alpha*v
    sigma = itl::dot(rtilde, v);

    if (sigma==0.) { iter.fail(5, "tfqmr breakdown: sigma=0"); break; }
    alpha = rho / sigma;
    
    //y0 = y1 - alpha * v;
    itl::add(y1, itl::scaled(v, -alpha), y0);

    // 12. h=Q1*A*Q2*y
    itl::solve(M2, y0, h);
    itl::mult(A, h, tmp);
    itl::solve(M1, tmp, h);

    //split the loop of "for m = 2k-1, 2k" 

    //The first one
    // 13. w=w-alpha*Q1AQ2y0
    //w = w - alpha * g;
    itl::add(w, itl::scaled(g, -alpha), w);
    // 18. d=y0+((theta0^2)*eta0/alpha)*d         //need check breakdown
    if (alpha==0.) { iter.fail(3, "tfqmr breakdown: alpha=0"); break; }
    //d = y1 + ( theta * theta * eta / alpha ) * d;
    itl::add(y1, itl::scaled(d, theta * theta * eta / alpha), d);

    // 14. theta=||w||_2/tau0       //need check breakdown
    if (tau==0.) { iter.fail(2, "tfqmr breakdown: tau=0"); break; }
    theta  = itl::two_norm(w) / tau;
    
    // 15. c=1/sqrt(1+theta^2)
    c = 1. / sqrt(1. + theta * theta);

    // 16. tau=tau0*theta*c
    tau = tau * c * theta;

    // 17.  eta=(c^2)*alpha
    eta = c * c * alpha;

    // 19. x=x+eta*d
    //x += eta * d;
    itl::add(x, itl::scaled(d, eta), x);
    // 20. kappa=tau*sqrt(m+1)
    kappa = tau * sqrt( 2.* (iter.iterations()+1) );

    // 21. check stopping criterion
    if ( iter.finished(kappa) ) {
      //before return, transform x to the solution of Ax = b
      solve(M2, x, tmp);
      itl::copy(tmp, x);
      break;
    }
    //g = h;
    itl::copy(h, g);
    //The second one

    // 13. w=w-alpha*Q1AQ2y0
    //w = w - alpha * g;
    itl::add(w, itl::scaled(g, -alpha), w);
    // 18. d=y0+((theta0^2)*eta0/alpha)*d
    if (alpha==0.) { iter.fail(3,"tfqmr breakdown: alpha=0"); break; }
    //d = y0 + ( theta * theta * eta / alpha ) * d;
    itl::add(y0, itl::scaled(d,  theta * theta * eta / alpha), d);
    // 14. theta=||w||_2/tau0
    if (tau==0.) { iter.fail(2, "tfqmr breakdown: tau=0"); break; }
    theta = itl::two_norm(w) / tau;
    
    // 15. c=1/sqrt(1+theta^2)
    c = 1. / sqrt(1. + theta * theta);

    // 16. tau=tau0*theta*c
    tau = tau * c * theta;

    // 17.  eta=(c^2)*alpha
    eta = c * c * alpha;

    // 19. x=x+eta*d
    //x += eta * d;
    itl::add(x, itl::scaled(d, eta), x);

    // 20. kappa=tau*sqrt(m+1)
    kappa = tau * sqrt(2.* (iter.iterations()+1)  + 1.);

    // 21. check stopping criterion
    if ( iter.finished(kappa) ) {
      itl::solve(M2, x, tmp);
      itl::copy(tmp, x);
      break;
    }    

    // 22. rho=dot(rtilde,w)
    // 23. beta=rho/rho0                     //need check breakdown
    
    rho0 = rho;
    rho = itl::dot(rtilde, w);
    if (rho0==0.) { iter.fail(4, "tfqmr breakdown: beta=0"); break; }
    beta=rho/rho0;

    // 24. y=w+beta*y0
    //y1 = w + beta * y0;
    itl::add(w, itl::scaled(y0, beta), y1);

    // 25. g=Q1AQ2y
    //g = Q1 * ( A * ( Q2 * y1) );
    itl::solve(M2, y1, g);
    itl::mult(A, g, tmp);
    itl::solve(M1, tmp, g);

    // 26. v=Q1AQ2y+beta*(Q1AQ2y0+beta*v)

    //v = g + beta * ( h + beta * v );
    itl::add(g, itl::scaled(h, beta), itl::scaled(v, beta*beta), v);

    ++iter;
  }

  return iter.error_code();
}

}

#endif
