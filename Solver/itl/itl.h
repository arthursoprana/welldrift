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
#ifndef ITL_ITL_H
#define ITL_ITL_H


/**@name Iterative Template Library

The following is the requirements for the parameter types 
used in the ITL template functions.

The Matrix object must either follow the MTL requirements for a
matrix, be derived from the multiplier class (for matrix-free
multiplication), or have a specialized matvec::mult function defined
for it.

\begin{verbatim}
class Vector {
forward_iterator begin();
forward_iterator end();
};
\end{verbatim}

forward iterator refers to the iterator requirement
defined in the STL.

\begin{verbatim}
class Preconditioner {
void solve(const VectorX& x, VectorZ& z);
void trans_solve(const VectorX& x, VectorZ& z);
};
\end{verbatim}

The Preconditioner object performs a preconditioning operation based
on vector x and stores the result in vector z.  The trans solve()
method only need be defined when the preconditioner is used with an
iterative solver that requires it.

\begin{verbatim}
class Iteration {
bool first();
bool finished(const VectorX& r);
bool finished(const Real& r);
bool converged(const VectorX& r);
bool converged(const Real& r);
void operator++();
void failed();
int error_code();
int iterations();
Real resid();
};
\end{verbatim}

The Iteration object calculates whether the solution has reached the
desired accuracy, or whether the maximum number of iterations has
been reached. The method finished() checks both convergence and
number of iterations. The method converged() only checks
convergence. The error code() method is used to determine the return
value for the this iterative solver function. The first() method is
used to determine the first iteration of the loop.


For all algorithms, if the error\_code() is 0, it suggests the algorithm 
converges. Otherwise, if the error\_code() returns 1, it means the maximum 
number of iteration has been reached but the desired accuacy is not reached.
For other return codes, see the respective document.
*/

/**

@name  Utilities
@memo  Utility classes and routines for ITL

*/

#include <itl/itl_config.h>
#include <iostream>
#include <complex>
#include <string>

namespace itl {


	template <class Real>
	class basic_iteration {
	public:


		typedef Real real;


		template <class Vector>
			basic_iteration(const Vector& b, int max_iter_, Real t, Real a = Real(0))
			: error(0), i(0), normb_(std::fabs(itl::two_norm(b))), 
			max_iter(max_iter_), rtol_(t), atol_(a) { }

			basic_iteration(Real nb, int max_iter_, Real t, Real a = Real(0))
				: error(0), i(0), normb_(nb), max_iter(max_iter_), rtol_(t), atol_(a) {}


				template <class Vector>
					bool finished(const Vector& r) {
						Real normr_ = std::fabs(itl::two_norm(r)); 
						if (converged(normr_))
							return true;
						else if (i < max_iter)
							return false;
						else {
							error = 1;
							return true;
						}
					}


					bool finished(const Real& r) {
						if (converged(r))
							return true;
						else if (i < max_iter)
							return false;
						else {
							error = 1;
							return true;
						}
					}

					template <typename T>
						bool finished(const std::complex<T>& r) { 
							if (converged(std::fabs(r)))
								return true;
							else if (i < max_iter)
								return false;
							else {
								error = 1;
								return true;
							}
						}

						inline bool converged(const Real& r) {
							if (normb_ == 0)
								return r < atol_;  // ignore relative tolerance if |b| is zero
							resid_ = r / normb_;
							return (resid_ <= rtol_ || r < atol_); // relative or absolute tolerance.
						}

						inline void operator++() { ++i; }

						inline bool first() { return i == 0; }

						inline int error_code() { return error; }

						inline int iterations() { return i; }

						inline Real resid() { return resid_ * normb_; }

						inline Real normb() const { return normb_; }

						inline Real tol() { return rtol_; }
						inline Real atol() { return atol_; } 

						inline void fail(int err_code) { error = err_code; }

						inline void fail(int err_code, const std::string& msg)
						{ error = err_code; err_msg = msg; }

						inline void set(Real v) { normb_ = v; }

	protected:
		int error;
		int i;
		const Real normb_;
		int max_iter;
		Real rtol_;
		Real atol_;
		Real resid_;
		std::string err_msg;
	};

	template <class Real>
	class noisy_iteration : public basic_iteration<Real> {
		typedef basic_iteration<Real> super;
	public:

		template <class Vector>
			noisy_iteration(const Vector& b, int max_iter_, 
			Real tol_, Real atol_ = Real(0))
			: super(b, max_iter_, tol_, atol_) { }

			template <class Vector>
				bool finished(const Vector& r) {
					using std::cout;
					using std::endl;

					Real normr_ = std::fabs(itl::two_norm(r)); // CHANGED
					bool ret;
					if (converged(normr_))
						ret = true;
					else if (this->i < this->max_iter)
						ret = false;
					else {
						this->error = 1;
						ret = true;
					}
					//cout << "iteration " << this->i << ": resid " 
					//<< this->resid()
					//<< endl;
					
					return ret;
				}


				bool finished(const Real& r) {
					using std::cout;
					using std::endl;

					bool ret;
					if (converged(r))
						ret = true;
					else if (this->i < this->max_iter)
						ret = false;
					else {
						this->error = 1;
						ret = true;
					}
					//cout << "iteration " << this->i << ": resid " 
					//<< this->resid()
					//<< endl;
					
					return ret;
				}

				template <typename T>
					bool finished(const std::complex<T>& r) { //for the case of complex
						using std::cout;
						using std::endl;

						bool ret;
						if (converged(std::fabs(r)))
							ret = true;
						else if (this->ii < this->imax_iter)
							ret = false;
						else {
							this->error = 1;
							ret = true;
						}
						//cout << "iteration " << this->i << ": resid " 
						//<< this->resid() << endl;
						
						return ret;
					}

					int error_code() {
						using std::cout;
						using std::endl;

#ifdef DEBUG	
						cout << "finished! error code = " << this->error << endl;
						cout << this->iterations() << " iterations" << endl;
						cout << this->resid() << " is actual final residual. " << endl
							<< this->resid()/this->normb() << " is actual relative tolerance achieved. "
							<< endl;
						cout << "Relative tol: " << this->rtol_ << "  Absolute tol: " << this->atol_ << endl;
#endif
						return this->error;
					}

	};


	struct identity_preconditioner {
		identity_preconditioner operator()() const {
			identity_preconditioner p;
			return p;
		}

		identity_preconditioner left() const {
			identity_preconditioner p;
			return p;
		}

		identity_preconditioner right() const {
			identity_preconditioner p;
			return p;
		}

	};


	template <class VecX, class VecZ>
		inline void solve(const identity_preconditioner& M, const VecX& x, 
		const VecZ& z) {
			itl::copy(x, const_cast<VecZ&>(z));
		}

		template <class VecX, class VecZ>
			inline void trans_solve(const identity_preconditioner& M, 
			const VecX& x, const VecZ& z) {
				itl::copy(x, const_cast<VecZ&>(z));
			}

			template <class Preconditioner, class VecX, class VecZ>
				inline void 
				solve(const Preconditioner& M, const VecX& x, const VecZ& z) {
					M.solve(x, const_cast<VecZ&>(z));
				}

				template <class Preconditioner, class VecX, class VecZ>
					inline void 
					trans_solve(const Preconditioner& M, const VecX& x, const VecZ& z) {
						M.trans_solve(x, const_cast<VecZ&>(z));
					}

}

#endif
