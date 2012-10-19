#ifndef H_WellSimulator_WELLSOLVER
#define H_WellSimulator_WELLSOLVER

// Includes MTL -----------------------------------------------------------------------------
#include <mtl/matrix.h>
#include <mtl/compressed2D.h>
#include <mtl/dense1D.h>
#include <mtl/mtl.h>
#include <mtl/utils.h>
#include <mtl/lu.h>
#include <mtl/envelope2D.h>

// Includes ITL -----------------------------------------------------------------------------
#include <itl/interface/mtl.h>
#include <itl/preconditioner/ssor.h>
#include <itl/itl.h>
#include <itl/preconditioner/ilu.h>
#include <itl/krylov/gmres.h>
#include <itl/krylov/qmr.h>

#include <memory>
#include <exception>
// Includes Boost
#include <BoostWrapper/SmartPointer.h>


namespace WellSimulator{
    //using namespace std;
    //using namespace std::tr1;
	typedef mtl::dense1D< double >  svector_type;
	typedef mtl::matrix< double, mtl::rectangle<>, mtl::dense<>, mtl::row_major>::type 					   dmatrix_type; // MATRIZ DENSA
	typedef mtl::matrix< double, mtl::rectangle<>, mtl::array< mtl::compressed<> >, mtl::row_major >::type smatrix_type; // MATRIZ ESPARSA
    typedef SharedPointer<svector_type> vector_ptr;
    typedef SharedPointer<smatrix_type> matrix_ptr;
    
    
    

	// WellSolver Class:
	// WellSolver has all methods to solve the linear system A.x = b
	class WellSolver{

	protected:
		vector_ptr m_x;		// Variables
		vector_ptr m_b;		// Source
		matrix_ptr m_A;		// Sparse Matrix

	public:
		WellSolver(){}
		WellSolver(matrix_ptr p_A, vector_ptr p_x, vector_ptr p_b){
			m_A = p_A;
			m_x = p_x;
			m_b = p_b;
		}

		~WellSolver(){}
		//void gmres_solve(smatrix_type &A, svector_type &x, svector_type &b);
		//void direct_solve	(dmatrix_type &A, svector_type &x, svector_type &b);

	};


	//void WellSolver::gmres_solve( smatrix_type &A, svector_type &x, svector_type &b ){

	//	int CONVERGENCE;
	//	try{			

	//		itl::ILU<smatrix_type> precond(A);
	//		// SSOR preconditioner			
	//		//itl::SSOR<smatrix_type> precond(A);		
	//		svector_type b2( A.ncols() );			
	//		itl::solve(precond(), b, b2); //gmres needs the preconditioned b to pass into iter object.
	//		//iteration
	//		int max_iter = 1000000;	//ex: 1000			
	//		itl::noisy_iteration<double> iter(b2, max_iter, 0.0, 1E-6);
	//		int restart = 10; //restart constant: 10
	//		// modified_gram_schmidt				
	//		itl::modified_gram_schmidt<svector_type> orth( restart, x.size() );			
	//		//gmres algorithm				
	//		CONVERGENCE = itl::gmres(A, x, b, precond(), restart, iter, orth); 

	//		if( CONVERGENCE == 1 )
	//			throw std::bad_exception( "Failure to converge..." );

	//	} catch( std::exception &e )
	//	{
	//		cerr << e.what() << "\n";
	//	}	
	//}

	//void WellSolver::direct_solve( dmatrix_type &A, svector_type &x, svector_type &b ){
	//		// Direct sOlver
	//		dmatrix_type LU(A.nrows(), A.ncols());
	//		svector_type pvector(A.nrows());
	//		copy(A, LU);
	//		lu_factor(LU, pvector);
	//		// call lu_solve with as many times for the same A as you want
	//		lu_solve(LU, pvector, b, x);

	//}


} // namespace WellSimulator

#endif // H_WellSimulator_WELLSOLVER