#include <DriftFluxWell.h>



#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <ctime>
#include <BoostWrapper/Math.h>


//using namespace std;


// Namespace =======================================================================================
namespace WellSimulator {

	
	
	typedef char						string_type;
	typedef double						real_type;
	typedef unsigned int				uint_type;
	//typedef WellSimulator::WellVector   vector_type;
    typedef std::vector<double>				vector_type;
	typedef NodeCoordinates				coord_type;
//#define PI = 3.1415926535897932384626433832795;
    enum{I_DIRECTION, J_DIRECTION, K_DIRECTION};
    enum{WaterPhase, OilPhase, GasPhase, NumberOfPhases = 3};

	class Timer{
	public:
        Timer() : m_print_time(false) {}

		inline void start(){
			this->old = (double)clock()/CLOCKS_PER_SEC;
		}

		inline void stop(){
			this->now = (double)clock()/CLOCKS_PER_SEC;
		}

		inline double elapsed(){
			return this->now - this->old;
		}

        inline void print(std::string p_msg){
            if(m_print_time){
                std::cout << p_msg << this->elapsed();
            }            
        }

        inline void enable_print_time(){
            m_print_time = true;
        }

        inline void disable_print_time(){
            m_print_time = false;
        }

	protected:
        bool m_print_time;
		double now;
		double old;
	};

    real_type calculate_inclination_correction(real_type p_inclination){
        return - (  boost::math::cbrt(cos(p_inclination)) 
                  * pow( abs(cos(p_inclination)), 1.0/6.0 ) 
                  * pow(1.0 + sin(p_inclination),2.0) 
                 );
    }
	
	DriftFluxWell::DriftFluxWell()
	{
	}
	DriftFluxWell::DriftFluxWell(
								 const uint_type& p_nnodes,
								 const real_type& p_radius
								)
								: m_oil_velocity	( p_nnodes, 0 ),
								  m_water_velocity	( p_nnodes, 0 ),
								  m_gas_velocity	( p_nnodes, 0 ),
								  m_mean_velocity	( p_nnodes, 0 ),
								  m_oil_vol_frac	( p_nnodes, 0 ),
								  m_water_vol_frac	( p_nnodes, 0 ),
								  m_gas_vol_frac	( p_nnodes, 0 ),
								  m_gas_vol_frac_old	( p_nnodes, 0 ),
								  m_oil_vol_frac_old	( p_nnodes, 0 ),
								  m_water_vol_frac_old	( p_nnodes, 0 ),
								  m_pressure_old		( p_nnodes, 0 ),
								  m_mean_velocity_old	( p_nnodes, 0 ),								  
								  m_id				( p_nnodes ),
								  m_gravity			( 3, 0 ),
								  m_delta			( total_var, 0 ),
                                  m_matrix   (new smatrix_type(total_var*p_nnodes,total_var*p_nnodes)),
                                  m_variables(new svector_type(total_var*p_nnodes)),
                                  m_source   (new svector_type(total_var*p_nnodes)),
                                  m_has_inclination_correction(true)
	{					 
	    

		for( uint_type i = 0; i < m_id.size(); ++i )
		{
			this->m_id[ i ].resize( total_var );
			m_id[ i ][ P ]	     = total_var*i;
			m_id[ i ][ alpha_g ] = total_var*i + 1;
			m_id[ i ][ alpha_o ] = total_var*i + 2;
			m_id[ i ][ v ]       = total_var*i + 3;
		}
		this->m_coordinates.resize( p_nnodes );
		
		this->m_pressure.resize(p_nnodes, 100000.0);
		this->m_nnodes = p_nnodes;
		this->m_radius = p_radius;		
	}


	
	DriftFluxWell::~DriftFluxWell()
	{
	}

    // Initialize for RESERVOIR SOLVER ---->> FDarcy
    void DriftFluxWell::Initialize(uint_type p_number_of_completions, uint_type p_direction, real_type p_well_radius, real_type p_BHPressure, vector_type p_well_length){


        // first volume has no completions
/*                            WELL
             ___________________________________________
            |   |     |     |     |     |     |     |   |
            x   |  x  |  x  |  x  |  x  |  x  |  x  |   x
            |___|_____|_____|_____|_____|_____|_____|___|
                   ^     ^     ^     ^     ^     ^      ^
                   |     |     |     |     |     |      |  Lateral Mass Inflow
*/
        uint_type well_size = p_number_of_completions+1;
        m_oil_velocity.resize       ( well_size, 0.0 );	
        m_water_velocity.resize	    ( well_size, 0.0 );
        m_gas_velocity.resize	    ( well_size, 0.0 );
        m_mean_velocity.resize	    ( well_size, 0.0 );
        m_oil_vol_frac.resize	    ( well_size, 0.0 );
        m_water_vol_frac.resize	    ( well_size, 0.0 );
        m_gas_vol_frac.resize	    ( well_size, 0.0 );
        m_gas_vol_frac_old.resize	( well_size, 0.0 );
        m_oil_vol_frac_old.resize	( well_size, 0.0 );
        m_water_vol_frac_old.resize	( well_size, 0.0 );
        m_pressure_old.resize		( well_size, 0.0 );
        m_mean_velocity_old.resize  ( well_size, 0.0 );								  
        m_id.resize				    ( well_size );
        m_gravity.resize			( 3, 0.0 );
        m_delta.resize			    ( total_var, 0 );
        m_total_production.resize   ( NumberOfPhases, 0);
       
        m_oil_flow.resize(well_size, MakeShared<ConstantInflow>(0.0));
        m_gas_flow.resize(well_size, MakeShared<ConstantInflow>(0.0));
        m_water_flow.resize(well_size, MakeShared<ConstantInflow>(0.0));

        m_matrix	= SharedPointer<smatrix_type>( new smatrix_type(total_var*well_size,total_var*well_size) );
        m_variables = SharedPointer<svector_type>( new svector_type(total_var*well_size) );
        m_source	= SharedPointer<svector_type>( new svector_type(total_var*well_size) );

        for( uint_type i = 0; i < m_id.size(); ++i )
        {
            this->m_id[ i ].resize( total_var );
            m_id[ i ][ P ]	     = total_var*i;
            m_id[ i ][ alpha_g ] = total_var*i + 1;
            m_id[ i ][ alpha_o ] = total_var*i + 2;
            m_id[ i ][ v ]       = total_var*i + 3;
        }
        this->m_coordinates.resize( well_size );

        this->m_pressure.resize(well_size , p_BHPressure);
        
        this->set_bottom_pressure(p_BHPressure);
        this->m_nnodes = well_size;
        this->m_radius = p_well_radius;

        this->set_with_gas (true);
        this->set_mass_flux(false);
        unsigned NODES = well_size;
       
        ////// CREATING WELL COORDINATE VECTOR //
        std::vector<coord_type> COORD_VECTOR(NODES);
        COORD_VECTOR[ 0 ][0] = 0;
        COORD_VECTOR[ 0 ][1] = 0;
        COORD_VECTOR[ 0 ][2] = 0;
        if(p_direction == I_DIRECTION){
            COORD_VECTOR[ 1 ][0] = COORD_VECTOR[ 0 ][0] + 0.1;
            COORD_VECTOR[ 1 ][1] = 0.0;
            COORD_VECTOR[ 1 ][2] = 0.0;
        }
        if(p_direction == J_DIRECTION){
            COORD_VECTOR[ 1 ][0] = 0.0;
            COORD_VECTOR[ 1 ][1] = COORD_VECTOR[ 0 ][1] + 0.1;
            COORD_VECTOR[ 1 ][2] = 0.0;
        }
        if(p_direction == K_DIRECTION){
            COORD_VECTOR[ 1 ][0] = 0.0;
            COORD_VECTOR[ 1 ][1] = 0.0;
            COORD_VECTOR[ 1 ][2] = COORD_VECTOR[ 0 ][2] + 0.1;
        }

        for( unsigned i = 2; i < NODES; ++i ){
            double  LENGTH = p_well_length[ i-1 ]; // [ m ]
            double  ds = LENGTH;
            if(p_direction == I_DIRECTION){
                COORD_VECTOR[ i ][0] = COORD_VECTOR[ i-1 ][0] + ds;
                COORD_VECTOR[ i ][1] = 0.0;
                COORD_VECTOR[ i ][2] = 0.0;
            }
            if(p_direction == J_DIRECTION){
                COORD_VECTOR[ i ][0] = 0.0;
                COORD_VECTOR[ i ][1] = COORD_VECTOR[ i-1 ][1] + ds;
                COORD_VECTOR[ i ][2] = 0.0;
            }
            if(p_direction == K_DIRECTION){
                COORD_VECTOR[ i ][0] = 0.0;
                COORD_VECTOR[ i ][1] = 0.0;
                COORD_VECTOR[ i ][2] = COORD_VECTOR[ i-1 ][2] + ds;
            }
           		
        }
        ////// CREATED....

        double PROFILE_PARAM_C_0 = 1.2;
        double delta_t = 1.0;
        double TOLERANCE = 1.0e-4;

        double alpha_G = 0.000001;
        double alpha_O = 0.999998;
        double alpha_W = 1.0 - (alpha_G+alpha_O);
        this->set_constant_vol_frac( alpha_O, alpha_G, alpha_W );	
        this->set_dt( delta_t );
        this->set_final_timestep( 1000 );
        this->set_C_0( PROFILE_PARAM_C_0 );

        this->set_delta( 0.0001, 0.0001, 0.0001, 0.0001 );
        this->set_constant_pressure( p_BHPressure );
        this->set_heel_pressure( p_BHPressure );
        this->set_constant_velocity( 0.0 );


        this->set_boundary_velocity( 0.0 );        
        this->set_newton_criteria( TOLERANCE );

        inflow_vector_type inflow_gas(NODES, MakeShared<ConstantInflow>(0.0));
        inflow_vector_type inflow_oil(NODES, MakeShared<ConstantInflow>(0.0));
        inflow_vector_type inflow_water(NODES, MakeShared<ConstantInflow>(0.0));

        this->initialize_flow(inflow_oil,inflow_water,inflow_gas);

        this->set_coordinates(COORD_VECTOR);
        this->set_gravity( 0., 0., 9.8 ); 
    }
	void DriftFluxWell::set_delta( real_type p_delta_P, real_type p_delta_alpha_g, real_type p_delta_alpha_o, real_type p_delta_v )
	{
		this->m_delta[ P ]			= p_delta_P;
		this->m_delta[ alpha_g ]	= p_delta_alpha_g;
		this->m_delta[ alpha_o ]	= p_delta_alpha_o;
		this->m_delta[ v ]			= p_delta_v;
	}
	void DriftFluxWell::set_dt( real_type p_dt )
	{
		this->m_dt = p_dt;
	}

    real_type norm(vector_type p_vector){
        real_type Norm = 0;
        for( uint_type i = 0; i < p_vector.size(); ++i )
            Norm += p_vector[ i ]*p_vector[ i ];
        return sqrt(Norm);
    }
	void DriftFluxWell::set_bottom_pressure( real_type p_pressure ){
		this->m_pressure[ 0 ] = p_pressure;
	}

	
	void DriftFluxWell::set_with_gas(bool p_choice){
		this->m_with_gas = p_choice;
	}


	real_type DriftFluxWell::dt(){
		return this->m_dt;
	}

	real_type DriftFluxWell::ksi( real_type p_velocity ){
		return p_velocity >= 0 ? 0.5 : -0.5;
        //return 0.0;
	}
	
	real_type DriftFluxWell::Volume( real_type dS ){
		return this->area()*dS;
	}

	real_type DriftFluxWell::area(){
		return this->m_radius*this->m_radius*PI;
	}

	real_type DriftFluxWell::gas_density(real_type p_pressure){
		//return 5.9733267069e-6*p_pressure; // Let's assume that for the moment... rho = P/(R.T), 
       // return p_pressure > 0.0 ? p_pressure/(316.0*316.0) : this->m_HEEL_PRESSURE/(316.0*316.0);
        //return p_pressure/(316.0*316.0);
        return m_gas_density_model->compute_density(p_pressure);
        //return 10.0;
		//return 1.1245;
	}									   // where R = 518.3 J/(Kg.K)and T = 323 K	

	real_type DriftFluxWell::liquid_density(
											real_type p_oil_vol_frac,
											real_type p_water_vol_frac,
											real_type p_pressure
											)
	{    
        if(p_water_vol_frac < 0.0) p_water_vol_frac = 0.0;
        if(p_oil_vol_frac   < 0.0) p_oil_vol_frac   = 0.0;
        real_type den = p_oil_vol_frac + p_water_vol_frac; // denominator
        if( abs(den) < 1.0e-12){
            return 0.5*this->oil_density  ( p_pressure ) + 0.5*this->water_density( p_pressure );
        }
        else{
            
            return (
                    ( p_oil_vol_frac  *this->oil_density  ( p_pressure ) 
                    + p_water_vol_frac*this->water_density( p_pressure )) / den				
                   );
        }
		
		
	}

	real_type DriftFluxWell::oil_density( real_type p_pressure ){
		//return 1000. + (p_pressure - 100000.0)/1.0e6; // that as well...
        //return 800.0;
        return m_oil_density_model->compute_density(p_pressure);
	}

	real_type DriftFluxWell::water_density( real_type p_pressure ){
		//return 1000. + (p_pressure - 100000.0)/1.0e6; // that as well...
        //return 1000.0;
        return m_water_density_model->compute_density(p_pressure);
	}


	real_type DriftFluxWell::mean_density( 
										  const real_type& p_oil_vol_frac,
										  const real_type& p_water_vol_frac,
										  const real_type& p_gas_vol_frac,
										  const real_type& p_pressure
										  )
    {         
       
		return (
			    p_oil_vol_frac   * this->oil_density  (p_pressure) +  
				p_water_vol_frac * this->water_density(p_pressure) +  
				p_gas_vol_frac   * this->gas_density  (p_pressure) 
			   );           
		
	}

	real_type DriftFluxWell::friction_factor( real_type p_reynolds ){
		//if( p_reynolds == 0.0 )
		//	return 0.0;
  //      else{
  //          double e = 1.0e-4; // FoFo
  //          double D = 2.0*m_radius;                
  //          double l = log(pow(7.0/p_reynolds,0.9) + 0.27*e/D);             
  //          double A = pow(-2.547*l,16.0);              
  //          double B = pow(37530.0/p_reynolds,16.0);                 
  //          double f0 = 8.0*pow(pow(8.0/p_reynolds,12.0)+1.0/pow((A+B),1.5),1.0/12.0);             
  //          return f0;
  //      }
        return p_reynolds == 0.0 ? 0.0 : abs(64/p_reynolds); // Laminar AtTheMoment...			
	}

	real_type DriftFluxWell::mean_velocity( uint_type p_index ){
		return this->m_mean_velocity[ p_index ];		
	}
	
	real_type DriftFluxWell::C_0(){
		return this->m_profile_parameter_C_0;
	}

	void DriftFluxWell::set_C_0(real_type p_C_0){
		this->m_profile_parameter_C_0 = p_C_0;
	}


	

	real_type DriftFluxWell::v_drift_flux( real_type p_gas_vol_frac, real_type p_pressure, real_type p_liquid_density )
	{
		real_type C  = this->m_profile_parameter_C_0;
		real_type rho_g = this->gas_density( p_pressure );
		real_type rho_l = p_liquid_density;
		real_type alphaC_0 = p_gas_vol_frac*C;
        
        

		return 0.05; // SImple scheme
		//return 0.;  // homogeneous
	}

	real_type DriftFluxWell::mod_v_drift_flux(											
											  real_type p_mean_velocity, 
											  real_type p_gas_vol_frac, 
											  real_type p_oil_vol_frac,
											  real_type p_water_vol_frac,
											  real_type p_pressure  
											  )
	{	
        if(p_pressure < 0.0 || p_gas_vol_frac < 0.0 || p_oil_vol_frac < 0.0 || p_gas_vol_frac > 1.0 || p_oil_vol_frac > 1.0){
            m_convergence_status = true;
            return 0.0;
        }
		real_type rho_m = this->mean_density( p_oil_vol_frac, p_water_vol_frac, p_gas_vol_frac, p_pressure );
		real_type rho_l = this->liquid_density( p_oil_vol_frac, p_water_vol_frac, p_pressure);
		real_type rho_g = this->gas_density( p_pressure );

        float64 sigma_go = 0.0;
        float64 sigma_gw = 0.0;
        float64 interfacial_tension = 0.0;         

        sigma_go = m_gas_oil_interfacial_tension_model->compute_interfacial_tension  ( ( p_pressure < 0.0) ? 0.0 : p_pressure );
        sigma_gw = m_gas_water_interfacial_tension_model->compute_interfacial_tension( ( p_pressure < 0.0) ? 0.0 : p_pressure );
         
        real_type den = p_oil_vol_frac + p_water_vol_frac;        
        if( abs(den) < 1.0e-12 ){
            interfacial_tension   = 0.5*(sigma_go + sigma_gw);
        }
        else{
            interfacial_tension = (p_oil_vol_frac*sigma_go + p_water_vol_frac*sigma_gw )/den;
        }

        float64 min_interfacial_tension = WellConstants::convert_Dynes_per_cm_to_Pa_m();
        if( interfacial_tension < min_interfacial_tension ){
            interfacial_tension = min_interfacial_tension;
        }

        
        real_type D_hat = sqrt( gravity()*(rho_l - rho_g)/interfacial_tension )*2.0*m_radius;
        real_type Ku;
        if(D_hat <= 2.0){
            Ku = 0.0;
        }else if(D_hat >= 50){
            Ku = 3.2;
        }else{
            Ku = 2.684*exp(0.003669*D_hat) - 3.847*exp(-0.1853*D_hat);
        }
        real_type Vc = pow( interfacial_tension*gravity()*(rho_l - rho_g)/(rho_l*rho_l) , 0.25 );
        m_gas_liquid_profile_parameter_model->set_flooding_velocity( Ku*sqrt(rho_l/rho_g)*Vc );
        m_gas_liquid_profile_parameter_model->set_mixture_velocity(p_mean_velocity);
        m_gas_liquid_profile_parameter_model->set_volume_fraction(p_gas_vol_frac);
        real_type C_0_gl = m_gas_liquid_profile_parameter_model->compute_profile_parameter();
        m_gas_liquid_drift_velocity_model->set_volume_fraction(p_gas_vol_frac);
        m_gas_liquid_drift_velocity_model->set_dispersed_density(rho_g);
        m_gas_liquid_drift_velocity_model->set_not_dispersed_density(rho_l);
        m_gas_liquid_drift_velocity_model->set_Ku_critical(Ku);
        m_gas_liquid_drift_velocity_model->set_characteristic_velocity(Vc);
        m_gas_liquid_drift_velocity_model->set_profile_parameter(C_0_gl);
		real_type v_d   = m_gas_liquid_drift_velocity_model->compute_drift_velocity();

        if(m_has_inclination_correction){
            v_d *= calculate_inclination_correction(m_well_inclination);
        }        

		return (v_d+(C_0_gl-1)*p_mean_velocity)/(1-(C_0_gl-1)*p_gas_vol_frac*(rho_l-rho_g)/rho_m );		
	}

	real_type DriftFluxWell::mod_v_drift_flux_ow(											
			real_type p_mean_velocity, 
			real_type p_gas_vol_frac, 
			real_type p_oil_vol_frac,
			real_type p_water_vol_frac,
			real_type p_pressure  
			)
	{	
        if(p_pressure < 0.0 || p_gas_vol_frac < 0.0 || p_oil_vol_frac < 0.0 || p_gas_vol_frac > 1.0 || p_oil_vol_frac > 1.0){
            m_convergence_status = true;
            return 0.0;
        }
		real_type rho_o = this->oil_density( p_pressure );
		real_type rho_w = this->water_density( p_pressure );
        real_type rho_l = this->liquid_density( p_oil_vol_frac, p_water_vol_frac, p_pressure );
        real_type alpha_ol = p_oil_vol_frac/(p_oil_vol_frac + p_water_vol_frac + 1.0e-20);
        m_oil_water_profile_parameter_model->set_volume_fraction(alpha_ol);       
        real_type C_0_ow = m_oil_water_profile_parameter_model->compute_profile_parameter();
        m_oil_water_drift_velocity_model->set_volume_fraction(alpha_ol);

        float64 sigma_go = 0.0;
        float64 sigma_gw = 0.0;
        float64 interfacial_tension = 0.0;         

        sigma_go = m_gas_oil_interfacial_tension_model->compute_interfacial_tension  ( ( p_pressure < 0.0) ? 0.0 : p_pressure );
        sigma_gw = m_gas_water_interfacial_tension_model->compute_interfacial_tension( ( p_pressure < 0.0) ? 0.0 : p_pressure );
        
        float64 min_interfacial_tension = WellConstants::convert_Dynes_per_cm_to_Pa_m();
        interfacial_tension = sigma_gw - sigma_go;
        if( interfacial_tension < min_interfacial_tension ){
            interfacial_tension = min_interfacial_tension;
        }

        real_type Vc = pow( interfacial_tension*gravity()*(rho_w - rho_o)/(rho_w*rho_w) , 0.25 );
        m_oil_water_drift_velocity_model->set_characteristic_velocity(Vc);
		real_type v_d   = m_oil_water_drift_velocity_model->compute_drift_velocity();

        if(m_has_inclination_correction){
            v_d *= calculate_inclination_correction(m_well_inclination);
        }

        float64 a3 = 0.017*exp( pow(m_well_inclination,3.28) );
        if (p_gas_vol_frac < a3){
            v_d = (1.0 - p_gas_vol_frac/a3)*v_d;
        }
        else{
            v_d = 0.0;
        }

        return (v_d + (C_0_ow - 1.0)*p_mean_velocity)/( 1.0 - (C_0_ow - 1.0)*p_oil_vol_frac*(rho_w-rho_o)/rho_l );
	}



	void DriftFluxWell::set_gravity( real_type p_valueX = 0., real_type p_valueY = 0., real_type p_valueZ = 9.8 )
	{
		this->m_gravity[ 0 ] = p_valueX;
		this->m_gravity[ 1 ] = p_valueY;
		this->m_gravity[ 2 ] = p_valueZ;
	}
	real_type DriftFluxWell::gravity(){
		return sqrt( m_gravity[ 0 ]*m_gravity[ 0 ] + m_gravity[ 1 ]*m_gravity[ 1 ] + m_gravity[ 2 ]*m_gravity[ 2 ] );
	}

	void DriftFluxWell::set_mean_velocity()
	{		
		for( uint_type i = 0; i < this->m_nnodes; ++i )
		{			
			real_type   P = this->m_pressure[ i ];
			real_type a_g = this->m_gas_vol_frac[ i ];
			real_type a_o = this->m_oil_vol_frac[ i ];
			real_type a_w = this->m_water_vol_frac[ i ];
			real_type d_g = this->gas_density( P );
			real_type d_o = this->oil_density( P );
			real_type d_w = this->water_density( P );
			real_type v_g = this->m_gas_velocity[ i ];
			real_type v_o = this->m_oil_velocity[ i ];
			real_type v_w = this->m_water_velocity [ i ];
			real_type d_m = this->mean_density( a_o, a_w, a_g, P );		
			this->m_mean_velocity[ i ] = (a_g*d_g*v_g + a_o*d_o*v_o + a_w*d_w*v_w) / d_m; 										  
		}
	}

	void DriftFluxWell::set_constant_oil_vol_frac( real_type p_value )
	{
		for( uint_type i = 0; i < this->m_nnodes; ++i )
		{		
			this->m_oil_vol_frac[ i ] = p_value;								  
		}
	}

	void DriftFluxWell::set_constant_vol_frac( 
											  real_type p_oil_vol_frac, 
											  real_type p_gas_vol_frac, 
											  real_type p_water_vol_frac 
											  )
	{
		for( uint_type i = 0; i < this->m_nnodes; ++i )
		{		
			this->m_oil_vol_frac[ i ]   = p_oil_vol_frac;
			this->m_gas_vol_frac[ i ]   = p_gas_vol_frac;
			this->m_water_vol_frac[ i ] = p_water_vol_frac;
		}
	}

	void DriftFluxWell::set_constant_pressure( real_type p_pressure )
	{
		for( uint_type i = 0; i < this->m_nnodes; ++i )
		{		
			this->m_pressure[ i ] = p_pressure;
		}
	}

	void DriftFluxWell::set_constant_velocity( real_type p_velocity )
	{
		for( uint_type i = 0; i < this->m_nnodes; ++i )
		{		
			this->m_mean_velocity[ i ] = p_velocity;
		}
	}
	
	uint_type DriftFluxWell::id( uint_type p_node , uint_type p_variable ){
		return this->m_id[ p_node ][ p_variable ];
	}
	
	real_type DriftFluxWell::segment_length( coord_type p_coord_i, coord_type p_coord_j )
	{		
		return sqrt( 
					(p_coord_i.getX()-p_coord_j.getX())*(p_coord_i.getX()-p_coord_j.getX())
				   +(p_coord_i.getY()-p_coord_j.getY())*(p_coord_i.getY()-p_coord_j.getY()) 
				   +(p_coord_i.getZ()-p_coord_j.getZ())*(p_coord_i.getZ()-p_coord_j.getZ())
				   );
	}
	real_type DriftFluxWell::dot( vector_type& p_vec1, vector_type& p_vec2 )
	{
		real_type sum = 0;
		if (p_vec1.size() != p_vec2.size()){
			throw std::runtime_error("Vector's length do not match...");
		}
		for( uint_type i = 0; i < p_vec1.size(); ++i )
			sum += p_vec1[ i ]*p_vec2[ i ];
		return sum;           
	}

	void DriftFluxWell::set_newton_criteria(real_type p_tolerance){
		this->NEWTON_CRIT = p_tolerance;
	}

	void DriftFluxWell::set_final_timestep(uint_type p_final_timestep ){
		this->m_FINAL_TIMESTEP = p_final_timestep;
	}
	real_type DriftFluxWell::gas_viscosity( real_type p_pressure ){
		//return 5.0e-6;
        return m_gas_viscosity_model->compute_viscosity(p_pressure);
	}

	real_type DriftFluxWell::oil_viscosity( real_type p_pressure ){
		//return 0.05;
        return m_oil_viscosity_model->compute_viscosity(p_pressure);
	}

	real_type DriftFluxWell::water_viscosity( real_type p_pressure ){
		//return 0.05;
        return m_water_viscosity_model->compute_viscosity(p_pressure);
	}


	real_type DriftFluxWell::R_m(
								 real_type p_pressureW,
								 real_type p_pressureP,
								 real_type p_pressureE,
								 real_type p_gas_vol_fracW,
								 real_type p_gas_vol_fracP,
								 real_type p_gas_vol_fracE,
								 real_type p_oil_vol_fracW,
								 real_type p_oil_vol_fracP,
								 real_type p_oil_vol_fracE,
								 real_type p_velocityW,
								 real_type p_velocityP,
								 uint_type p_node,
								 string_type	   position = 'C'
								 )
	{
		
		switch ( position )
		{
		
		case 'L':
			{
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.;
				real_type dS  = dSw + dSe;
				real_type dV = this->Volume( dS );

				real_type water_vol_fracOld = 1.0 - (m_oil_vol_frac_old[ p_node ] + m_gas_vol_frac_old[ p_node ]);
				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);

				real_type rho_P_old = this->mean_density(m_oil_vol_frac_old[ p_node ], water_vol_fracOld, m_gas_vol_frac_old[ p_node ], m_pressure_old[ p_node ]);
                
                real_type rho_P = this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);			
				real_type rho_W = this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);

                real_type rhoG_P		= this->gas_density( p_pressureP );			
                real_type rhoG_W		= this->gas_density( p_pressureW );

                real_type rhoW_P		= this->water_density( p_pressureP );			
                real_type rhoW_W		= this->water_density( p_pressureW );

                real_type rhoO_P_old	= this->oil_density( m_pressure_old[ p_node ] );	
                real_type rhoO_P		= this->oil_density( p_pressureP );			
                real_type rhoO_W		= this->oil_density( p_pressureW );

                real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );			
                real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
                
                real_type mixture_inlet;

                real_type Qoil = m_oil_flow[ p_node ]->get_current_value();
                real_type Qwater = m_water_flow[ p_node ]->get_current_value();
                real_type Qgas = m_gas_flow[ p_node ]->get_current_value();

                if(m_mass_flux){   
                    mixture_inlet = Qoil + Qwater + Qgas;
                }                      
                else
                {
				    real_type rhoGas_P	     = this->gas_density	( p_pressureP );
				    real_type rhoWater_P	 = this->water_density	( p_pressureP );
				    real_type rhoOil_P	     = this->oil_density	( p_pressureP );
				    mixture_inlet  = rhoOil_P*Qoil + rhoWater_P*Qwater + rhoGas_P*Qgas;
                }  				

				


                real_type mod_Vow_w	= this->mod_v_drift_flux_ow( 
                    p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                    0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                    );


                real_type mod_Vgj_w	= this->mod_v_drift_flux( 
                    p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                    0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                    );

                real_type gas_vol_frac_w   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
                real_type oil_vol_frac_w   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracW);
                real_type water_vol_frac_w = 0.5*(water_vol_fracP + water_vol_fracW);
                real_type rho_g_w = 0.5*(rhoG_W + rhoG_P);
                real_type rho_o_w = 0.5*(rhoO_W + rhoO_P);
                real_type rho_w_w = 0.5*(rhoW_W + rhoW_P);
                real_type rho_l_w = 0.5*(rhoL_W + rhoL_P);
                real_type rho_w = 0.5*(rho_W + rho_P);

                real_type gas_velocity_w   = p_velocityW + rho_l_w/rho_w*mod_Vgj_w;
                real_type liquid_velocity_w = p_velocityW - gas_vol_frac_w/(1 - gas_vol_frac_w + 1.0e-20)*rho_g_w/rho_w*mod_Vgj_w;
                real_type oil_velocity_w   = liquid_velocity_w + rho_w_w/rho_l_w*mod_Vow_w;               
                real_type water_velocity_w = liquid_velocity_w - (oil_vol_frac_w/(water_vol_frac_w + 1.0e-20))*(rho_o_w/rho_l_w)*mod_Vow_w;               
                
              

                real_type ksi_w_oil     = this->ksi( oil_velocity_w   );
                real_type ksi_w_gas     = this->ksi( gas_velocity_w   );
                real_type ksi_w_water   = this->ksi( water_velocity_w );

                real_type m_w_oil   = oil_velocity_w  *( (0.5+ksi_w_oil  )*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w_oil  )*rhoO_P*p_oil_vol_fracP );
                real_type m_w_water = water_velocity_w*( (0.5+ksi_w_water)*rhoW_W*water_vol_fracW + (0.5-ksi_w_water)*rhoW_P*water_vol_fracP );
                real_type m_w_gas   = gas_velocity_w  *( (0.5+ksi_w_gas  )*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w_gas  )*rhoG_P*p_gas_vol_fracP );

                return (rho_P-rho_P_old)*dV/dt() - mixture_inlet 
                    +	area()*( p_velocityP*rho_P - (m_w_oil + m_w_water + m_w_gas) );

                //real_type ksi_w = this->ksi( p_velocityW );	
			/*return (rho_P-rho_P_old)*dV/dt() - mixture_inlet 
				  +	area()*( p_velocityP*rho_P 
						   - p_velocityW*( (0.5+ksi_w)*rho_W + (0.5-ksi_w)*rho_P ) );*/	
			
			}			
		
		case 'C':
			{
				
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.5*this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
				real_type dS  = dSw + dSe;
				real_type dV = this->Volume( dS );

				real_type water_vol_fracOld = 1.0 - (m_oil_vol_frac_old[ p_node ] + m_gas_vol_frac_old[ p_node ]);
				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
				real_type water_vol_fracE	= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);

				real_type rho_P_old = this->mean_density(m_oil_vol_frac_old[ p_node ], water_vol_fracOld, m_gas_vol_frac_old[ p_node ], m_pressure_old[ p_node ]);
                


				real_type rho_P = this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);
				real_type rho_E = this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
				real_type rho_W = this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);

                real_type rhoG_P		= this->gas_density( p_pressureP );
                real_type rhoG_E		= this->gas_density( p_pressureE );
                real_type rhoG_W		= this->gas_density( p_pressureW );

                real_type rhoW_P		= this->water_density( p_pressureP );
                real_type rhoW_E		= this->water_density( p_pressureE );
                real_type rhoW_W		= this->water_density( p_pressureW );

                real_type rhoO_P_old	= this->oil_density( m_pressure_old[ p_node ] );	
                real_type rhoO_P		= this->oil_density( p_pressureP );
                real_type rhoO_E		= this->oil_density( p_pressureE );
                real_type rhoO_W		= this->oil_density( p_pressureW );

                real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
                real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
                real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
                
                real_type mixture_inlet;

                real_type Qoil = m_oil_flow[ p_node ]->get_current_value();
                real_type Qwater = m_water_flow[ p_node ]->get_current_value();
                real_type Qgas = m_gas_flow[ p_node ]->get_current_value();

                if(m_mass_flux){   
                    mixture_inlet = Qoil + Qwater + Qgas;
                }                      
                else
                {
                    real_type rhoGas_P	     = this->gas_density	( p_pressureP );
                    real_type rhoWater_P	 = this->water_density	( p_pressureP );
                    real_type rhoOil_P	     = this->oil_density	( p_pressureP );
                    mixture_inlet  = rhoOil_P*Qoil + rhoWater_P*Qwater + rhoGas_P*Qgas;
                }  


                real_type mod_Vow_e	= this->mod_v_drift_flux_ow( 
                    p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                    0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                    );

                real_type mod_Vow_w	= this->mod_v_drift_flux_ow( 
                    p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                    0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                    );

                real_type mod_Vgj_e		= this->mod_v_drift_flux( 
                    p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                    0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                    );

                real_type mod_Vgj_w		= this->mod_v_drift_flux( 
                    p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                    0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                    );

                real_type gas_vol_frac_w   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
                real_type oil_vol_frac_w   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracW);
                real_type water_vol_frac_w = 0.5*(water_vol_fracP + water_vol_fracW);
                real_type rho_g_w = 0.5*(rhoG_W + rhoG_P);
                real_type rho_o_w = 0.5*(rhoO_W + rhoO_P);
                real_type rho_w_w = 0.5*(rhoW_W + rhoW_P);
                real_type rho_l_w = 0.5*(rhoL_W + rhoL_P);
                real_type rho_w = 0.5*(rho_W + rho_P);

                real_type gas_velocity_w   = p_velocityW + rho_l_w/rho_w*mod_Vgj_w;
                real_type liquid_velocity_w = p_velocityW - gas_vol_frac_w/(1 - gas_vol_frac_w + 1.0e-20)*rho_g_w/rho_w*mod_Vgj_w;	 
                real_type water_velocity_w = liquid_velocity_w - (oil_vol_frac_w/(water_vol_frac_w + 1.0e-20))*(rho_o_w/rho_l_w)*mod_Vow_w;               
                real_type oil_velocity_w = liquid_velocity_w + rho_w_w/rho_l_w*mod_Vow_w;  

                
                
                real_type gas_vol_frac_e   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE);
                real_type oil_vol_frac_e   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracE);
                real_type water_vol_frac_e = 0.5*(water_vol_fracP + water_vol_fracE);
                real_type rho_g_e = 0.5*(rhoG_E + rhoG_P);
                real_type rho_o_e = 0.5*(rhoO_E + rhoO_P);
                real_type rho_w_e = 0.5*(rhoW_E + rhoW_P);
                real_type rho_l_e = 0.5*(rhoL_E + rhoL_P);
                real_type rho_e = 0.5*(rho_E + rho_P);

                real_type gas_velocity_e = p_velocityP + rho_l_e/rho_e*mod_Vgj_e;
                real_type liquid_velocity_e = p_velocityP - gas_vol_frac_e/(1 - gas_vol_frac_e + 1.0e-20)*rho_g_e/rho_e*mod_Vgj_e;	 
                real_type water_velocity_e = liquid_velocity_e - (oil_vol_frac_e/(water_vol_frac_e + 1.0e-20))*(rho_o_e/rho_l_e)*mod_Vow_e;               
                real_type oil_velocity_e = liquid_velocity_e + rho_w_e/rho_l_e*mod_Vow_e;

              

                real_type ksi_e_oil     = this->ksi( oil_velocity_e   );
                real_type ksi_e_gas     = this->ksi( gas_velocity_e   );
                real_type ksi_e_water   = this->ksi( water_velocity_e );

                real_type ksi_w_oil     = this->ksi( oil_velocity_w   );
                real_type ksi_w_gas     = this->ksi( gas_velocity_w   );
                real_type ksi_w_water   = this->ksi( water_velocity_w );	

                real_type m_w_oil   = oil_velocity_w  *( (0.5+ksi_w_oil  )*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w_oil  )*rhoO_P*p_oil_vol_fracP );
                real_type m_w_water = water_velocity_w*( (0.5+ksi_w_water)*rhoW_W*water_vol_fracW + (0.5-ksi_w_water)*rhoW_P*water_vol_fracP );
                real_type m_w_gas   = gas_velocity_w  *( (0.5+ksi_w_gas  )*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w_gas  )*rhoG_P*p_gas_vol_fracP );

                real_type m_e_oil   = oil_velocity_e  *( (0.5+ksi_e_oil  )*rhoO_P*p_oil_vol_fracP + (0.5-ksi_e_oil  )*rhoO_E*p_oil_vol_fracE );
                real_type m_e_water = water_velocity_e*( (0.5+ksi_e_water)*rhoW_P*water_vol_fracP + (0.5-ksi_e_water)*rhoW_E*water_vol_fracE );
                real_type m_e_gas   = gas_velocity_e  *( (0.5+ksi_e_gas  )*rhoG_P*p_gas_vol_fracP + (0.5-ksi_e_gas  )*rhoG_E*p_gas_vol_fracE );


                return (rho_P-rho_P_old)*dV/dt() - mixture_inlet 
                    + area()*( m_e_oil+m_e_water+m_e_gas - ( m_w_oil+m_w_water+m_w_gas ) );


				/*real_type ksi_e = this->ksi( p_velocityP );
				real_type ksi_w = this->ksi( p_velocityW );			

			return (rho_P-rho_P_old)*dV/dt() - mixture_inlet 
				  + area()*( p_velocityP*( (0.5+ksi_e)*rho_P + (0.5-ksi_e)*rho_E )
						   - p_velocityW*( (0.5+ksi_w)*rho_W + (0.5-ksi_w)*rho_P ) );*/
						
			}
		default:
			return 0.;	
		}
		
	}

	real_type DriftFluxWell::R_g(
								 real_type p_pressureW,
								 real_type p_pressureP,
								 real_type p_pressureE,
								 real_type p_gas_vol_fracW,
								 real_type p_gas_vol_fracP,
								 real_type p_gas_vol_fracE,
								 real_type p_oil_vol_fracW,
								 real_type p_oil_vol_fracP,
								 real_type p_oil_vol_fracE,
								 real_type p_velocityW,
								 real_type p_velocityP,
								 uint_type p_node,
								 string_type	   position = 'C'
								 )
	{
		switch( position )
		{

		case 'L':
			{
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.;
				real_type dS  = dSw + dSe;
				real_type dV = this->Volume( dS );

				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);

				real_type rho_P		= this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);			
				real_type rho_W		= this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);

				real_type rhoG_P_old	= this->gas_density( m_pressure_old[ p_node ] );	
				real_type rhoG_P		= this->gas_density( p_pressureP );			
				real_type rhoG_W		= this->gas_density( p_pressureW );

				real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );			
				real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
                
                real_type gas_inlet;

                real_type Qgas = m_gas_flow[ p_node ]->get_current_value();
                if(m_mass_flux){
                    gas_inlet = Qgas;
                }                      
                else
                {
                    gas_inlet = rhoG_P*Qgas;
                }



				real_type mod_Vgj_w	= this->mod_v_drift_flux( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);

                real_type gas_velocity_w = p_velocityW + 0.5*(rhoL_W + rhoL_P)/(0.5*(rho_W + rho_P))*mod_Vgj_w;
				real_type ksi_w = this->ksi( gas_velocity_w );	

                return (p_gas_vol_fracP*rhoG_P - m_gas_vol_frac_old[ p_node ]*rhoG_P_old)*dV/dt() - gas_inlet				 
				 + p_velocityP*area()*rhoG_P*p_gas_vol_fracP
				 - gas_velocity_w*area()*( (0.5+ksi_w)*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*p_gas_vol_fracP );				 

			//return (p_gas_vol_fracP*rhoG_P - m_gas_vol_frac_old[ p_node ]*rhoG_P_old)*dV/dt() - gas_inlet				 
			//	 + p_velocityP*area()*rhoG_P*p_gas_vol_fracP
			//	 - p_velocityW*area()*( (0.5+ksi_w)*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*p_gas_vol_fracP )				 
			//	 //+ mod_Vgj_e*area()*( rhoG_P*rhoL_P/rho_P * p_gas_vol_fracP )
			//	 - mod_Vgj_w*area()*( (0.5+ksi_w)*rhoG_W*rhoL_W/rho_W * p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*rhoL_P/rho_P * p_gas_vol_fracP );			
			
			}
			

		case 'C':
			{
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.5*this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
				real_type dS  = dSw + dSe;
				real_type dV  = this->Volume( dS );			

				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
				real_type water_vol_fracE	= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);

				real_type rho_P		= this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);	
				real_type rho_E		= this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
				real_type rho_W		= this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);

				real_type rhoG_P_old = this->gas_density( m_pressure_old[ p_node ] );	
				real_type rhoG_P		= this->gas_density( p_pressureP );
				real_type rhoG_E		= this->gas_density( p_pressureE );
				real_type rhoG_W		= this->gas_density( p_pressureW );

				real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
				real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
				real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );

                real_type gas_inlet;
                real_type Qgas = m_gas_flow[ p_node ]->get_current_value();
                if(m_mass_flux){
                    gas_inlet = Qgas;
                }                      
                else
                {
                    gas_inlet = rhoG_P*Qgas;
                }

				real_type mod_Vgj_e		= this->mod_v_drift_flux( 
					p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
					0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
					);

				real_type mod_Vgj_w		= this->mod_v_drift_flux( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);

                real_type gas_velocity_w = p_velocityW + 0.5*(rhoL_W + rhoL_P)/(0.5*(rho_W + rho_P))*mod_Vgj_w;
                real_type gas_velocity_e = p_velocityP + 0.5*(rhoL_P + rhoL_E)/(0.5*(rho_P + rho_E))*mod_Vgj_e;

				real_type ksi_e = this->ksi( gas_velocity_e );
				real_type ksi_w = this->ksi( gas_velocity_w );
                return (p_gas_vol_fracP*rhoG_P - m_gas_vol_frac_old[ p_node ]*rhoG_P_old)*dV/dt() - gas_inlet
                    + gas_velocity_e*area()*( (0.5+ksi_e)*rhoG_P*p_gas_vol_fracP + (0.5-ksi_e)*rhoG_E*p_gas_vol_fracE )
                    - gas_velocity_w*area()*( (0.5+ksi_w)*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*p_gas_vol_fracP );
                   
			/*return (p_gas_vol_fracP*rhoG_P - m_gas_vol_frac_old[ p_node ]*rhoG_P_old)*dV/dt() - gas_inlet
				 + p_velocityP*area()*( (0.5+ksi_e)*rhoG_P*p_gas_vol_fracP + (0.5-ksi_e)*rhoG_E*p_gas_vol_fracE )
				 - p_velocityW*area()*( (0.5+ksi_w)*rhoG_W*p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*p_gas_vol_fracP )
				 + mod_Vgj_e*area()*( (0.5+ksi_e)*rhoG_P*rhoL_P/rho_P * p_gas_vol_fracP + (0.5-ksi_e)*rhoG_E*rhoL_E/rho_E * p_gas_vol_fracE )
				 - mod_Vgj_w*area()*( (0.5+ksi_w)*rhoG_W*rhoL_W/rho_W * p_gas_vol_fracW + (0.5-ksi_w)*rhoG_P*rhoL_P/rho_P * p_gas_vol_fracP );*/
			}

		default:
			return 0.;
			
		}
	}

	real_type DriftFluxWell::R_o(
								 real_type p_pressureW,
								 real_type p_pressureP,
								 real_type p_pressureE,
								 real_type p_gas_vol_fracW,
								 real_type p_gas_vol_fracP,
								 real_type p_gas_vol_fracE,
								 real_type p_oil_vol_fracW,
								 real_type p_oil_vol_fracP,
								 real_type p_oil_vol_fracE,
								 real_type p_velocityW,
								 real_type p_velocityP,
								 uint_type p_node,
								 string_type	   position = 'C'
								 )
	{
		switch( position )
		{

		case 'L':
			{
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.;
				real_type dS  = dSw + dSe;
				real_type dV = this->Volume( dS );

				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP); 
                
				real_type rho_P		= this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);	                
				real_type rho_W		= this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);

									
				real_type rhoG_P		= this->gas_density( p_pressureP );			
				real_type rhoG_W		= this->gas_density( p_pressureW );

				real_type rhoW_P		= this->water_density( p_pressureP );			
				real_type rhoW_W		= this->water_density( p_pressureW );

				real_type rhoO_P_old	= this->oil_density( m_pressure_old[ p_node ] );	
				real_type rhoO_P		= this->oil_density( p_pressureP );			
				real_type rhoO_W		= this->oil_density( p_pressureW );

				real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );			
				real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
                
                real_type oil_inlet;
                real_type Qoil = m_oil_flow[ p_node ]->get_current_value();
                if(m_mass_flux){
                    oil_inlet = Qoil;
                }                      
                else
                {
                    oil_inlet = rhoO_P*Qoil;
                }
				

				real_type mod_Vow_w	= this->mod_v_drift_flux_ow( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);
               

				real_type mod_Vgj_w	= this->mod_v_drift_flux( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);

				//real_type mod_Vgj_e	= this->mod_v_drift_flux(p_velocityP, p_gas_vol_fracP, p_oil_vol_fracP,	water_vol_fracP, p_pressureP);
                real_type gas_vol_frac_w = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
                real_type rho_g_w = 0.5*(rhoG_W + rhoG_P);
                real_type rho_w_w = 0.5*(rhoW_W + rhoW_P);
                real_type rho_l_w = 0.5*(rhoL_W + rhoL_P);
                real_type rho_w = 0.5*(rho_W + rho_P);
                real_type liquid_velocity_w = p_velocityW - gas_vol_frac_w/(1 - gas_vol_frac_w + 1.0e-20)*rho_g_w/rho_w*mod_Vgj_w;	 
                
                real_type oil_velocity_w = liquid_velocity_w + rho_w_w/rho_l_w*mod_Vow_w;               
				real_type ksi_w = this->ksi( oil_velocity_w );

                return (p_oil_vol_fracP*rhoO_P - m_oil_vol_frac_old[ p_node ]*rhoO_P_old)*dV/dt() - oil_inlet				 
                    + p_velocityP*area()*rhoO_P*p_oil_vol_fracP 
                    - oil_velocity_w*area()*( (0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w)*rhoO_P*p_oil_vol_fracP );				 
                                
			//return (p_oil_vol_fracP*rhoO_P - m_oil_vol_frac_old[ p_node ]*rhoO_P_old)*dV/dt() - oil_inlet				 
			//	 + p_velocityP*area()*rhoO_P*p_oil_vol_fracP 
   //              - p_velocityW*area()*( (0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w)*rhoO_P*p_oil_vol_fracP )				 
			//	 - mod_Vow_w*area()*( rhoW_W/rhoL_W*(0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + rhoW_P/rhoL_P*(0.5-ksi_w)*rhoO_P*p_oil_vol_fracP )	
			//	 //- mod_Vgj_e*area()*( rhoG_P/rho_P*p_gas_vol_fracP/(1.0-p_gas_vol_fracP + 1.0e-20)*rhoO_P*p_oil_vol_fracP )
			//	 + mod_Vgj_w*area()*( rhoG_W/rho_W*p_gas_vol_fracW/(1.0-p_gas_vol_fracW + 1.0e-20)*(0.5+ksi_w)*rhoO_W*p_oil_vol_fracW 
			//						+ rhoG_P/rho_P*p_gas_vol_fracP/(1.0-p_gas_vol_fracP + 1.0e-20)*(0.5-ksi_w)*rhoO_P*p_oil_vol_fracP );
			
			}
			

		case 'C':
			{
				real_type dSw = 0.5*this->segment_length( m_coordinates[ p_node-1 ], m_coordinates[ p_node ] );
				real_type dSe = 0.5*this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
				real_type dS  = dSw + dSe;
				real_type dV  = this->Volume( dS );	

				real_type water_vol_fracW	= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
				real_type water_vol_fracP	= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
				real_type water_vol_fracE	= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);

				real_type rho_P		= this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);	
				real_type rho_E		= this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
				real_type rho_W		= this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW, p_pressureW);
								
				real_type rhoG_P		= this->gas_density( p_pressureP );
				real_type rhoG_E		= this->gas_density( p_pressureE );
				real_type rhoG_W		= this->gas_density( p_pressureW );

				real_type rhoW_P		= this->water_density( p_pressureP );
				real_type rhoW_E		= this->water_density( p_pressureE );
				real_type rhoW_W		= this->water_density( p_pressureW );

				real_type rhoO_P_old	= this->oil_density( m_pressure_old[ p_node ] );	
				real_type rhoO_P		= this->oil_density( p_pressureP );
				real_type rhoO_E		= this->oil_density( p_pressureE );
				real_type rhoO_W		= this->oil_density( p_pressureW );

				real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
				real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
				real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );

                real_type oil_inlet;
                real_type Qoil = m_oil_flow[ p_node ]->get_current_value();
                if(m_mass_flux){
                    oil_inlet = Qoil;
                }                      
                else
                {
                    oil_inlet = rhoO_P*Qoil;
                }

				real_type mod_Vow_e	= this->mod_v_drift_flux_ow( 
					p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
					0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
					);

				real_type mod_Vow_w	= this->mod_v_drift_flux_ow( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);

				real_type mod_Vgj_e		= this->mod_v_drift_flux( 
					p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
					0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
					);

				real_type mod_Vgj_w		= this->mod_v_drift_flux( 
					p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
					0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
					);

                real_type gas_vol_frac_w = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
                real_type rho_g_w = 0.5*(rhoG_W + rhoG_P);
                real_type rho_w_w = 0.5*(rhoW_W + rhoW_P);
                real_type rho_l_w = 0.5*(rhoL_W + rhoL_P);
                real_type rho_w = 0.5*(rho_W + rho_P);
                real_type liquid_velocity_w = p_velocityW - gas_vol_frac_w/(1 - gas_vol_frac_w + 1.0e-20)*rho_g_w/rho_w*mod_Vgj_w;	 

                real_type oil_velocity_w = liquid_velocity_w + rho_w_w/rho_l_w*mod_Vow_w;  

                real_type gas_vol_frac_e = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE);
                real_type rho_g_e = 0.5*(rhoG_E + rhoG_P);
                real_type rho_w_e = 0.5*(rhoW_E + rhoW_P);
                real_type rho_l_e = 0.5*(rhoL_E + rhoL_P);
                real_type rho_e = 0.5*(rho_E + rho_P);
                real_type liquid_velocity_e = p_velocityP - gas_vol_frac_e/(1 - gas_vol_frac_e + 1.0e-20)*rho_g_e/rho_e*mod_Vgj_e;	 

                real_type oil_velocity_e = liquid_velocity_e + rho_w_e/rho_l_e*mod_Vow_e;  
                
				real_type ksi_e = this->ksi( oil_velocity_e );
				real_type ksi_w = this->ksi( oil_velocity_w );

                return (p_oil_vol_fracP*rhoO_P - m_oil_vol_frac_old[ p_node ]*rhoO_P_old)*dV/dt() - oil_inlet
                    + oil_velocity_e*area()*( (0.5+ksi_e)*rhoO_P*p_oil_vol_fracP + (0.5-ksi_e)*rhoO_E*p_oil_vol_fracE )
                    - oil_velocity_w*area()*( (0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w)*rhoO_P*p_oil_vol_fracP );
                    
               
			/*return (p_oil_vol_fracP*rhoO_P - m_oil_vol_frac_old[ p_node ]*rhoO_P_old)*dV/dt() - oil_inlet
				 + p_velocityP*area()*( (0.5+ksi_e)*rhoO_P*p_oil_vol_fracP + (0.5-ksi_e)*rhoO_E*p_oil_vol_fracE )
				 - p_velocityW*area()*( (0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + (0.5-ksi_w)*rhoO_P*p_oil_vol_fracP )
				 + mod_Vow_e*area()*( rhoW_P/rhoL_P*(0.5+ksi_e)*rhoO_P*p_oil_vol_fracP + rhoW_E/rhoL_E*(0.5-ksi_e)*rhoO_E*p_oil_vol_fracE )
				 - mod_Vow_w*area()*( rhoW_W/rhoL_W*(0.5+ksi_w)*rhoO_W*p_oil_vol_fracW + rhoW_P/rhoL_P*(0.5-ksi_w)*rhoO_P*p_oil_vol_fracP )
				 - mod_Vgj_e*area()*( rhoG_P/rho_P*p_gas_vol_fracP/(1.0-p_gas_vol_fracP + 1.0e-20)*(0.5+ksi_e)*rhoO_P*p_oil_vol_fracP 
									+ rhoG_E/rho_E*p_gas_vol_fracE/(1.0-p_gas_vol_fracE + 1.0e-20)*(0.5-ksi_e)*rhoO_E*p_oil_vol_fracE )
				 + mod_Vgj_w*area()*( rhoG_W/rho_W*p_gas_vol_fracW/(1.0-p_gas_vol_fracW + 1.0e-20)*(0.5+ksi_w)*rhoO_W*p_oil_vol_fracW 
									+ rhoG_P/rho_P*p_gas_vol_fracP/(1.0-p_gas_vol_fracP + 1.0e-20)*(0.5-ksi_w)*rhoO_P*p_oil_vol_fracP );*/
			}

		default:
			return 0.;
			
		}
	}

	



	real_type DriftFluxWell::R_v(
                                 real_type p_pressureW,
								 real_type p_pressureP,
								 real_type p_pressureE,
                                 real_type p_pressureEE,
                                 real_type p_gas_vol_fracW,
								 real_type p_gas_vol_fracP,
								 real_type p_gas_vol_fracE,
                                 real_type p_gas_vol_fracEE,
                                 real_type p_oil_vol_fracW,
								 real_type p_oil_vol_fracP,
								 real_type p_oil_vol_fracE,
                                 real_type p_oil_vol_fracEE,
								 real_type p_velocityW,
								 real_type p_velocityP,
								 real_type p_velocityE,
								 uint_type p_node,
								 string_type	   position = 'C'
								 )
	{		
		switch( position )
		{
		case 'F':
			{
			real_type dS  = this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
			real_type dSe = this->segment_length( m_coordinates[ p_node+1 ], m_coordinates[ p_node+2 ] );
			real_type dSw = 0.;
			real_type dV = this->Volume( dS );

			
			real_type water_vol_fracP_old	= 1.0 - (m_oil_vol_frac_old[ p_node ]   + m_gas_vol_frac_old[ p_node ]	);
			real_type water_vol_fracE_old	= 1.0 - (m_oil_vol_frac_old[ p_node+1 ] + m_gas_vol_frac_old[ p_node+1 ]);
            real_type water_vol_fracW		= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
            real_type water_vol_fracP		= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
            real_type water_vol_fracE		= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);	
            real_type water_vol_fracEE		= 1.0 - (p_gas_vol_fracEE + p_oil_vol_fracEE);

			real_type rho_P_old = this->mean_density(m_oil_vol_frac_old[ p_node ], water_vol_fracP_old, m_gas_vol_frac_old[ p_node ], m_pressure_old[ p_node ]);
			real_type rho_E_old = this->mean_density(m_oil_vol_frac_old[ p_node+1 ], water_vol_fracE_old, m_gas_vol_frac_old[ p_node+1 ], m_pressure_old[ p_node+1 ]);

			vector_type S( 3 );
			S[ 0 ] = m_coordinates[ p_node ].getX() - m_coordinates[ p_node+1 ].getX();
			S[ 1 ] = m_coordinates[ p_node ].getY() - m_coordinates[ p_node+1 ].getY();
			S[ 2 ] = m_coordinates[ p_node ].getZ() - m_coordinates[ p_node+1 ].getZ();

			real_type d_e = dS/(dS+dSe);
			real_type d_w = dS/(dS+dSw);
						
			real_type angle = get_inclination() - PI/2;//PI/2 - 0*acos( dot( m_gravity, S )/(norm(m_gravity)*norm(S)) );
			
            real_type rho_W   = this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW , p_pressureW );
            real_type rho_P   = this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);
            real_type rho_E   = this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
            real_type rho_EE  = this->mean_density(p_oil_vol_fracEE, water_vol_fracEE, p_gas_vol_fracEE , p_pressureEE);

            real_type rhoG_W		= this->gas_density( p_pressureW );
            real_type rhoG_P		= this->gas_density( p_pressureP );
            real_type rhoG_E		= this->gas_density( p_pressureE );
            real_type rhoG_EE		= this->gas_density( p_pressureEE );

            real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
            real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
            real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
            real_type rhoL_EE		= this->liquid_density( p_oil_vol_fracEE, water_vol_fracEE, p_pressureEE );

            real_type rhoW_W		= this->water_density( p_pressureW );
            real_type rhoW_P		= this->water_density( p_pressureP );
            real_type rhoW_E		= this->water_density( p_pressureE );
            real_type rhoW_EE		= this->water_density( p_pressureEE );

            real_type rhoO_W		= this->oil_density( p_pressureW );
            real_type rhoO_P		= this->oil_density( p_pressureP );
            real_type rhoO_E		= this->oil_density( p_pressureE );
            real_type rhoO_EE		= this->oil_density( p_pressureEE ); 
			
			real_type mean_pressure = 0.5*(p_pressureP + p_pressureE);
			real_type viscosity = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE)*gas_viscosity	 ( mean_pressure ) 
								+ 0.5*(p_oil_vol_fracP + p_oil_vol_fracE)*oil_viscosity	 ( mean_pressure ) 
								+ 0.5*(water_vol_fracP + water_vol_fracE)*water_viscosity( mean_pressure );


            real_type mod_Vow_W	= this->mod_v_drift_flux_ow( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vow_P	= this->mod_v_drift_flux_ow( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vow_E	= this->mod_v_drift_flux_ow( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );

            real_type mod_Vgj_W	= this->mod_v_drift_flux( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vgj_P	= this->mod_v_drift_flux( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vgj_E	= this->mod_v_drift_flux( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );


            //real_type Vc  = p_velocityP;			
            //real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2*m_radius/viscosity);
            //real_type f_P = this->friction_factor( Re );

            // OTHER Vc = j
            real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
            real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2.0*m_radius/viscosity);
            real_type f_P = this->friction_factor( Re );


            real_type gas_vol_frac_W   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
            real_type oil_vol_frac_W   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracW);
            real_type water_vol_frac_W = 0.5*(water_vol_fracP + water_vol_fracW);
            real_type rho_g_W = 0.5*(rhoG_W + rhoG_P);
            real_type rho_o_W = 0.5*(rhoO_W + rhoO_P);
            real_type rho_w_W = 0.5*(rhoW_W + rhoW_P);
            real_type rho_l_W = 0.5*(rhoL_W + rhoL_P);
            real_type rho_m_W = 0.5*(rho_W + rho_P);

            real_type gas_velocity_W    = p_velocityW + rho_l_W/rho_g_W*mod_Vgj_W;
            real_type liquid_velocity_W = p_velocityW - gas_vol_frac_W/(1 - gas_vol_frac_W + 1.0e-20)*rho_g_W/rho_m_W*mod_Vgj_W;	 
            real_type water_velocity_W  = liquid_velocity_W - (oil_vol_frac_W/(water_vol_frac_W + 1.0e-20))*(rho_o_W/rho_l_W)*mod_Vow_W;               
            real_type oil_velocity_W    = liquid_velocity_W + rho_w_W/rho_l_W*mod_Vow_W;  



            real_type gas_vol_frac_P   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE);
            real_type oil_vol_frac_P   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracE);
            real_type water_vol_frac_P = 0.5*(water_vol_fracP + water_vol_fracE);
            real_type rho_g_P = 0.5*(rhoG_E + rhoG_P);
            real_type rho_o_P = 0.5*(rhoO_E + rhoO_P);
            real_type rho_w_P = 0.5*(rhoW_E + rhoW_P);
            real_type rho_l_P = 0.5*(rhoL_E + rhoL_P);
            real_type rho_m_P = 0.5*(rho_E + rho_P);

            real_type gas_velocity_P    = p_velocityP + rho_l_P/rho_g_P*mod_Vgj_P;
            real_type liquid_velocity_P = p_velocityP - gas_vol_frac_P/(1 - gas_vol_frac_P + 1.0e-20)*rho_g_P/rho_m_P*mod_Vgj_P;	 
            real_type water_velocity_P  = liquid_velocity_P - (oil_vol_frac_P/(water_vol_frac_P + 1.0e-20))*(rho_o_P/rho_l_P)*mod_Vow_P;               
            real_type oil_velocity_P    = liquid_velocity_P + rho_w_P/rho_l_P*mod_Vow_P;

            real_type gas_vol_frac_E   = 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE);
            real_type oil_vol_frac_E   = 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE);
            real_type water_vol_frac_E = 0.5*(water_vol_fracE + water_vol_fracEE);
            real_type rho_g_E = 0.5*(rhoG_E + rhoG_EE);
            real_type rho_o_E = 0.5*(rhoO_E + rhoO_EE);
            real_type rho_w_E = 0.5*(rhoW_E + rhoW_EE);
            real_type rho_l_E = 0.5*(rhoL_E + rhoL_EE);
            real_type rho_m_E = 0.5*(rho_E + rho_EE);

            real_type gas_velocity_E    = p_velocityE + rho_l_E/rho_g_E*mod_Vgj_E;
            real_type liquid_velocity_E = p_velocityE - gas_vol_frac_E/(1 - gas_vol_frac_E + 1.0e-20)*rho_g_E/rho_m_E*mod_Vgj_E;	 
            real_type water_velocity_E  = liquid_velocity_E - (oil_vol_frac_E/(water_vol_frac_E + 1.0e-20))*(rho_o_E/rho_l_E)*mod_Vow_E;               
            real_type oil_velocity_E    = liquid_velocity_E + rho_w_E/rho_l_E*mod_Vow_E;



            real_type ksi_W_oil     = this->ksi( oil_velocity_W   );
            real_type ksi_W_gas     = this->ksi( gas_velocity_W   );
            real_type ksi_W_water   = this->ksi( water_velocity_W );

            real_type ksi_P_oil     = this->ksi( oil_velocity_P   );
            real_type ksi_P_gas     = this->ksi( gas_velocity_P   );
            real_type ksi_P_water   = this->ksi( water_velocity_P );

            real_type ksi_E_oil     = this->ksi( oil_velocity_E   );
            real_type ksi_E_gas     = this->ksi( gas_velocity_E   );
            real_type ksi_E_water   = this->ksi( water_velocity_E );

            real_type m_W_oil   = oil_velocity_W  *( (0.5+ksi_W_oil  )*rhoO_W*p_oil_vol_fracW + (0.5-ksi_W_oil  )*rhoO_P*p_oil_vol_fracP );
            real_type m_W_water = water_velocity_W*( (0.5+ksi_W_water)*rhoW_W*water_vol_fracW + (0.5-ksi_W_water)*rhoW_P*water_vol_fracP );
            real_type m_W_gas   = gas_velocity_W  *( (0.5+ksi_W_gas  )*rhoG_W*p_gas_vol_fracW + (0.5-ksi_W_gas  )*rhoG_P*p_gas_vol_fracP );

            real_type m_P_oil   = oil_velocity_P  *( (0.5+ksi_P_oil  )*rhoO_P*p_oil_vol_fracP + (0.5-ksi_P_oil  )*rhoO_E*p_oil_vol_fracE );
            real_type m_P_water = water_velocity_P*( (0.5+ksi_P_water)*rhoW_P*water_vol_fracP + (0.5-ksi_P_water)*rhoW_E*water_vol_fracE );
            real_type m_P_gas   = gas_velocity_P  *( (0.5+ksi_P_gas  )*rhoG_P*p_gas_vol_fracP + (0.5-ksi_P_gas  )*rhoG_E*p_gas_vol_fracE );

            real_type m_E_oil   = oil_velocity_E  *( (0.5+ksi_E_oil  )*rhoO_P*p_oil_vol_fracE + (0.5-ksi_E_oil  )*rhoO_E*p_oil_vol_fracEE );
            real_type m_E_water = water_velocity_E*( (0.5+ksi_E_water)*rhoW_P*water_vol_fracE + (0.5-ksi_E_water)*rhoW_E*water_vol_fracEE );
            real_type m_E_gas   = gas_velocity_E  *( (0.5+ksi_E_gas  )*rhoG_P*p_gas_vol_fracE + (0.5-ksi_E_gas  )*rhoG_E*p_gas_vol_fracEE );

            real_type ksi_e_oil = this->ksi( (1-d_e)*oil_velocity_P + d_e*oil_velocity_E );
            real_type ksi_w_oil = this->ksi( (1-d_w)*oil_velocity_P + d_w*oil_velocity_W );
            real_type ksi_e_water = this->ksi( (1-d_e)*water_velocity_P + d_e*water_velocity_E );
            real_type ksi_w_water = this->ksi( (1-d_w)*water_velocity_P + d_w*water_velocity_W );
            real_type ksi_e_gas = this->ksi( (1-d_e)*gas_velocity_P + d_e*gas_velocity_E );
            real_type ksi_w_gas = this->ksi( (1-d_w)*gas_velocity_P + d_w*gas_velocity_W );

            real_type m_e = 0.5*(m_E_oil  +m_P_oil)  *( (0.5+ksi_e_oil  )*oil_velocity_P   + (0.5-ksi_e_oil  )*oil_velocity_E )
                          + 0.5*(m_E_water+m_P_water)*( (0.5+ksi_e_water)*water_velocity_P + (0.5-ksi_e_water)*water_velocity_E )
                          + 0.5*(m_E_gas  +m_P_gas)  *( (0.5+ksi_e_gas  )*gas_velocity_P   + (0.5-ksi_e_gas  )*gas_velocity_E );

            /*real_type m_w = 0.5*(m_W_oil  +m_P_oil)  *( (0.5+ksi_w_oil  )*oil_velocity_W   + (0.5-ksi_w_oil  )*oil_velocity_P )
                          + 0.5*(m_W_water+m_P_water)*( (0.5+ksi_w_water)*water_velocity_W + (0.5-ksi_w_water)*water_velocity_P )
                          + 0.5*(m_W_gas  +m_P_gas)  *( (0.5+ksi_w_gas  )*gas_velocity_W   + (0.5-ksi_w_gas  )*gas_velocity_P ); */ 
            real_type m_w = m_P_oil*oil_velocity_P +m_P_water*water_velocity_P + m_P_gas*gas_velocity_P ; 

            real_type m_t =  oil_velocity_P  *(rhoO_P*p_oil_vol_fracP + rhoO_E*p_oil_vol_fracE)
                          +  water_velocity_P*(rhoW_P*water_vol_fracP + rhoW_E*water_vol_fracE)
                          +  gas_velocity_P  *(rhoG_P*p_gas_vol_fracP + rhoG_E*p_gas_vol_fracE);


            return ( m_t - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
                + (m_e - m_w)*area()                
                + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc);



            // New friction factor wells
            //OUYANG
            //real_type q_w = ((*m_oil_flow)[ p_node ] + (*m_gas_flow)[ p_node ] + (*m_water_flow)[ p_node ])/dS;
            //real_type v_eq = q_w/(PI*2*m_radius);
            //real_type Re_w = abs(0.5*(rho_P + rho_E)*v_eq*2*m_radius/viscosity);
			//real_type f_P = this->friction_factor( Re )*(1+0.04304*pow(Re_w,0.6142));
            // ASHEIM
            //real_type q_m = p_velocityP*area();
            //real_type q_i = ((*m_oil_flow)[ p_node ] + (*m_gas_flow)[ p_node ] + (*m_water_flow)[ p_node ]);
            //real_type f_complet = q_m == 0? 0.0 : 4*2*m_radius*q_i/dS/q_m + 2*m_radius*q_i/dS/q_m*q_i/dS/q_m;
            //real_type f_P = this->friction_factor( Re ) + f_complet;
            

			//real_type mod_Vgj_P		= this->mod_v_drift_flux( 
			//												 (1-d_w)*p_velocityP + d_w*p_velocityW, p_gas_vol_fracP, p_oil_vol_fracP,
			//												 water_vol_fracP, p_pressureP
			//												 );
			//real_type mod_Vgj_E		= this->mod_v_drift_flux( 
			//												 (1-d_e)*p_velocityP + d_e*p_velocityE, p_gas_vol_fracE, p_oil_vol_fracE,
			//												 water_vol_fracE, p_pressureE
			//												 );

			//// OTHER Vc = j
			////real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
			////real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2*m_radius/viscosity);
			////real_type f_P = this->friction_factor( Re );

			//real_type ksi_e = this->ksi( (1-d_e)*p_velocityP + d_e*p_velocityE );
			//real_type ksi_w = this->ksi( (1-d_w)*p_velocityP + d_w*p_velocityW );
			//
			//
			////*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
			//	 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*p_velocityP*abs(p_velocityP)
			//	 + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
			//	 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);*/
			//
			//real_type ksi_E = this->ksi( p_velocityE );
			//real_type ksi_P = this->ksi( p_velocityP );			

			//real_type m_E	= ((0.5+ksi_E)*rho_E + (0.5-ksi_E)*rho_EE )*p_velocityE*area();
			//real_type m_P	= ((0.5+ksi_P)*rho_P + (0.5-ksi_P)*rho_E  )*p_velocityP*area();
			//real_type m_W	= rho_P*p_velocityP*area();

			//return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + 0.5*(m_E+m_P)*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
			//	 - m_W*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc)
			//     + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE + 1.0e-20)
			//	 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP + 1.0e-20);


		
			
			/*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
				 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
				 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
				 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*m_j[ p_node ]*abs(m_j[ p_node ]);
				*/
				 
				 // + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
				// - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);
			
			}
			


		case 'L':
			{
			real_type dS  = this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
			real_type dSe = 0.;
			real_type dSw = this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node-1 ] );
			real_type dV = this->Volume( dS );

			real_type water_vol_fracP_old	= 1.0 - (m_oil_vol_frac_old[ p_node ]   + m_gas_vol_frac_old[ p_node ]	);
			real_type water_vol_fracE_old	= 1.0 - (m_oil_vol_frac_old[ p_node+1 ] + m_gas_vol_frac_old[ p_node+1 ]);
            real_type water_vol_fracW		= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
            real_type water_vol_fracP		= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
            real_type water_vol_fracE		= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);	
            real_type water_vol_fracEE		= 1.0 - (p_gas_vol_fracEE + p_oil_vol_fracEE);		

			real_type rho_P_old = this->mean_density(m_oil_vol_frac_old[ p_node ], water_vol_fracP_old, m_gas_vol_frac_old[ p_node ], m_pressure_old[ p_node ]);
			real_type rho_E_old = this->mean_density(m_oil_vol_frac_old[ p_node+1 ], water_vol_fracE_old, m_gas_vol_frac_old[ p_node+1 ], m_pressure_old[ p_node+1 ]);

			vector_type S( 3 );
			S[ 0 ] = m_coordinates[ p_node ].getX() - m_coordinates[ p_node+1 ].getX();
			S[ 1 ] = m_coordinates[ p_node ].getY() - m_coordinates[ p_node+1 ].getY();
			S[ 2 ] = m_coordinates[ p_node ].getZ() - m_coordinates[ p_node+1 ].getZ();

			real_type d_e = dS/(dS+dSe);
			real_type d_w = dS/(dS+dSw);

			real_type angle = get_inclination() - PI/2; //- 0*acos( dot( m_gravity, S )/(norm(m_gravity)*norm(S)) );						
			
            real_type rho_W   = this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW , p_pressureW );
            real_type rho_P   = this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);
            real_type rho_E   = this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
            real_type rho_EE  = this->mean_density(p_oil_vol_fracEE, water_vol_fracEE, p_gas_vol_fracEE , p_pressureEE);

            real_type rhoG_W		= this->gas_density( p_pressureW );
            real_type rhoG_P		= this->gas_density( p_pressureP );
            real_type rhoG_E		= this->gas_density( p_pressureE );
            real_type rhoG_EE		= this->gas_density( p_pressureEE );

            real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
            real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
            real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
            real_type rhoL_EE		= this->liquid_density( p_oil_vol_fracEE, water_vol_fracEE, p_pressureEE );

            real_type rhoW_W		= this->water_density( p_pressureW );
            real_type rhoW_P		= this->water_density( p_pressureP );
            real_type rhoW_E		= this->water_density( p_pressureE );
            real_type rhoW_EE		= this->water_density( p_pressureEE );

            real_type rhoO_W		= this->oil_density( p_pressureW );
            real_type rhoO_P		= this->oil_density( p_pressureP );
            real_type rhoO_E		= this->oil_density( p_pressureE );
            real_type rhoO_EE		= this->oil_density( p_pressureEE ); 
			
			real_type mean_pressure = 0.5*(p_pressureP + p_pressureE);
			real_type viscosity = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE)*gas_viscosity	 ( mean_pressure ) 
								+ 0.5*(p_oil_vol_fracP + p_oil_vol_fracE)*oil_viscosity	 ( mean_pressure ) 
								+ 0.5*(water_vol_fracP + water_vol_fracE)*water_viscosity( mean_pressure );


            real_type mod_Vow_W	= this->mod_v_drift_flux_ow( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vow_P	= this->mod_v_drift_flux_ow( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vow_E	= this->mod_v_drift_flux_ow( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );

            real_type mod_Vgj_W	= this->mod_v_drift_flux( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vgj_P	= this->mod_v_drift_flux( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vgj_E	= this->mod_v_drift_flux( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );


            //real_type Vc  = p_velocityP;			
            //real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2.0*m_radius/viscosity);
            //real_type f_P = this->friction_factor( Re );

            // OTHER Vc = j
            real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
            real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2.0*m_radius/viscosity);
            real_type f_P = this->friction_factor( Re );



            real_type gas_vol_frac_W   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
            real_type oil_vol_frac_W   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracW);
            real_type water_vol_frac_W = 0.5*(water_vol_fracP + water_vol_fracW);
            real_type rho_g_W = 0.5*(rhoG_W + rhoG_P);
            real_type rho_o_W = 0.5*(rhoO_W + rhoO_P);
            real_type rho_w_W = 0.5*(rhoW_W + rhoW_P);
            real_type rho_l_W = 0.5*(rhoL_W + rhoL_P);
            real_type rho_m_W = 0.5*(rho_W + rho_P);

            real_type gas_velocity_W    = p_velocityW + rho_l_W/rho_m_W*mod_Vgj_W;
            real_type liquid_velocity_W = p_velocityW - gas_vol_frac_W/(1 - gas_vol_frac_W + 1.0e-20)*rho_g_W/rho_m_W*mod_Vgj_W;	 
            real_type water_velocity_W  = liquid_velocity_W - (oil_vol_frac_W/(water_vol_frac_W + 1.0e-20))*(rho_o_W/rho_l_W)*mod_Vow_W;               
            real_type oil_velocity_W    = liquid_velocity_W + rho_w_W/rho_l_W*mod_Vow_W;  



            real_type gas_vol_frac_P   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE);
            real_type oil_vol_frac_P   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracE);
            real_type water_vol_frac_P = 0.5*(water_vol_fracP + water_vol_fracE);
            real_type rho_g_P = 0.5*(rhoG_E + rhoG_P);
            real_type rho_o_P = 0.5*(rhoO_E + rhoO_P);
            real_type rho_w_P = 0.5*(rhoW_E + rhoW_P);
            real_type rho_l_P = 0.5*(rhoL_E + rhoL_P);
            real_type rho_m_P = 0.5*(rho_E + rho_P);

            real_type gas_velocity_P    = p_velocityP + rho_l_P/rho_m_P*mod_Vgj_P;
            real_type liquid_velocity_P = p_velocityP - gas_vol_frac_P/(1 - gas_vol_frac_P + 1.0e-20)*rho_g_P/rho_m_P*mod_Vgj_P;	 
            real_type water_velocity_P  = liquid_velocity_P - (oil_vol_frac_P/(water_vol_frac_P + 1.0e-20))*(rho_o_P/rho_l_P)*mod_Vow_P;               
            real_type oil_velocity_P    = liquid_velocity_P + rho_w_P/rho_l_P*mod_Vow_P;

            real_type gas_vol_frac_E   = 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE);
            real_type oil_vol_frac_E   = 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE);
            real_type water_vol_frac_E = 0.5*(water_vol_fracE + water_vol_fracEE);
            real_type rho_g_E = 0.5*(rhoG_E + rhoG_EE);
            real_type rho_o_E = 0.5*(rhoO_E + rhoO_EE);
            real_type rho_w_E = 0.5*(rhoW_E + rhoW_EE);
            real_type rho_l_E = 0.5*(rhoL_E + rhoL_EE);
            real_type rho_m_E = 0.5*(rho_E + rho_EE);

            real_type gas_velocity_E    = p_velocityE + rho_l_E/rho_m_E*mod_Vgj_E;
            real_type liquid_velocity_E = p_velocityE - gas_vol_frac_E/(1 - gas_vol_frac_E + 1.0e-20)*rho_g_E/rho_m_E*mod_Vgj_E;	 
            real_type water_velocity_E  = liquid_velocity_E - (oil_vol_frac_E/(water_vol_frac_E + 1.0e-20))*(rho_o_E/rho_l_E)*mod_Vow_E;               
            real_type oil_velocity_E    = liquid_velocity_E + rho_w_E/rho_l_E*mod_Vow_E;



            real_type ksi_W_oil     = this->ksi( oil_velocity_W   );
            real_type ksi_W_gas     = this->ksi( gas_velocity_W   );
            real_type ksi_W_water   = this->ksi( water_velocity_W );

            real_type ksi_P_oil     = this->ksi( oil_velocity_P   );
            real_type ksi_P_gas     = this->ksi( gas_velocity_P   );
            real_type ksi_P_water   = this->ksi( water_velocity_P );

            real_type ksi_E_oil     = this->ksi( oil_velocity_E   );
            real_type ksi_E_gas     = this->ksi( gas_velocity_E   );
            real_type ksi_E_water   = this->ksi( water_velocity_E );

            real_type m_W_oil   = oil_velocity_W  *( (0.5+ksi_W_oil  )*rhoO_W*p_oil_vol_fracW + (0.5-ksi_W_oil  )*rhoO_P*p_oil_vol_fracP );
            real_type m_W_water = water_velocity_W*( (0.5+ksi_W_water)*rhoW_W*water_vol_fracW + (0.5-ksi_W_water)*rhoW_P*water_vol_fracP );
            real_type m_W_gas   = gas_velocity_W  *( (0.5+ksi_W_gas  )*rhoG_W*p_gas_vol_fracW + (0.5-ksi_W_gas  )*rhoG_P*p_gas_vol_fracP );

            real_type m_P_oil   = oil_velocity_P  *( (0.5+ksi_P_oil  )*rhoO_P*p_oil_vol_fracP + (0.5-ksi_P_oil  )*rhoO_E*p_oil_vol_fracE );
            real_type m_P_water = water_velocity_P*( (0.5+ksi_P_water)*rhoW_P*water_vol_fracP + (0.5-ksi_P_water)*rhoW_E*water_vol_fracE );
            real_type m_P_gas   = gas_velocity_P  *( (0.5+ksi_P_gas  )*rhoG_P*p_gas_vol_fracP + (0.5-ksi_P_gas  )*rhoG_E*p_gas_vol_fracE );

            real_type m_E_oil   = oil_velocity_E  *( (0.5+ksi_E_oil  )*rhoO_P*p_oil_vol_fracE + (0.5-ksi_E_oil  )*rhoO_E*p_oil_vol_fracEE );
            real_type m_E_water = water_velocity_E*( (0.5+ksi_E_water)*rhoW_P*water_vol_fracE + (0.5-ksi_E_water)*rhoW_E*water_vol_fracEE );
            real_type m_E_gas   = gas_velocity_E  *( (0.5+ksi_E_gas  )*rhoG_P*p_gas_vol_fracE + (0.5-ksi_E_gas  )*rhoG_E*p_gas_vol_fracEE );

            real_type ksi_e_oil = this->ksi( (1-d_e)*oil_velocity_P + d_e*oil_velocity_E );
            real_type ksi_w_oil = this->ksi( (1-d_w)*oil_velocity_P + d_w*oil_velocity_W );
            real_type ksi_e_water = this->ksi( (1-d_e)*water_velocity_P + d_e*water_velocity_E );
            real_type ksi_w_water = this->ksi( (1-d_w)*water_velocity_P + d_w*water_velocity_W );
            real_type ksi_e_gas = this->ksi( (1-d_e)*gas_velocity_P + d_e*gas_velocity_E );
            real_type ksi_w_gas = this->ksi( (1-d_w)*gas_velocity_P + d_w*gas_velocity_W );

            real_type m_e = 0.5*(m_E_oil  +m_P_oil)  *( (0.5+ksi_e_oil  )*oil_velocity_P   + (0.5-ksi_e_oil  )*oil_velocity_E )
                + 0.5*(m_E_water+m_P_water)*( (0.5+ksi_e_water)*water_velocity_P + (0.5-ksi_e_water)*water_velocity_E )
                + 0.5*(m_E_gas  +m_P_gas)  *( (0.5+ksi_e_gas  )*gas_velocity_P   + (0.5-ksi_e_gas  )*gas_velocity_E );

            real_type m_w = 0.5*(m_W_oil  +m_P_oil)  *( (0.5+ksi_w_oil  )*oil_velocity_W   + (0.5-ksi_w_oil  )*oil_velocity_P )
                + 0.5*(m_W_water+m_P_water)*( (0.5+ksi_w_water)*water_velocity_W + (0.5-ksi_w_water)*water_velocity_P )
                + 0.5*(m_W_gas  +m_P_gas)  *( (0.5+ksi_w_gas  )*gas_velocity_W   + (0.5-ksi_w_gas  )*gas_velocity_P );  


            real_type m_t =  oil_velocity_P  *(rhoO_P*p_oil_vol_fracP + rhoO_E*p_oil_vol_fracE)
                          +  water_velocity_P*(rhoW_P*water_vol_fracP + rhoW_E*water_vol_fracE)
                          +  gas_velocity_P  *(rhoG_P*p_gas_vol_fracP + rhoG_E*p_gas_vol_fracE);

            return ( m_t - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
                + (m_e - m_w)*area()                
                + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc);





            // New friction factor wells
            // OUYANG
            //real_type q_w = 0.5*((*m_oil_flow)[ p_node ]+(*m_oil_flow)[ p_node +1] + (*m_gas_flow)[ p_node ]+(*m_gas_flow)[ p_node +1] + (*m_water_flow)[ p_node ]+(*m_water_flow)[ p_node+1 ])/dS;
            //real_type v_eq = q_w/(PI*2*m_radius);
            //real_type Re_w = abs(0.5*(rho_P + rho_E)*v_eq*2*m_radius/viscosity);
            //real_type f_P = this->friction_factor( Re )*(1+0.04304*pow(Re_w,0.6142));
            // ASHEIM
            //real_type q_m = p_velocityP*area();
            //real_type q_i = 0.5*((*m_oil_flow)[ p_node ] + (*m_oil_flow)[ p_node+1 ] + (*m_gas_flow)[ p_node ]+ (*m_gas_flow)[ p_node+1 ] + (*m_water_flow)[ p_node ] + (*m_water_flow)[ p_node+1 ]);
            //real_type f_complet = q_m == 0? 0.0 : 4*2*m_radius*q_i/dS/q_m + 2*m_radius*q_i/dS/q_m*q_i/dS/q_m;
            //real_type f_P = this->friction_factor( Re ) + f_complet;
            
			//real_type mod_Vgj_P		= this->mod_v_drift_flux( 
			//												 (1-d_w)*p_velocityP + d_w*p_velocityW, p_gas_vol_fracP, p_oil_vol_fracP,
			//												 water_vol_fracP, p_pressureP
			//												 );
			//real_type mod_Vgj_E		= this->mod_v_drift_flux( 
			//												 (1-d_e)*p_velocityP + d_e*p_velocityE, p_gas_vol_fracE, p_oil_vol_fracE,
			//												 water_vol_fracE, p_pressureE
			//												 );

			//
			//// OTHER Vc = j
			////real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
			////real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2*m_radius/viscosity);
			////real_type f_P = this->friction_factor( Re );

			//real_type ksi_e = this->ksi( (1-d_e)*p_velocityP + d_e*p_velocityE );
			//real_type ksi_w = this->ksi( (1-d_w)*p_velocityP + d_w*p_velocityW );
			//
			//
			////*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE ) 
			//	 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*p_velocityP*abs(p_velocityP)
			//	 + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
			//	 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);*/
			//
			//
			//real_type ksi_P = this->ksi( p_velocityP );
			//real_type ksi_W = this->ksi( p_velocityW );

			//real_type m_E	= rho_E*p_velocityE*area();
			//real_type m_P	= ((0.5+ksi_P)*rho_P + (0.5-ksi_P)*rho_E  )*p_velocityP*area();
			//real_type m_W	= ((0.5+ksi_W)*rho_W + (0.5-ksi_W)*rho_P  )*p_velocityW*area();

   //                    
			//return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + 0.5*(m_E+m_P)*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
			//	 - 0.5*(m_P+m_W)*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc)
			//	 + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE + 1.0e-20)
			//	 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP + 1.0e-20);
			//
			////*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE ) 
			//	 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*m_j[ p_node ]*abs(m_j[ p_node ]);
			//	 */
			//	 
			//	 //+ rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
			//	// - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);

            
			}
			


		case 'C':
			{
			real_type dS  = this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node+1 ] );
			real_type dSe = this->segment_length( m_coordinates[ p_node+1 ], m_coordinates[ p_node+2 ] );
			real_type dSw = this->segment_length( m_coordinates[ p_node ], m_coordinates[ p_node-1 ] );
			real_type dV = this->Volume( dS );

			real_type water_vol_fracP_old	= 1.0 - (m_oil_vol_frac_old[ p_node ]   + m_gas_vol_frac_old[ p_node ]	);
			real_type water_vol_fracE_old	= 1.0 - (m_oil_vol_frac_old[ p_node+1 ] + m_gas_vol_frac_old[ p_node+1 ]);
			real_type water_vol_fracW		= 1.0 - (p_gas_vol_fracW + p_oil_vol_fracW);
			real_type water_vol_fracP		= 1.0 - (p_gas_vol_fracP + p_oil_vol_fracP);
			real_type water_vol_fracE		= 1.0 - (p_gas_vol_fracE + p_oil_vol_fracE);	
			real_type water_vol_fracEE		= 1.0 - (p_gas_vol_fracEE + p_oil_vol_fracEE);

			real_type rho_P_old = this->mean_density(m_oil_vol_frac_old[ p_node ], water_vol_fracP_old, m_gas_vol_frac_old[ p_node ], m_pressure_old[ p_node ]);
			real_type rho_E_old = this->mean_density(m_oil_vol_frac_old[ p_node+1 ], water_vol_fracE_old, m_gas_vol_frac_old[ p_node+1 ], m_pressure_old[ p_node+1 ]);
            real_type rhoG_P_old = this->gas_density( m_pressure_old[ p_node ] );
            real_type rhoO_P_old = this->oil_density( m_pressure_old[ p_node ] );
            real_type rhoW_P_old = this->water_density( m_pressure_old[ p_node ] );

			vector_type S( 3 );
			S[ 0 ] = m_coordinates[ p_node ].getX() - m_coordinates[ p_node+1 ].getX();
			S[ 1 ] = m_coordinates[ p_node ].getY() - m_coordinates[ p_node+1 ].getY();
			S[ 2 ] = m_coordinates[ p_node ].getZ() - m_coordinates[ p_node+1 ].getZ();
			
			real_type d_e = dS/(dS+dSe);
			real_type d_w = dS/(dS+dSw);

			real_type angle = get_inclination() - PI/2;// - 0*acos( dot( m_gravity, S )/(norm(m_gravity)*norm(S)) );		
			
			real_type rho_W   = this->mean_density(p_oil_vol_fracW, water_vol_fracW, p_gas_vol_fracW , p_pressureW );
			real_type rho_P   = this->mean_density(p_oil_vol_fracP, water_vol_fracP, p_gas_vol_fracP, p_pressureP);
			real_type rho_E   = this->mean_density(p_oil_vol_fracE, water_vol_fracE, p_gas_vol_fracE, p_pressureE);
			real_type rho_EE  = this->mean_density(p_oil_vol_fracEE, water_vol_fracEE, p_gas_vol_fracEE , p_pressureEE);
			
            real_type rhoG_W		= this->gas_density( p_pressureW );
			real_type rhoG_P		= this->gas_density( p_pressureP );
			real_type rhoG_E		= this->gas_density( p_pressureE );
            real_type rhoG_EE		= this->gas_density( p_pressureEE );

            real_type rhoL_W		= this->liquid_density( p_oil_vol_fracW, water_vol_fracW, p_pressureW );
			real_type rhoL_P		= this->liquid_density( p_oil_vol_fracP, water_vol_fracP, p_pressureP );
			real_type rhoL_E		= this->liquid_density( p_oil_vol_fracE, water_vol_fracE, p_pressureE );
            real_type rhoL_EE		= this->liquid_density( p_oil_vol_fracEE, water_vol_fracEE, p_pressureEE );

            real_type rhoW_W		= this->water_density( p_pressureW );
            real_type rhoW_P		= this->water_density( p_pressureP );
            real_type rhoW_E		= this->water_density( p_pressureE );
            real_type rhoW_EE		= this->water_density( p_pressureEE );

            real_type rhoO_W		= this->oil_density( p_pressureW );
            real_type rhoO_P		= this->oil_density( p_pressureP );
            real_type rhoO_E		= this->oil_density( p_pressureE );
            real_type rhoO_EE		= this->oil_density( p_pressureEE );            
			
			real_type mean_pressure = 0.5*(p_pressureP + p_pressureE);
			real_type viscosity = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE)*gas_viscosity	 ( mean_pressure ) 
								+ 0.5*(p_oil_vol_fracP + p_oil_vol_fracE)*oil_viscosity	 ( mean_pressure ) 
								+ 0.5*(water_vol_fracP + water_vol_fracE)*water_viscosity( mean_pressure );	


            // New friction factor wells
            // OUYANG
            //real_type q_w = 0.5*((*m_oil_flow)[ p_node ]+(*m_oil_flow)[ p_node +1] + (*m_gas_flow)[ p_node ]+(*m_gas_flow)[ p_node +1] + (*m_water_flow)[ p_node ]+(*m_water_flow)[ p_node+1 ])/dS;
            //real_type v_eq = q_w/(PI*2*m_radius);
            //real_type Re_w = abs(0.5*(rho_P + rho_E)*v_eq*2*m_radius/viscosity);
            //real_type f_P = this->friction_factor( Re )*(1+0.04304*pow(Re_w,0.6142));
            // ASHEIM
            //real_type q_m = p_velocityP*area();
            //real_type q_i = 0.5*((*m_oil_flow)[ p_node ] + (*m_oil_flow)[ p_node+1 ] + (*m_gas_flow)[ p_node ]+ (*m_gas_flow)[ p_node+1 ] + (*m_water_flow)[ p_node ] + (*m_water_flow)[ p_node+1 ]);
            //real_type f_complet = q_m == 0? 0.0 : 4*2*m_radius*q_i/dS/q_m + 2*m_radius*q_i/dS/q_m*q_i/dS/q_m;
            //real_type f_P = this->friction_factor( Re ) + f_complet;
           

			/*real_type mod_Vgj_P		= this->mod_v_drift_flux( 
															(1-d_w)*p_velocityP + d_w*p_velocityW, p_gas_vol_fracP, p_oil_vol_fracP,
															 water_vol_fracP, p_pressureP
															 );
			real_type mod_Vgj_E		= this->mod_v_drift_flux( 
															 (1-d_e)*p_velocityP + d_e*p_velocityE, p_gas_vol_fracE, p_oil_vol_fracE,
															 water_vol_fracE, p_pressureE
															 );
            real_type mod_Vow_P	= this->mod_v_drift_flux_ow( 
															(1-d_w)*p_velocityP + d_w*p_velocityW, p_gas_vol_fracP, p_oil_vol_fracP,
															 water_vol_fracP, p_pressureP
															 );
			real_type mod_Vow_E	= this->mod_v_drift_flux_ow( 
															 (1-d_e)*p_velocityP + d_e*p_velocityE, p_gas_vol_fracE, p_oil_vol_fracE,
															 water_vol_fracE, p_pressureE
															 );*/


            real_type mod_Vow_W	= this->mod_v_drift_flux_ow( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vow_P	= this->mod_v_drift_flux_ow( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vow_E	= this->mod_v_drift_flux_ow( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );

            real_type mod_Vgj_W	= this->mod_v_drift_flux( 
                p_velocityW, 0.5*(p_gas_vol_fracP + p_gas_vol_fracW), 0.5*(p_oil_vol_fracP + p_oil_vol_fracW),
                0.5*(water_vol_fracP + water_vol_fracW), 0.5*(p_pressureP + p_pressureW)
                );
            real_type mod_Vgj_P	= this->mod_v_drift_flux( 
                p_velocityP, 0.5*(p_gas_vol_fracP + p_gas_vol_fracE), 0.5*(p_oil_vol_fracP + p_oil_vol_fracE),
                0.5*(water_vol_fracP + water_vol_fracE), 0.5*(p_pressureP + p_pressureE)
                );
            real_type mod_Vgj_E	= this->mod_v_drift_flux( 
                p_velocityE, 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE), 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE),
                0.5*(water_vol_fracE + water_vol_fracEE), 0.5*(p_pressureE + p_pressureEE)
                );
                

            //real_type Vc  = p_velocityP;			
            //real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2.0*m_radius/viscosity);
            //real_type f_P = this->friction_factor( Re );

            // OTHER Vc = j
            real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
            real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2.0*m_radius/viscosity);
            real_type f_P = this->friction_factor( Re );



            real_type gas_vol_frac_W   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracW);
            real_type oil_vol_frac_W   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracW);
            real_type water_vol_frac_W = 0.5*(water_vol_fracP + water_vol_fracW);
            real_type rho_g_W = 0.5*(rhoG_W + rhoG_P);
            real_type rho_o_W = 0.5*(rhoO_W + rhoO_P);
            real_type rho_w_W = 0.5*(rhoW_W + rhoW_P);
            real_type rho_l_W = 0.5*(rhoL_W + rhoL_P);
            real_type rho_m_W = 0.5*(rho_W + rho_P);

            real_type gas_velocity_W    = p_velocityW + rho_l_W/rho_m_W*mod_Vgj_W;
            real_type liquid_velocity_W = p_velocityW - gas_vol_frac_W/(1 - gas_vol_frac_W + 1.0e-20)*rho_g_W/rho_m_W*mod_Vgj_W;	 
            real_type water_velocity_W  = liquid_velocity_W - (oil_vol_frac_W/(water_vol_frac_W + 1.0e-20))*(rho_o_W/rho_l_W)*mod_Vow_W;               
            real_type oil_velocity_W    = liquid_velocity_W + rho_w_W/rho_l_W*mod_Vow_W;  



            real_type gas_vol_frac_P   = 0.5*(p_gas_vol_fracP + p_gas_vol_fracE);
            real_type oil_vol_frac_P   = 0.5*(p_oil_vol_fracP + p_oil_vol_fracE);
            real_type water_vol_frac_P = 0.5*(water_vol_fracP + water_vol_fracE);
            real_type rho_g_P = 0.5*(rhoG_E + rhoG_P);
            real_type rho_o_P = 0.5*(rhoO_E + rhoO_P);
            real_type rho_w_P = 0.5*(rhoW_E + rhoW_P);
            real_type rho_l_P = 0.5*(rhoL_E + rhoL_P);
            real_type rho_m_P = 0.5*(rho_E + rho_P);

            real_type gas_velocity_P    = p_velocityP + rho_l_P/rho_m_P*mod_Vgj_P;
            real_type liquid_velocity_P = p_velocityP - gas_vol_frac_P/(1 - gas_vol_frac_P + 1.0e-20)*rho_g_P/rho_m_P*mod_Vgj_P;	 
            real_type water_velocity_P  = liquid_velocity_P - (oil_vol_frac_P/(water_vol_frac_P + 1.0e-20))*(rho_o_P/rho_l_P)*mod_Vow_P;               
            real_type oil_velocity_P    = liquid_velocity_P + rho_w_P/rho_l_P*mod_Vow_P;

            real_type gas_vol_frac_E   = 0.5*(p_gas_vol_fracE + p_gas_vol_fracEE);
            real_type oil_vol_frac_E   = 0.5*(p_oil_vol_fracE + p_oil_vol_fracEE);
            real_type water_vol_frac_E = 0.5*(water_vol_fracE + water_vol_fracEE);
            real_type rho_g_E = 0.5*(rhoG_E + rhoG_EE);
            real_type rho_o_E = 0.5*(rhoO_E + rhoO_EE);
            real_type rho_w_E = 0.5*(rhoW_E + rhoW_EE);
            real_type rho_l_E = 0.5*(rhoL_E + rhoL_EE);
            real_type rho_m_E = 0.5*(rho_E + rho_EE);

            real_type gas_velocity_E    = p_velocityE + rho_l_E/rho_m_E*mod_Vgj_E;
            real_type liquid_velocity_E = p_velocityE - gas_vol_frac_E/(1 - gas_vol_frac_E + 1.0e-20)*rho_g_E/rho_m_E*mod_Vgj_E;	 
            real_type water_velocity_E  = liquid_velocity_E - (oil_vol_frac_E/(water_vol_frac_E + 1.0e-20))*(rho_o_E/rho_l_E)*mod_Vow_E;               
            real_type oil_velocity_E    = liquid_velocity_E + rho_w_E/rho_l_E*mod_Vow_E;



            real_type ksi_W_oil     = this->ksi( oil_velocity_W   );
            real_type ksi_W_gas     = this->ksi( gas_velocity_W   );
            real_type ksi_W_water   = this->ksi( water_velocity_W );

            real_type ksi_P_oil     = this->ksi( oil_velocity_P   );
            real_type ksi_P_gas     = this->ksi( gas_velocity_P   );
            real_type ksi_P_water   = this->ksi( water_velocity_P );

            real_type ksi_E_oil     = this->ksi( oil_velocity_E   );
            real_type ksi_E_gas     = this->ksi( gas_velocity_E   );
            real_type ksi_E_water   = this->ksi( water_velocity_E );

            real_type m_W_oil   = oil_velocity_W  *( (0.5+ksi_W_oil  )*rhoO_W*p_oil_vol_fracW + (0.5-ksi_W_oil  )*rhoO_P*p_oil_vol_fracP );
            real_type m_W_water = water_velocity_W*( (0.5+ksi_W_water)*rhoW_W*water_vol_fracW + (0.5-ksi_W_water)*rhoW_P*water_vol_fracP );
            real_type m_W_gas   = gas_velocity_W  *( (0.5+ksi_W_gas  )*rhoG_W*p_gas_vol_fracW + (0.5-ksi_W_gas  )*rhoG_P*p_gas_vol_fracP );

            real_type m_P_oil   = oil_velocity_P  *( (0.5+ksi_P_oil  )*rhoO_P*p_oil_vol_fracP + (0.5-ksi_P_oil  )*rhoO_E*p_oil_vol_fracE );
            real_type m_P_water = water_velocity_P*( (0.5+ksi_P_water)*rhoW_P*water_vol_fracP + (0.5-ksi_P_water)*rhoW_E*water_vol_fracE );
            real_type m_P_gas   = gas_velocity_P  *( (0.5+ksi_P_gas  )*rhoG_P*p_gas_vol_fracP + (0.5-ksi_P_gas  )*rhoG_E*p_gas_vol_fracE );

            real_type m_E_oil   = oil_velocity_E  *( (0.5+ksi_E_oil  )*rhoO_P*p_oil_vol_fracE + (0.5-ksi_E_oil  )*rhoO_E*p_oil_vol_fracEE );
            real_type m_E_water = water_velocity_E*( (0.5+ksi_E_water)*rhoW_P*water_vol_fracE + (0.5-ksi_E_water)*rhoW_E*water_vol_fracEE );
            real_type m_E_gas   = gas_velocity_E  *( (0.5+ksi_E_gas  )*rhoG_P*p_gas_vol_fracE + (0.5-ksi_E_gas  )*rhoG_E*p_gas_vol_fracEE );

            real_type ksi_e_oil = this->ksi( (1-d_e)*oil_velocity_P + d_e*oil_velocity_E );
            real_type ksi_w_oil = this->ksi( (1-d_w)*oil_velocity_P + d_w*oil_velocity_W );
            real_type ksi_e_water = this->ksi( (1-d_e)*water_velocity_P + d_e*water_velocity_E );
            real_type ksi_w_water = this->ksi( (1-d_w)*water_velocity_P + d_w*water_velocity_W );
            real_type ksi_e_gas = this->ksi( (1-d_e)*gas_velocity_P + d_e*gas_velocity_E );
            real_type ksi_w_gas = this->ksi( (1-d_w)*gas_velocity_P + d_w*gas_velocity_W );

            real_type m_e = 0.5*(m_E_oil  +m_P_oil)  *( (0.5+ksi_e_oil  )*oil_velocity_P   + (0.5-ksi_e_oil  )*oil_velocity_E )
                          + 0.5*(m_E_water+m_P_water)*( (0.5+ksi_e_water)*water_velocity_P + (0.5-ksi_e_water)*water_velocity_E )
                          + 0.5*(m_E_gas  +m_P_gas)  *( (0.5+ksi_e_gas  )*gas_velocity_P   + (0.5-ksi_e_gas  )*gas_velocity_E );

            real_type m_w = 0.5*(m_W_oil  +m_P_oil)  *( (0.5+ksi_w_oil  )*oil_velocity_W   + (0.5-ksi_w_oil  )*oil_velocity_P )
                          + 0.5*(m_W_water+m_P_water)*( (0.5+ksi_w_water)*water_velocity_W + (0.5-ksi_w_water)*water_velocity_P )
                          + 0.5*(m_W_gas  +m_P_gas)  *( (0.5+ksi_w_gas  )*gas_velocity_W   + (0.5-ksi_w_gas  )*gas_velocity_P );

          /*  real_type m_e = rhoO_E*p_oil_vol_fracE*((1-d_e)*oil_velocity_P + d_e*oil_velocity_E )*( (0.5+ksi_e_oil  )*oil_velocity_P   + (0.5-ksi_e_oil  )*oil_velocity_E )
                          + rhoW_E*water_vol_fracE*((1-d_e)*water_velocity_P + d_e*water_velocity_E )*( (0.5+ksi_e_water)*water_velocity_P + (0.5-ksi_e_water)*water_velocity_E )
                          + rhoG_E*p_gas_vol_fracE*((1-d_e)*gas_velocity_P + d_e*gas_velocity_E )*( (0.5+ksi_e_gas  )*gas_velocity_P   + (0.5-ksi_e_gas  )*gas_velocity_E );

            real_type m_w = rhoO_P*p_oil_vol_fracP*((1-d_w)*oil_velocity_P + d_w*oil_velocity_W )*( (0.5+ksi_w_oil  )*oil_velocity_W   + (0.5-ksi_w_oil  )*oil_velocity_P )
                          + rhoW_P*water_vol_fracP*((1-d_w)*water_velocity_P + d_w*water_velocity_W )*( (0.5+ksi_w_water)*water_velocity_W + (0.5-ksi_w_water)*water_velocity_P )
                          + rhoG_P*p_gas_vol_fracP*((1-d_w)*gas_velocity_P + d_w*gas_velocity_W )*( (0.5+ksi_w_gas  )*gas_velocity_W   + (0.5-ksi_w_gas  )*gas_velocity_P );
        */    
            real_type m_t =  oil_velocity_P  *(rhoO_P*p_oil_vol_fracP + rhoO_E*p_oil_vol_fracE)
                          +  water_velocity_P*(rhoW_P*water_vol_fracP + rhoW_E*water_vol_fracE)
                          +  gas_velocity_P  *(rhoG_P*p_gas_vol_fracP + rhoG_E*p_gas_vol_fracE);

            return ( m_t - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
                + (m_e - m_w)*area()                
                + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc);
           


			// OTHER Vc = j
			//real_type Vc = p_velocityP + 0.5*(p_gas_vol_fracP*(rhoL_P - rhoG_P)/rho_P*mod_Vgj_P + p_gas_vol_fracE*(rhoL_E - rhoG_E)/rho_E*mod_Vgj_E);
			//real_type Re  = abs(0.5*(rho_P + rho_E)*Vc*2*m_radius/viscosity);
			//real_type f_P = this->friction_factor( Re );

		
		/*	real_type ksi_e = this->ksi( (1-d_e)*p_velocityP + d_e*p_velocityE );
			real_type ksi_w = this->ksi( (1-d_w)*p_velocityP + d_w*p_velocityW );*/
			
			
			/*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
				 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
				 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
				 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*p_velocityP*abs(p_velocityP)
				 + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
				 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);*/
			//real_type ksi_E = this->ksi( p_velocityE );
			//real_type ksi_P = this->ksi( p_velocityP );
			//real_type ksi_W = this->ksi( p_velocityW );

			//real_type m_E	= ((0.5+ksi_E)*rho_E + (0.5-ksi_E)*rho_EE )*p_velocityE*area();
			//real_type m_P	= ((0.5+ksi_P)*rho_P + (0.5-ksi_P)*rho_E  )*p_velocityP*area();
			//real_type m_W	= ((0.5+ksi_W)*rho_W + (0.5-ksi_W)*rho_P  )*p_velocityW*area();

			//return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
			//	 + 0.5*(m_E+m_P)*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
			//	 - 0.5*(m_P+m_W)*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
			//	 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*Vc*abs(Vc)
			//	 + rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE + 1.0e-20)
			//	 - rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP + 1.0e-20);

			/*return ( (rho_P+rho_E)*p_velocityP - (rho_P_old+rho_E_old)*m_mean_velocity_old[ p_node ] )*0.5*dV/dt()
				 + rho_E*area()*( (1-d_e)*p_velocityP + d_e*p_velocityE )*( (0.5+ksi_e)*p_velocityP + (0.5-ksi_e)*p_velocityE )
				 - rho_P*area()*( (1-d_w)*p_velocityP + d_w*p_velocityW )*( (0.5+ksi_w)*p_velocityW + (0.5-ksi_w)*p_velocityP )
				 + (p_pressureE-p_pressureP)*area() + 0.5*(rho_P+rho_E)*gravity()*sin( angle )*dV + 0.125/m_radius*f_P*(rho_P + rho_E)*dV*m_j[ p_node ]*abs(m_j[ p_node ]);
			*/
				 
				 //+ rhoG_E*rhoL_E/rho_E*area()*mod_Vgj_E*mod_Vgj_E*p_gas_vol_fracE/(1-p_gas_vol_fracE)
				 //- rhoG_P*rhoL_P/rho_P*area()*mod_Vgj_P*mod_Vgj_P*p_gas_vol_fracP/(1-p_gas_vol_fracP);
			}

		default:
			return 0.;
			
		}
	}

	

	void DriftFluxWell::GMRES_Solve( smatrix_type &A, svector_type &x, svector_type &b )
	{
		m_convergence_status = true;
            			
        itl::ILU<smatrix_type> precond(A);
        // SSOR preconditioner			
        //itl::SSOR<smatrix_type> precond(A);		
        svector_type b2( A.ncols() );			
        itl::solve(precond(), b, b2); //gmres needs the preconditioned b to pass into iter object.
        //iteration
        int max_iter = 1000;	//ex: 1000			
        itl::noisy_iteration<double> iter(b2, max_iter, 0.0, 1E-6);
        int restart = 10; //restart constant: 10
        // modified_gram_schmidt				
        itl::modified_gram_schmidt<svector_type> orth( restart, x.size() );			
        //gmres algorithm	            
        m_convergence_status = itl::gmres(A, x, b, precond(), restart, iter, orth); 
			
	
	}

	void DriftFluxWell::compute_Jacobian()
	{
		bool WITH_GAS = this->m_with_gas;
		bool WITH_MOMENTUM = true;

		real_type s_R_m;			real_type s_R_g;			real_type s_R_o;			real_type s_R_v;
		real_type R_m_dPW;			real_type R_g_dPW;			real_type R_o_dPW;			real_type R_v_dPP;  
		real_type R_m_dPP;			real_type R_g_dPP;			real_type R_o_dPP;			real_type R_v_dPE; 
		real_type R_m_dPE;			real_type R_g_dPE;			real_type R_o_dPE;			real_type R_v_dalphaGasP;
		real_type R_m_dalphaGasW;	real_type R_g_dalphaGasW;	real_type R_o_dalphaGasW;	real_type R_v_dalphaGasE;  
		real_type R_m_dalphaGasP;	real_type R_g_dalphaGasP;	real_type R_o_dalphaGasP;	real_type R_v_dalphaOilP; 
		real_type R_m_dalphaGasE;	real_type R_g_dalphaGasE;	real_type R_o_dalphaGasE;	real_type R_v_dalphaOilE;
		real_type R_m_dalphaOilW;	real_type R_g_dalphaOilW;	real_type R_o_dalphaOilW;	real_type R_v_dvW;  
		real_type R_m_dalphaOilP;	real_type R_g_dalphaOilP;	real_type R_o_dalphaOilP;	real_type R_v_dvP; 
		real_type R_m_dalphaOilE;	real_type R_g_dalphaOilE;	real_type R_o_dalphaOilE;	real_type R_v_dvE;
		real_type R_m_dvW;			real_type R_g_dvW;			real_type R_o_dvW;			
		real_type R_m_dvP;			real_type R_g_dvP;			real_type R_o_dvP;			
		
        // EXTRA DERIVATIVES 
		real_type R_v_dPW;
        real_type R_v_dPEE;
        real_type R_v_dalphaGasW;
        real_type R_v_dalphaGasEE;
        real_type R_v_dalphaOilW;
        real_type R_v_dalphaOilEE;

		real_type delta_PP;	
		real_type delta_alphaGasP;
		real_type delta_alphaOilP;
		real_type delta_vP;
		real_type delta_PW;	
		real_type delta_alphaGasW;
		real_type delta_alphaOilW;
		real_type delta_vW;	
		real_type delta_PE;	
		real_type delta_alphaGasE;
		real_type delta_alphaOilE;
		real_type delta_vE;	
        real_type delta_PEE;	
        real_type delta_alphaGasEE;
        real_type delta_alphaOilEE;
		
		uint_type CENT = 0;
		uint_type EAST = CENT + 1;
		uint_type WEST = CENT;		
				
		

		delta_PP		= m_delta[ P ]*m_pressure[ CENT ];
		delta_alphaGasP	= m_gas_vol_frac[ CENT ] > 1e-8 ? m_delta[ alpha_g ]*m_gas_vol_frac[ CENT ] : 1e-4*m_delta[ alpha_g ];
		delta_alphaOilP	= m_oil_vol_frac[ CENT ] > 1e-8 ? m_delta[ alpha_o ]*m_oil_vol_frac[ CENT ] : 1e-4*m_delta[ alpha_o ];
		delta_vP		= abs(m_mean_velocity[ CENT ]) > 1e-8 ? m_delta[ v ]*m_mean_velocity[ CENT ] : 1e-4*m_delta[ v ];		
		delta_PE		= m_delta[ P ]*m_pressure[ EAST ];
		delta_alphaGasE	= m_gas_vol_frac[ EAST ] > 1e-8 ? m_delta[ alpha_g ]*m_gas_vol_frac[ EAST ] : 1e-4*m_delta[ alpha_g ];
		delta_alphaOilE	= m_oil_vol_frac[ EAST ] > 1e-8 ? m_delta[ alpha_o ]*m_oil_vol_frac[ EAST ] : 1e-4*m_delta[ alpha_o ];
		delta_vE		= abs(m_mean_velocity[ EAST ]) > 1e-8 ? m_delta[ v ]*m_mean_velocity[ EAST ] : 1e-4*m_delta[ v ];
        delta_PEE		    = m_delta[ P ]*m_pressure[ EAST+1 ];
        delta_alphaGasEE	= m_gas_vol_frac[ EAST+1 ] > 1e-8 ? m_delta[ alpha_g ]*m_gas_vol_frac[ EAST+1 ] : 1e-4*m_delta[ alpha_g ];
        delta_alphaOilEE	= m_oil_vol_frac[ EAST+1 ] > 1e-8 ? m_delta[ alpha_o ]*m_oil_vol_frac[ EAST+1 ] : 1e-4*m_delta[ alpha_o ];

	//MIXTURE - PRESSURE == 0		
		(*this->m_matrix)( id(0,P), id(0,P) )		= 1.;		
		
	//GAS
		if( WITH_GAS ){			
			(*this->m_matrix)( id(0,alpha_g), id(0,alpha_g) ) = 1.;			
		}
		else{
			(*this->m_matrix)( id(0,alpha_g), id(0,alpha_g) ) = 1.;
		}
	//OIL
		(*this->m_matrix)( id(0,alpha_o), id(0,alpha_o) ) = 1.;

		if( WITH_MOMENTUM ){
			//MOMENTUM						
			s_R_v		   = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dPP		   = this->R_v(m_pressure[ CENT ]  + delta_PP, m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ],
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dPE		   = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
            R_v_dPEE	   = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ] + delta_PEE, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaGasP = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaGasE = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE, m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaGasEE= this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ] + delta_alphaGasEE, 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaOilP = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaOilE = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dalphaOilEE= this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ] + delta_alphaOilEE, m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], 0, 'F');
			R_v_dvP		   = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ] + delta_vP, m_mean_velocity[ EAST ], 0, 'F');
			R_v_dvE		   = this->R_v(m_pressure[ CENT ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
									   m_oil_vol_frac[ CENT ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ CENT ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ] + delta_vE, 0, 'F');
			//CENTRAL		
			(*this->m_matrix)( id(0,v), id(0,P) )	    = (R_v_dPP - s_R_v)/delta_PP;		
			(*this->m_matrix)( id(0,v), id(0,alpha_g) ) = (R_v_dalphaGasP - s_R_v)/delta_alphaGasP;
			(*this->m_matrix)( id(0,v), id(0,alpha_o) ) = (R_v_dalphaOilP - s_R_v)/delta_alphaOilP;
			(*this->m_matrix)( id(0,v), id(0,v) )	    = (R_v_dvP - s_R_v)/delta_vP;
			//EAST	
			(*this->m_matrix)( id(0,v), id(0,P) + total_var )	    = (R_v_dPE - s_R_v)/delta_PE;
			(*this->m_matrix)( id(0,v), id(0,alpha_g) + total_var ) = (R_v_dalphaGasE - s_R_v)/delta_alphaGasE;
			(*this->m_matrix)( id(0,v), id(0,alpha_o) + total_var ) = (R_v_dalphaOilE - s_R_v)/delta_alphaOilE;
			(*this->m_matrix)( id(0,v), id(0,v) + total_var )	    = (R_v_dvE - s_R_v)/delta_vE;
			(*this->m_source)[ id(0, v) ]						    = -s_R_v;
            //EEAST	
            (*this->m_matrix)( id(0,v), id(0,P) + 2*total_var )	      = (R_v_dPEE - s_R_v)/delta_PEE;
            (*this->m_matrix)( id(0,v), id(0,alpha_g) + 2*total_var ) = (R_v_dalphaGasEE - s_R_v)/delta_alphaGasEE;
            (*this->m_matrix)( id(0,v), id(0,alpha_o) + 2*total_var ) = (R_v_dalphaOilEE - s_R_v)/delta_alphaOilEE;
		}
		else{
			(*this->m_matrix)( id(0,v), id(0,P) + total_var ) = 1.;			
		}

		

		++CENT;
		++EAST;
		uint_type LAST = this->number_of_nodes()-1;
		for( uint_type i = 1; i < LAST - 1; ++i )
		{

			delta_PP		= m_delta[ P ]*m_pressure[ CENT ];
			delta_alphaGasP	= m_gas_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ CENT ] : m_delta[ alpha_g ];
			delta_alphaOilP	= m_oil_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ CENT ] : m_delta[ alpha_o ];
			delta_vP		= abs(m_mean_velocity[ CENT ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ CENT ] : m_delta[ v ];
						
			delta_PE		= m_delta[ P ]*m_pressure[ EAST ];
			delta_alphaGasE	= m_gas_vol_frac[ EAST ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ EAST ] : m_delta[ alpha_g ];
			delta_alphaOilE	= m_oil_vol_frac[ EAST ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ EAST ] : m_delta[ alpha_o ];
			delta_vE		= abs(m_mean_velocity[ EAST ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ EAST ] : m_delta[ v ];

            delta_PEE		    = m_delta[ P ]*m_pressure[ EAST+1 ];
            delta_alphaGasEE	= m_gas_vol_frac[ EAST+1 ] > 1e-8 ? m_delta[ alpha_g ]*m_gas_vol_frac[ EAST+1 ] : 1e-4*m_delta[ alpha_g ];
            delta_alphaOilEE	= m_oil_vol_frac[ EAST+1 ] > 1e-8 ? m_delta[ alpha_o ]*m_oil_vol_frac[ EAST+1 ] : 1e-4*m_delta[ alpha_o ];

			delta_PW		= m_delta[ P ]*m_pressure[ WEST ];
			delta_alphaGasW	= m_gas_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ WEST ] : m_delta[ alpha_g ];
			delta_alphaOilW	= m_oil_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ WEST ] : m_delta[ alpha_o ];
			delta_vW		= abs(m_mean_velocity[ WEST ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ WEST ] : m_delta[ v ];
			
		//MIXTURE
			s_R_m		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dPW		   = this->R_m(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dPP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dPE		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaGasW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaGasP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaGasE = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaOilW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaOilP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dalphaOilE = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
			R_m_dvW		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], i, 'C');
			R_m_dvP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, i, 'C');
			//WEST			
			(*this->m_matrix)( id(i,P), id(i,P) - total_var )		= (R_m_dPW - s_R_m)/delta_PW;
			(*this->m_matrix)( id(i,P), id(i,alpha_g) - total_var )	= (R_m_dalphaGasW - s_R_m)/delta_alphaGasW;
			(*this->m_matrix)( id(i,P), id(i,alpha_o) - total_var )	= (R_m_dalphaOilW - s_R_m)/delta_alphaOilW;
			(*this->m_matrix)( id(i,P), id(i,v)  - total_var )		= (R_m_dvW - s_R_m)/delta_vW;
			//CENTRAL			
			(*this->m_matrix)( id(i,P), id(i,P) )		= (R_m_dPP - s_R_m)/delta_PP;
			(*this->m_matrix)( id(i,P), id(i,alpha_g) )	= (R_m_dalphaGasP - s_R_m)/delta_alphaGasP;
			(*this->m_matrix)( id(i,P), id(i,alpha_o) )	= (R_m_dalphaOilP - s_R_m)/delta_alphaOilP;
			(*this->m_matrix)( id(i,P), id(i,v) )		= (R_m_dvP - s_R_m)/delta_vP;
			//EAST						
			(*this->m_matrix)( id(i,P), id(i,P) + total_var )		= (R_m_dPE - s_R_m)/delta_PE;
			(*this->m_matrix)( id(i,P), id(i,alpha_g) + total_var )	= (R_m_dalphaGasE - s_R_m)/delta_alphaGasE;
			(*this->m_matrix)( id(i,P), id(i,alpha_o) + total_var )	= (R_m_dalphaOilE - s_R_m)/delta_alphaOilE;
			//SOURCE
			(*this->m_source)[ id(i, P) ]	 = -s_R_m;		

			if( WITH_GAS ){
			//GAS
				s_R_g		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
										   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dPW		   = this->R_g(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dPP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dPE		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaGasW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaGasP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaGasE = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaOilW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaOilP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dalphaOilE = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_g_dvW		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], i, 'C');
				R_g_dvP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, i, 'C');

				//WEST			
				(*this->m_matrix)( id(i,alpha_g), id(i,P) - total_var )			= (R_g_dPW - s_R_g)/delta_PW;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_g) - total_var )	= (R_g_dalphaGasW - s_R_g)/delta_alphaGasW;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_o) - total_var )	= (R_g_dalphaOilW - s_R_g)/delta_alphaOilW;
				(*this->m_matrix)( id(i,alpha_g), id(i,v) - total_var )			= (R_g_dvW - s_R_g)/delta_vW;
				//CENTRAL			
				(*this->m_matrix)( id(i,alpha_g), id(i,P) )			= (R_g_dPP - s_R_g)/delta_PP;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_g) )	= (R_g_dalphaGasP - s_R_g)/delta_alphaGasP;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_o) )	= (R_g_dalphaOilP - s_R_g)/delta_alphaOilP;
				(*this->m_matrix)( id(i,alpha_g), id(i,v) )			= (R_g_dvP - s_R_g)/delta_vP;
				//EAST				
				(*this->m_matrix)( id(i,alpha_g), id(i,P) + total_var )			= (R_g_dPE - s_R_g)/delta_PE;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_g) + total_var )	= (R_g_dalphaGasE - s_R_g)/delta_alphaGasE;
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_o) + total_var )	= (R_g_dalphaOilE - s_R_g)/delta_alphaOilE;

				// SOURCE
				(*this->m_source)[ id(i, alpha_g) ] = -s_R_g;
			}
			else{
				(*this->m_matrix)( id(i,alpha_g), id(i,alpha_g) ) = 1.;
			}

			
			//Oil
				s_R_o		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
										   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dPW		   = this->R_o(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dPP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dPE		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaGasW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaGasP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaGasE = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaOilW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaOilP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dalphaOilE = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], i, 'C');
				R_o_dvW		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], i, 'C');
				R_o_dvP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, i, 'C');

				//WEST			
				(*this->m_matrix)( id(i,alpha_o), id(i,P) - total_var )			= (R_o_dPW - s_R_o)/delta_PW;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_g) - total_var )	= (R_o_dalphaGasW - s_R_o)/delta_alphaGasW;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_o) - total_var )	= (R_o_dalphaOilW - s_R_o)/delta_alphaOilW;
				(*this->m_matrix)( id(i,alpha_o), id(i,v) - total_var )			= (R_o_dvW - s_R_o)/delta_vW;
				//CENTRAL			
				(*this->m_matrix)( id(i,alpha_o), id(i,P) )			= (R_o_dPP - s_R_o)/delta_PP;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_g) )	= (R_o_dalphaGasP - s_R_o)/delta_alphaGasP;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_o) )	= (R_o_dalphaOilP - s_R_o)/delta_alphaOilP;
				(*this->m_matrix)( id(i,alpha_o), id(i,v) )			= (R_o_dvP - s_R_o)/delta_vP;
				//EAST				
				(*this->m_matrix)( id(i,alpha_o), id(i,P) + total_var )			= (R_o_dPE - s_R_o)/delta_PE;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_g) + total_var )	= (R_o_dalphaGasE - s_R_o)/delta_alphaGasE;
				(*this->m_matrix)( id(i,alpha_o), id(i,alpha_o) + total_var )	= (R_o_dalphaOilE - s_R_o)/delta_alphaOilE;

				// SOURCE
				(*this->m_source)[ id(i, alpha_o) ] = -s_R_o;
			

			if( WITH_MOMENTUM ){
                s_R_v		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dPW		   = this->R_v(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ],
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');              
                R_v_dPP		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ],
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dPE		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dPEE	   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ] + delta_PEE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaGasW = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');               
                R_v_dalphaGasP = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaGasE = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE, m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaGasEE= this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ] + delta_alphaGasEE, 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaOilW = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaOilP = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaOilE = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dalphaOilEE= this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ] + delta_alphaOilEE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dvW		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ],  
					m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], i, 'C');
                R_v_dvP		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, m_mean_velocity[ EAST ], i, 'C');
                R_v_dvE		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST+1 ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST+1 ], 
                    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST+1 ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ] + delta_vE, i, 'C');

				//WEST	
                (*this->m_matrix)( id(i,v), id(i,P) - total_var )		    = (R_v_dPW - s_R_v)/delta_PW;	
                (*this->m_matrix)( id(i,v), id(i,alpha_g) - total_var )     = (R_v_dalphaGasW - s_R_v)/delta_alphaGasW;
                (*this->m_matrix)( id(i,v), id(i,alpha_o) - total_var  )    = (R_v_dalphaOilW - s_R_v)/delta_alphaOilW;
				(*this->m_matrix)( id(i,v), id(i,v) - total_var )           = (R_v_dvW - s_R_v)/delta_vW;
				//CENTRAL			
				(*this->m_matrix)( id(i,v), id(i,P) )		 = (R_v_dPP - s_R_v)/delta_PP;	
				(*this->m_matrix)( id(i,v), id(i,alpha_g) )  = (R_v_dalphaGasP - s_R_v)/delta_alphaGasP;
				(*this->m_matrix)( id(i,v), id(i,alpha_o) )  = (R_v_dalphaOilP - s_R_v)/delta_alphaOilP;
				(*this->m_matrix)( id(i,v), id(i,v) )		 = (R_v_dvP - s_R_v)/delta_vP;
				//EAST				
				(*this->m_matrix)( id(i,v), id(i,P) + total_var )		 = (R_v_dPE - s_R_v)/delta_PE;	
				(*this->m_matrix)( id(i,v), id(i,alpha_g) + total_var )  = (R_v_dalphaGasE - s_R_v)/delta_alphaGasE;
				(*this->m_matrix)( id(i,v), id(i,alpha_o) + total_var )  = (R_v_dalphaOilE - s_R_v)/delta_alphaOilE;
				(*this->m_matrix)( id(i,v), id(i,v) + total_var )		 = (R_v_dvE - s_R_v)/delta_vE;
                //EEAST				
                (*this->m_matrix)( id(i,v), id(i,P) + 2*total_var )		   = (R_v_dPEE - s_R_v)/delta_PEE;	
                (*this->m_matrix)( id(i,v), id(i,alpha_g) + 2*total_var )  = (R_v_dalphaGasEE - s_R_v)/delta_alphaGasEE;
                (*this->m_matrix)( id(i,v), id(i,alpha_o) + 2*total_var )  = (R_v_dalphaOilEE - s_R_v)/delta_alphaOilEE;
				//SOURCE
				(*this->m_source)[ id(i, v) ]	 = -s_R_v;
			}
			else{
				(*this->m_matrix)( id(i,v), id(i,P) + total_var ) = 1.;
			}
		
						
			++WEST;	++CENT;	++EAST;
		}

		delta_PP		= m_delta[ P ]*m_pressure[ CENT ];
		delta_alphaGasP	= m_gas_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ CENT ] : m_delta[ alpha_g ];
		delta_alphaOilP	= m_oil_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ CENT ] : m_delta[ alpha_o ];
		delta_vP		= abs(m_mean_velocity[ CENT ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ CENT ] : m_delta[ v ];					
		delta_PE		= m_delta[ P ]*m_pressure[ EAST ];
		delta_alphaGasE	= m_gas_vol_frac[ EAST ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ EAST ] : m_delta[ alpha_g ];
		delta_alphaOilE	= m_oil_vol_frac[ EAST ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ EAST ] : m_delta[ alpha_o ];
		delta_vE		= abs(m_mean_velocity[ EAST ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ EAST ] : m_delta[ v ];			
		delta_PW		= m_delta[ P ]*m_pressure[ WEST ];
		delta_alphaGasW	= m_gas_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ WEST ] : m_delta[ alpha_g ];
		delta_alphaOilW	= m_oil_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ WEST ] : m_delta[ alpha_o ];
		delta_vW		= abs(m_mean_velocity[ WEST ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ WEST ] : m_delta[ v ];
		
		
	//MIXTURE  
		s_R_m		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dPW		   = this->R_m(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dPP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dPE		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaGasW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaGasP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaGasE = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaOilW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaOilP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dalphaOilE = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dvW		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST - 1, 'C');
		R_m_dvP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST - 1, 'C');
		//WEST			
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,P) - total_var )			= (R_m_dPW - s_R_m)/delta_PW;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_g) - total_var )	= (R_m_dalphaGasW - s_R_m)/delta_alphaGasW;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_o) - total_var )	= (R_m_dalphaOilW - s_R_m)/delta_alphaOilW;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,v)  - total_var )		= (R_m_dvW - s_R_m)/delta_vW;
		//CENTRAL			
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,P) )			= (R_m_dPP - s_R_m)/delta_PP;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_g) )	= (R_m_dalphaGasP - s_R_m)/delta_alphaGasP;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_o) )	= (R_m_dalphaOilP - s_R_m)/delta_alphaOilP;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,v) )			= (R_m_dvP - s_R_m)/delta_vP;
		//EAST						
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,P) + total_var )			= (R_m_dPE - s_R_m)/delta_PE;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_g) + total_var )	= (R_m_dalphaGasE - s_R_m)/delta_alphaGasE;
		(*this->m_matrix)( id(LAST - 1,P), id(LAST - 1,alpha_o) + total_var )	= (R_m_dalphaOilE - s_R_m)/delta_alphaOilE;
		//SOURCE
		(*this->m_source)[ id(LAST - 1, P) ]	 = -s_R_m;		

		
		if( WITH_GAS ){
		//GAS	
			s_R_g		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
										m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dPW		   = this->R_g(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dPP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dPE		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaGasW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaGasP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaGasE = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaOilW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaOilP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dalphaOilE = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dvW		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_g_dvP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST - 1, 'C');

			//WEST			
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,P) - total_var )			= (R_g_dPW - s_R_g)/delta_PW;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_g) - total_var )		= (R_g_dalphaGasW - s_R_g)/delta_alphaGasW;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_o) - total_var )		= (R_g_dalphaOilW - s_R_g)/delta_alphaOilW;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,v) - total_var )			= (R_g_dvW - s_R_g)/delta_vW;
			//CENTRAL			
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,P) )			= (R_g_dPP - s_R_g)/delta_PP;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_g) )		= (R_g_dalphaGasP - s_R_g)/delta_alphaGasP;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_o) )		= (R_g_dalphaOilP - s_R_g)/delta_alphaOilP;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,v) )			= (R_g_dvP - s_R_g)/delta_vP;
			//EAST				
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,P) + total_var )			= (R_g_dPE - s_R_g)/delta_PE;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_g) + total_var )		= (R_g_dalphaGasE - s_R_g)/delta_alphaGasE;
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_o) + total_var )		= (R_g_dalphaOilE - s_R_g)/delta_alphaOilE;

			// SOURCE
			(*this->m_source)[ id(LAST - 1, alpha_g) ] = -s_R_g;
			
		}
		else{
			(*this->m_matrix)( id(LAST - 1,alpha_g), id(LAST - 1,alpha_g) ) = 1.;
		}


		//Oil
			s_R_o		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
										m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dPW		   = this->R_o(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dPP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dPE		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaGasW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaGasP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaGasE = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE,  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaOilW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaOilP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dalphaOilE = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dvW		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST - 1, 'C');
			R_o_dvP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ],  
									    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST - 1, 'C');

			//WEST			
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,P) - total_var )			= (R_o_dPW - s_R_o)/delta_PW;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_g) - total_var )		= (R_o_dalphaGasW - s_R_o)/delta_alphaGasW;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_o) - total_var )		= (R_o_dalphaOilW - s_R_o)/delta_alphaOilW;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,v) - total_var )			= (R_o_dvW - s_R_o)/delta_vW;
			//CENTRAL			
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,P) )			= (R_o_dPP - s_R_o)/delta_PP;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_g) )		= (R_o_dalphaGasP - s_R_o)/delta_alphaGasP;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_o) )		= (R_o_dalphaOilP - s_R_o)/delta_alphaOilP;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,v) )			= (R_o_dvP - s_R_o)/delta_vP;
			//EAST				
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,P) + total_var )			= (R_o_dPE - s_R_o)/delta_PE;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_g) + total_var )		= (R_o_dalphaGasE - s_R_o)/delta_alphaGasE;
			(*this->m_matrix)( id(LAST - 1,alpha_o), id(LAST - 1,alpha_o) + total_var )		= (R_o_dalphaOilE - s_R_o)/delta_alphaOilE;

			// SOURCE
			(*this->m_source)[ id(LAST - 1, alpha_o) ] = -s_R_o;

		if( WITH_MOMENTUM ){				
			//MOMENTUM

            s_R_v		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dPW		   = this->R_v(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ],
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');              
            R_v_dPP		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ],
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dPE		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ] + delta_PE, m_pressure[ EAST ] + delta_PE, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dalphaGasW = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');               
            R_v_dalphaGasP = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dalphaGasE = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ] + delta_alphaGasE, m_gas_vol_frac[ EAST ] + delta_alphaGasE, 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dalphaOilW = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dalphaOilP = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dalphaOilE = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_oil_vol_frac[ EAST ] + delta_alphaOilE, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dvW		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ],  
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dvP		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, m_mean_velocity[ EAST ], LAST - 1, 'L');
            R_v_dvE		   = this->R_v(m_pressure[ WEST ], m_pressure[ CENT ], m_pressure[ EAST ], m_pressure[ EAST ], m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], m_gas_vol_frac[ EAST ], m_gas_vol_frac[ EAST ], 
                m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], m_oil_vol_frac[ EAST ], m_oil_vol_frac[ EAST ], m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], m_mean_velocity[ EAST ] + delta_vE, LAST - 1, 'L');



				//WEST
                (*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,P) - total_var )		    = (R_v_dPW - s_R_v)/delta_PW;	
                (*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_g) - total_var )   = (R_v_dalphaGasW - s_R_v)/delta_alphaGasW;
                (*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_o) - total_var  )  = (R_v_dalphaOilW - s_R_v)/delta_alphaOilW;
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,v) - total_var )         = (R_v_dvW - s_R_v)/delta_vW;
				//CENTRAL			
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,P) )			= (R_v_dPP - s_R_v)/delta_PP;	
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_g) )   = (R_v_dalphaGasP - s_R_v)/delta_alphaGasP;
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_o) )   = (R_v_dalphaOilP - s_R_v)/delta_alphaOilP;
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,v) )			= (R_v_dvP - s_R_v)/delta_vP;
				//EAST				
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,P) + total_var )		   = (R_v_dPE - s_R_v)/delta_PE;	
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_g) + total_var )  = (R_v_dalphaGasE - s_R_v)/delta_alphaGasE;
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,alpha_o) + total_var )  = (R_v_dalphaOilE - s_R_v)/delta_alphaOilE;
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,v) + total_var )		   = (R_v_dvE - s_R_v)/delta_vE;
				//SOURCE
				(*this->m_source)[ id(LAST - 1, v) ]	 = -s_R_v;				
			}
			else{
				(*this->m_matrix)( id(LAST - 1,v), id(LAST - 1,P) + total_var ) = 1.;
			}
	
	
		

		++WEST;	++CENT;

		delta_PP		= m_delta[ P ]*m_pressure[ CENT ];
		delta_alphaGasP	= m_gas_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ CENT ] : m_delta[ alpha_g ];
		delta_alphaOilP	= m_oil_vol_frac[ CENT ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ CENT ] : m_delta[ alpha_o ];
		delta_vP		= abs(m_mean_velocity[ CENT ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ CENT ] : m_delta[ v ];					
		
			
		delta_PW		= m_delta[ P ]*m_pressure[ WEST ];
		delta_alphaGasW	= m_gas_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_g ]*m_gas_vol_frac[ WEST ] : m_delta[ alpha_g ];
		delta_alphaOilW	= m_oil_vol_frac[ WEST ] > 1e-12 ? m_delta[ alpha_o ]*m_oil_vol_frac[ WEST ] : m_delta[ alpha_o ];
		delta_vW		= abs(m_mean_velocity[ WEST ]) > 1e-12 ? m_delta[ v ]*m_mean_velocity[ WEST ] : m_delta[ v ];
		

		//MIXTURE				
			s_R_m		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dPW		   = this->R_m(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dPP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dalphaGasW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dalphaGasP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dalphaOilW = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dalphaOilP = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dvW		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
								       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST, 'L');
			R_m_dvP		   = this->R_m(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST, 'L');
			//WEST			
			(*this->m_matrix)( id(LAST,P), id(LAST,P) - total_var )			= (R_m_dPW - s_R_m)/delta_PW;
			(*this->m_matrix)( id(LAST,P), id(LAST,alpha_g) - total_var )	= (R_m_dalphaGasW - s_R_m)/delta_alphaGasW;
			(*this->m_matrix)( id(LAST,P), id(LAST,alpha_o) - total_var )	= (R_m_dalphaOilW - s_R_m)/delta_alphaOilW;
			(*this->m_matrix)( id(LAST,P), id(LAST,v)  - total_var )		= (R_m_dvW - s_R_m)/delta_vW;
			//CENTRAL					
			(*this->m_matrix)( id(LAST,P), id(LAST,P) )			= (R_m_dPP - s_R_m)/delta_PP;
			(*this->m_matrix)( id(LAST,P), id(LAST,alpha_g) )	= (R_m_dalphaGasP - s_R_m)/delta_alphaGasP;
			(*this->m_matrix)( id(LAST,P), id(LAST,alpha_o) )	= (R_m_dalphaOilP - s_R_m)/delta_alphaOilP;
			(*this->m_matrix)( id(LAST,P), id(LAST,v) )			= (R_m_dvP - s_R_m)/delta_vP;
			// SOURCE
			(*this->m_source)[ id(LAST, P) ]	 = -s_R_m;			
			
		
		if( WITH_GAS ){
		//GAS	
			s_R_g		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dPW		   = this->R_g(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dPP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dalphaGasW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dalphaGasP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dalphaOilW = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dalphaOilP = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dvW		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
								       m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST, 'L');
			R_g_dvP		   = this->R_g(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									   m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST, 'L');

			
			//WEST		
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,P) - total_var )	    = (R_g_dPW - s_R_g)/delta_PW;
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,alpha_g) - total_var ) = (R_g_dalphaGasW - s_R_g)/delta_alphaGasW;
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,alpha_o) - total_var ) = (R_g_dalphaOilW - s_R_g)/delta_alphaOilW;
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,v) - total_var )	    = (R_g_dvW - s_R_g)/delta_vW;
			//CENTRAL		
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,P) )	    = (R_g_dPP - s_R_g)/delta_PP;
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,alpha_g) ) = (R_g_dalphaGasP - s_R_g)/delta_alphaGasP;	
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,alpha_o) ) = (R_g_dalphaOilP - s_R_g)/delta_alphaOilP;	
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,v) )	    = (R_g_dvP - s_R_g)/delta_vP;
			// SOURCE
			(*this->m_source)[ id(LAST, alpha_g) ] = -s_R_g;			

		}
		else{
			(*this->m_matrix)( id(LAST,alpha_g), id(LAST,alpha_g) ) = 1.;
		}		
       
		s_R_o		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dPW		   = this->R_o(m_pressure[ WEST ] + delta_PW, m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dPP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ] + delta_PP, 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dalphaGasW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ] + delta_alphaGasW, m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dalphaGasP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ] + delta_alphaGasP, 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dalphaOilW = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ] + delta_alphaOilW, m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dalphaOilP = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ] + delta_alphaOilP, 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dvW		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
								    m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ] + delta_vW, m_mean_velocity[ CENT ], LAST, 'L');
		R_o_dvP		   = this->R_o(m_pressure[ WEST ], m_pressure[ CENT ], 0, m_gas_vol_frac[ WEST ], m_gas_vol_frac[ CENT ], 0,  
									m_oil_vol_frac[ WEST ], m_oil_vol_frac[ CENT ], 0, m_mean_velocity[ WEST ], m_mean_velocity[ CENT ] + delta_vP, LAST, 'L');

			
		//WEST		
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,P) - total_var )	    = (R_o_dPW - s_R_o)/delta_PW;
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,alpha_g) - total_var ) = (R_o_dalphaGasW - s_R_o)/delta_alphaGasW;
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,alpha_o) - total_var ) = (R_o_dalphaOilW - s_R_o)/delta_alphaOilW;
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,v) - total_var )	    = (R_o_dvW - s_R_o)/delta_vW;
		//CENTRAL		
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,P) )	    = (R_o_dPP - s_R_o)/delta_PP;
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,alpha_g) ) = (R_o_dalphaGasP - s_R_o)/delta_alphaGasP;	
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,alpha_o) ) = (R_o_dalphaOilP - s_R_o)/delta_alphaOilP;	
		(*this->m_matrix)( id(LAST,alpha_o), id(LAST,v) )	    = (R_o_dvP - s_R_o)/delta_vP;
		// SOURCE           
		(*this->m_source)[ id(LAST, alpha_o) ] = -s_R_o;			



		// MOMENTUM					
			(*this->m_matrix)( id(LAST,v), id(LAST,v) )	   = 1.;
		
		
		
			
         /* std::ofstream matrixm("WellData\\Matrixmatlab.dat");
            matrixm << std::setprecision(3);
			matrixm << "[ ";
			for( unsigned i =0 ; i<m_matrix->nrows(); ++i )
			{				
				for( unsigned j =0 ; j < m_matrix->ncols(); ++j )
				{
					matrixm << "  " << (*m_matrix)(i,j) << "  ";
				}				
				matrixm << "\n";				
			}
			matrixm << " ] \n\n";
			matrixm << "[ ";
			for( unsigned i =0 ; i<m_matrix->nrows(); ++i )
			{				
				matrixm << (*m_source)[ i ] << "\n";
			}
			matrixm << " ] ";*/

	}

    real_type DriftFluxWell::calculate_new_delta_t_size_diverged_solution(real_type delta_t_old){
        return 0.5*delta_t_old;
    }

    real_type DriftFluxWell::calculate_new_delta_t_size_converged_solution(real_type delta_t_old)
    {
        real_type delta_t_S;

        real_type delta_S_max;       

        delta_S_max = 0.0;      

        for( int i = 0; i < number_of_nodes(); ++i )
        {   

            real_type vol_frac_gas_variation = std::fabs( m_gas_vol_frac[i] - m_gas_vol_frac_old[i] );

            real_type vol_frac_oil_variation= std::fabs( m_oil_vol_frac[i] - m_oil_vol_frac_old[i] );

            real_type vol_frac_variation = std::max(vol_frac_gas_variation, vol_frac_oil_variation);


            if (vol_frac_variation > delta_S_max) {
                delta_S_max = vol_frac_variation;
            }            
        }

        delta_t_S = delta_t_old * 0.05 / delta_S_max;
       
        real_type delta_t = std::min(  m_max_delta_t, std::max( 1e-5, delta_t_S));
        if( m_current_time + delta_t > m_final_time){
           delta_t =  m_final_time - m_current_time;
        }
        return delta_t;

    }

    void DriftFluxWell::restore_initial_guess(){
        for( uint_type i = 0; i < number_of_nodes()-1; ++i )
        {
            this->m_pressure[ i ]		= m_pressure_old[ i ];
            this->m_gas_vol_frac[ i ]	= m_gas_vol_frac_old[ i ];
            this->m_oil_vol_frac[ i ]	= m_oil_vol_frac_old[ i ];
            this->m_mean_velocity[ i ]	= m_mean_velocity_old[ i ];
                                           
            (*this->m_variables)[ total_var*i ]             = 0.0;
            (*this->m_variables)[ total_var*i + alpha_g ]   = 0.0;
            (*this->m_variables)[ total_var*i + alpha_o ]   = 0.0;
            (*this->m_variables)[ total_var*i + v ]         = 0.0;

          /*  this->m_water_vol_frac[ i ] = 1.0 - (m_gas_vol_frac[ i ] + m_oil_vol_frac[ i ]);



            real_type KSI	   = this->ksi( m_mean_velocity[ i ] );
            real_type c_P	   = (0.5+KSI)*m_pressure[ i ] + (0.5-KSI)*m_pressure[ i+1 ];
            real_type c_alphaG = (0.5+KSI)*m_gas_vol_frac[ i ] + (0.5-KSI)*m_gas_vol_frac[ i+1 ];
            real_type c_alphaO = (0.5+KSI)*m_oil_vol_frac[ i ] + (0.5-KSI)*m_oil_vol_frac[ i+1 ];
            real_type c_alphaW = 1.0 - (c_alphaG+ c_alphaO);

            real_type rho  = this->mean_density  ( c_alphaO, c_alphaW, c_alphaG, c_P );
            real_type rhoL = this->liquid_density(c_alphaO, c_alphaW, c_P);
            real_type rhoG = this->gas_density   ( c_P );
            real_type rhoO = this->gas_density   ( c_P );
            real_type rhoW = this->gas_density   ( c_P );

            real_type mod_Vgj = this->mod_v_drift_flux( m_mean_velocity[ i ], 0.5*(m_gas_vol_frac[ i ]+m_gas_vol_frac[ i+1 ]), 0.5*(m_oil_vol_frac[ i ]+m_oil_vol_frac[ i+1 ]), 0.5*(m_water_vol_frac[ i ] + m_water_vol_frac[ i+1 ]), 0.5*(m_pressure[ i ]+m_pressure[ i+1 ]) );
            real_type mod_Vow = this->mod_v_drift_flux_ow( m_mean_velocity[ i ], 0.5*(m_gas_vol_frac[ i ]+m_gas_vol_frac[ i+1 ]), 0.5*(m_oil_vol_frac[ i ]+m_oil_vol_frac[ i+1 ]), 0.5*(m_water_vol_frac[ i ] + m_water_vol_frac[ i+1 ]), 0.5*(m_pressure[ i ]+m_pressure[ i+1 ]) );


            this->m_gas_velocity[ i ] = m_mean_velocity[ i ] + rhoL/rho*mod_Vgj;			
            real_type liquid_velocity = m_mean_velocity[ i ] - c_alphaG/(1 - c_alphaG + 1.0e-20)*rhoG/rho*mod_Vgj;	

            this->m_oil_velocity[ i ]   = liquid_velocity + rhoW/rhoL*mod_Vow;
            this->m_water_velocity[ i ] = liquid_velocity - (c_alphaO/c_alphaW)*(rhoO/rhoL)*mod_Vow;		*/

        }

        unsigned i = number_of_nodes()-1;

        this->m_pressure[ i ]		= m_pressure_old[ i ];
        this->m_gas_vol_frac[ i ]	= m_gas_vol_frac_old[ i ];
        this->m_oil_vol_frac[ i ]	= m_oil_vol_frac_old[ i ];
        this->m_mean_velocity[ i ]	= m_mean_velocity_old[ i ];

        (*this->m_variables)[ total_var*i ]             = 0.0;
        (*this->m_variables)[ total_var*i + alpha_g ]   = 0.0;
        (*this->m_variables)[ total_var*i + alpha_o ]   = 0.0;
        (*this->m_variables)[ total_var*i + v ]         = 0.0;

      /*  this->m_water_vol_frac[ i ] = 1.0 - (m_gas_vol_frac[ i ] + m_oil_vol_frac[ i ]);
        this->m_gas_velocity[ i ]	= m_mean_velocity[ i ];
        this->m_oil_velocity[ i ]   = m_mean_velocity[ i ];
        this->m_water_velocity[ i ] = m_mean_velocity[ i ];*/



        // Last volume fraction must be equal
        this->m_gas_vol_frac[ 0 ] = this->m_gas_vol_frac[ 1 ];
        this->m_oil_vol_frac[ 0 ] = this->m_oil_vol_frac[ 1 ];

    }


	void DriftFluxWell::update_variables()
	{
		
		for( uint_type i = 0; i < number_of_nodes()-1; ++i )
		{
			this->m_pressure[ i ]		+= (*this->m_variables)[ total_var*i ];
			this->m_gas_vol_frac[ i ]	+= (*this->m_variables)[ total_var*i + alpha_g ];
			this->m_oil_vol_frac[ i ]	+= (*this->m_variables)[ total_var*i + alpha_o ];
			this->m_mean_velocity[ i ]	+= (*this->m_variables)[ total_var*i + v ];

           // (*this->m_variables)[ total_var*i ]             = 0.0;
           // (*this->m_variables)[ total_var*i + alpha_g ]   = 0.0;
           // (*this->m_variables)[ total_var*i + alpha_o ]   = 0.0;
          //  (*this->m_variables)[ total_var*i + v ]         = 0.0;

			this->m_water_vol_frac[ i ] = 1.0 - (m_gas_vol_frac[ i ] + m_oil_vol_frac[ i ]);

						
			
			real_type KSI	   = this->ksi( m_mean_velocity[ i ] );
			real_type c_P	   = (0.5+KSI)*m_pressure[ i ] + (0.5-KSI)*m_pressure[ i+1 ];
			real_type c_alphaG = (0.5+KSI)*m_gas_vol_frac[ i ] + (0.5-KSI)*m_gas_vol_frac[ i+1 ];
			real_type c_alphaO = (0.5+KSI)*m_oil_vol_frac[ i ] + (0.5-KSI)*m_oil_vol_frac[ i+1 ];
			real_type c_alphaW = 1.0 - (c_alphaG+ c_alphaO);

			real_type rho  = this->mean_density  ( c_alphaO, c_alphaW, c_alphaG, c_P );
			real_type rhoL = this->liquid_density(c_alphaO, c_alphaW, c_P);
			real_type rhoG = this->gas_density   ( c_P );
			real_type rhoO = this->gas_density   ( c_P );
			real_type rhoW = this->gas_density   ( c_P );
			
			real_type mod_Vgj = this->mod_v_drift_flux( m_mean_velocity[ i ], 0.5*(m_gas_vol_frac[ i ]+m_gas_vol_frac[ i+1 ]), 0.5*(m_oil_vol_frac[ i ]+m_oil_vol_frac[ i+1 ]), 0.5*(m_water_vol_frac[ i ] + m_water_vol_frac[ i+1 ]), 0.5*(m_pressure[ i ]+m_pressure[ i+1 ]) );
			real_type mod_Vow = this->mod_v_drift_flux_ow( m_mean_velocity[ i ], 0.5*(m_gas_vol_frac[ i ]+m_gas_vol_frac[ i+1 ]), 0.5*(m_oil_vol_frac[ i ]+m_oil_vol_frac[ i+1 ]), 0.5*(m_water_vol_frac[ i ] + m_water_vol_frac[ i+1 ]), 0.5*(m_pressure[ i ]+m_pressure[ i+1 ]) );

			
			this->m_gas_velocity[ i ] = m_mean_velocity[ i ] + rhoL/rho*mod_Vgj;			
			real_type liquid_velocity = m_mean_velocity[ i ] - c_alphaG/(1 - c_alphaG + 1.0e-20)*rhoG/rho*mod_Vgj;	
			
			this->m_oil_velocity[ i ]   = liquid_velocity + rhoW/rhoL*mod_Vow;
			this->m_water_velocity[ i ] = liquid_velocity - (c_alphaO/(c_alphaW + 1.0e-20))*(rhoO/rhoL)*mod_Vow;		
			
		}

		unsigned i = number_of_nodes()-1;


		this->m_pressure[ i ]		+= (*this->m_variables)[ total_var*i ];
		this->m_gas_vol_frac[ i ]	+= (*this->m_variables)[ total_var*i + alpha_g ];
		this->m_oil_vol_frac[ i ]	+= (*this->m_variables)[ total_var*i + alpha_o ];
		this->m_mean_velocity[ i ]	+= (*this->m_variables)[ total_var*i + v ];

      //  (*this->m_variables)[ total_var*i ]             = 0.0;
      //  (*this->m_variables)[ total_var*i + alpha_g ]   = 0.0;
      //  (*this->m_variables)[ total_var*i + alpha_o ]   = 0.0;
     //   (*this->m_variables)[ total_var*i + v ]         = 0.0;

		this->m_water_vol_frac[ i ] = 1.0 - (m_gas_vol_frac[ i ] + m_oil_vol_frac[ i ]);
		this->m_gas_velocity[ i ]	= m_mean_velocity[ i ];
		this->m_oil_velocity[ i ]   = m_mean_velocity[ i ];
		this->m_water_velocity[ i ] = m_mean_velocity[ i ];

	

		// Last volume fraction must be equal
		this->m_gas_vol_frac[ 0 ] = this->m_gas_vol_frac[ 1 ];
		this->m_oil_vol_frac[ 0 ] = this->m_oil_vol_frac[ 1 ];
		
	}

	

	void DriftFluxWell::set_boundary_velocity( real_type p_velocity ){
		this->m_mean_velocity[ number_of_nodes() - 1 ] = p_velocity;
	}

    std::string make_filename( const std::string& basename, int index, const std::string& ext )
	{
        std::ostringstream result;
		result << basename << index << ext;
		return result.str();
	}

    void DriftFluxWell::update_variables_for_new_timestep(){
        for( uint_type i = 0; i < number_of_nodes(); ++i){

            m_gas_vol_frac_old[ i ]   = m_gas_vol_frac[ i ];			
            m_oil_vol_frac_old[ i ]   = m_oil_vol_frac[ i ];
            m_water_vol_frac_old[ i ] = m_water_vol_frac[ i ];
            m_pressure_old[ i ]		  = m_pressure[ i ];		
            m_mean_velocity_old[ i ]  = m_mean_velocity[ i ];

            m_water_flow[i]->calculate_value_at_time(m_current_time);
            m_oil_flow[i]->calculate_value_at_time(m_current_time);
            m_gas_flow[i]->calculate_value_at_time(m_current_time);               
        }    
    }                     

	void DriftFluxWell::solve()
	{
        bool check_time = false;
		m_current_time = 0;

		uint_type FINAL_TIMESTEP = this->m_FINAL_TIMESTEP;

		this->update_variables_for_new_timestep();
		
		this->set_bottom_pressure( m_HEEL_PRESSURE ); // Pressure at heel is set

        bool log_output_is_active = true;

        std::ofstream log_results_file;
        if(log_output_is_active)
        {                           
            log_results_file.open( "..\\WellData\\log_results.txt" );
            log_results_file << "Current time" << "\t"
                << "Newton Iter"	<< "\t"
                << "Final norm"     << "\n";
            log_results_file << std::setprecision(10);
        }
		
		for(uint_type TIMESTEP = 0; TIMESTEP < FINAL_TIMESTEP; ++TIMESTEP){
		//	(*m_water_flow)[ number_of_nodes()-1 ] = (TIMESTEP+1)*dt() > 10.0 ? 1.5 : 1.5*(TIMESTEP+1)*dt()/10;
		//	(*m_oil_flow)[ number_of_nodes()-1 ]   = (TIMESTEP+1)*dt() > 10.0 ? 1.5 : 1.5*(TIMESTEP+1)*dt()/10;
		//	(*m_gas_flow)[ number_of_nodes()-1 ]   = (TIMESTEP+1)*dt() > 10.0 ? 0.02 : 0.02*(TIMESTEP+1)*dt()/10;
			

			if(TIMESTEP){ 
				// Only enters loop for TIMESTEP > 0                    
                set_dt( calculate_new_delta_t_size_converged_solution( dt() ) );                 
				//for( uint_type i = 0; i < number_of_nodes(); ++i){
				//	//real_type dalpha    = m_gas_vol_frac[ i ]-m_gas_vol_frac_old[ i ];
				//	//real_type dpressure = m_pressure[ i ]-m_pressure_old[ i ];
				//	//real_type dvelocity = m_mean_velocity[ i ]-m_mean_velocity_old[ i ];				
				//	/*m_gas_vol_frac_old[ i ] = m_gas_vol_frac[ i ];
				//	m_gas_vol_frac[ i ] =  m_gas_vol_frac[ i ] + dalpha;

				//	m_oil_vol_frac_old[ i ] = m_oil_vol_frac[ i ];
				//	m_water_vol_frac_old[ i ] = m_water_vol_frac[ i ];

				//	m_pressure_old[ i ] = m_pressure[ i ];
				//	m_pressure[ i ]     = m_pressure[ i ] + 0*dpressure;

				//	m_mean_velocity_old[ i ] = m_mean_velocity[ i ];
				//	m_mean_velocity[ i ]     = m_mean_velocity[ i ] + 0*dvelocity;*/

				//	m_gas_vol_frac_old[ i ] = m_gas_vol_frac[ i ];	
				//	m_oil_vol_frac_old[ i ] = m_oil_vol_frac[ i ];
				//	//m_water_vol_frac_old[ i ] = m_water_vol_frac[ i ];
				//	m_pressure_old[ i ] = m_pressure[ i ];	
				//	m_mean_velocity_old[ i ] = m_mean_velocity[ i ];
				//
				//}

                this->update_variables_for_new_timestep();

                
			}
			
			uint_type r = 0;
            std::queue<real_type> norm_history;
			real_type norma;
			do
			{
				static Timer timer;
                //timer.enable_print_time();

                timer.start();                   
				this->compute_Jacobian();
                timer.stop();
                timer.print("\njacobian time = ");

				timer.start();
				GMRES_Solve( *m_matrix, *m_variables, *m_source );
				
				timer.stop();
                timer.print("\nsolver time = ");
								
				this->update_variables();			
				
                std::cout << std::setprecision(10);
				
				norma = itl::two_norm(*m_source);				
                std::cout << "\n----norma residuo: " << norma << "-----\n";
				
                norm_history.push(norma);
                if(r > 3){
                    norm_history.pop();
                }
				/*for( uint_type i = 0; i < number_of_nodes()-1; ++i ){
					cout << setprecision(10);				
					cout << "pressure[ "<< i <<" ] = " << this->m_pressure[ i ] <<"\t";
					if( i == 0)
						cout <<"\t";
					cout << "Gas_vol_frac[ "<< i <<" ] = " << this->m_gas_vol_frac[ i ] <<"\t";
					cout << "Oil_vol_frac[ "<< i <<" ] = " << this->m_oil_vol_frac[ i ] <<"\t";	
					if( i == number_of_nodes() - 1)
						cout <<"\t";				
					cout << "mean_velocity[ "<< i <<" ] = " << this->m_mean_velocity[ i ] << "\n";			

				}*/
                //if( norm_history.front() < norm_history.back() && r > 3)    m_convergence_status = true;
                m_convergence_status = false;
                if(norma > this->NEWTON_CRIT && r > 50 || m_convergence_status){
                    set_dt( calculate_new_delta_t_size_diverged_solution( dt() ) ); 
                    std::cout << "\n********* Breaking timestep = " << dt();
                    int r_inner = 0;
                    real_type new_norm = 0.0;
                    std::queue<real_type> new_norm_history;
                    restore_initial_guess();
                    do{                                
                        // Restart solution with half timestep 
                        timer.start();
                        this->compute_Jacobian();							
                        timer.stop();
                        timer.print("\njacobian time = ");

                        timer.start();
                        GMRES_Solve( *m_matrix, *m_variables, *m_source );

                        timer.stop();
                        timer.print("\nsolver time = ");
                       
                        this->update_variables();			

                        std::cout << std::setprecision(10);

                        new_norm = itl::two_norm(*m_source);				
                        std::cout << "\n----norma residuo: " << new_norm << "-----\n";

                        new_norm_history.push(new_norm);
                        if(r_inner > 2){
                            new_norm_history.pop();
                        }
                       // if( new_norm > norma && r_inner > 10) break;
                        if( new_norm_history.front() < new_norm_history.back() && r_inner > 2) {
                           // restore_initial_guess();
                            break;
                        }
                        ++r_inner;                                             
                    }while(new_norm > this->NEWTON_CRIT && r_inner < 30);
                    std::cout << "\n********* Returning to normal loop";
                }

				
				++r;
			}while(norma > this->NEWTON_CRIT && r < 1000);

            m_current_time += this->dt();
            std::cout << "TIME: " << m_current_time << " seconds\n\n\n"; 
                                           
            if(log_output_is_active)
            {          		
                log_results_file << m_current_time	<< "\t"
                    << r	    << "\t"
                    << norma	<< "\n";                                 
            }

            /*double transient_norm = 0;
            double transient_norm_P = 0;
            double transient_norm_alphaG = 0;
            double transient_norm_alphaO = 0;
            double transient_norm_v = 0;
            for(unsigned i = 0; i < m_nnodes; ++i){
                transient_norm_P += (m_pressure_old[i] - m_pressure[i])*(m_pressure_old[i] - m_pressure[i]);
                transient_norm_alphaG += (m_gas_vol_frac_old[i] - m_gas_vol_frac[i])*(m_gas_vol_frac_old[i] - m_gas_vol_frac[i]);	
                transient_norm_alphaO += (m_oil_vol_frac_old[i] - m_oil_vol_frac[i])*(m_oil_vol_frac_old[i] - m_oil_vol_frac[i]);
                transient_norm_v += (m_mean_velocity_old[i] - m_mean_velocity[i])*(m_mean_velocity_old[i] - m_mean_velocity[i]);
            }
            transient_norm_P = sqrt(transient_norm_P);
            transient_norm_alphaG = sqrt(transient_norm_alphaG);
            transient_norm_alphaO = sqrt(transient_norm_alphaO);
            transient_norm_v = sqrt(transient_norm_v);
            transient_norm = 0.0*transient_norm_P + transient_norm_alphaG + transient_norm_alphaO + 0.0*transient_norm_v;
            transient_norm = transient_norm/dt();*/

            //if(transient_norm < 1e-6 && TIMESTEP > 0 || TIMESTEP == FINAL_TIMESTEP-1 || abs(m_current_time - m_final_time) < 1.0e-8 )
            if(TIMESTEP == FINAL_TIMESTEP-1 || abs(m_current_time - m_final_time) < 1.0e-8 )
            {
                std::ofstream results_file;
                //results_file.open( make_filename( "results", TIMESTEP, ".dat" ).c_str() );
                results_file.open( "..\\WellData\\results.txt" );
                results_file << "Pressure [Pa]"	<< "\t"
                    << "Gas Volume Fraction [-]"    << "\t"
                    << "Oil Volume Fraction [-]"	<< "\t"
                    << "Water Volume Fraction [-]"	<< "\t"
                    << "Mixture Velocity [m/s]"	    << "\t"
                    << "Gas Velocity [m/s]"	        << "\t"
                    << "Oil Velocity [m/s]"	        << "\t"
                    << "Water Velocity [m/s]"	    << "\n";
                for( uint_type i = 0; i < number_of_nodes(); ++i )
                {
                    results_file << std::setprecision(10);				
                    results_file << this->m_pressure      [ i ]	<< "\t"
                        << this->m_gas_vol_frac  [ i ]	<< "\t"
                        << this->m_oil_vol_frac  [ i ]	<< "\t"
                        << this->m_water_vol_frac[ i ]	<< "\t"
                        << this->m_mean_velocity [ i ]	<< "\t"
                        << this->m_gas_velocity  [ i ]	<< "\t"
                        << this->m_oil_velocity  [ i ]	<< "\t"
                        << this->m_water_velocity[ i ]	<< "\n";												
                }
                results_file.close();
                break;
            }
		}

        if(log_output_is_active)
        {  
            log_results_file.close();                     
        }
	}




    void DriftFluxWell::solve(vector_type& p_pressure)
	{
		m_current_time = 0;
        static int STEPS = 0;
        m_total_production[OilPhase]    = 0.0;
        m_total_production[GasPhase]    = 0.0;
        m_total_production[WaterPhase]  = 0.0;

		uint_type FINAL_TIMESTEP = this->m_FINAL_TIMESTEP;
		for( uint_type i = 0; i < number_of_nodes(); ++i){
			
			m_gas_vol_frac_old[ i ]   = m_gas_vol_frac[ i ];			
			m_oil_vol_frac_old[ i ]   = m_oil_vol_frac[ i ];
			m_water_vol_frac_old[ i ] = m_water_vol_frac[ i ];
			m_pressure_old[ i ]		  = m_pressure[ i ];		
			m_mean_velocity_old[ i ]  = m_mean_velocity[ i ];
			
		}
		
		this->set_bottom_pressure( m_HEEL_PRESSURE ); // Pressure at heel is set
		
		for(uint_type TIMESTEP = 0; TIMESTEP < FINAL_TIMESTEP; ++TIMESTEP){
			if(TIMESTEP){ 
				// Only enters loop for TIMESTEP > 0

				for( uint_type i = 0; i < number_of_nodes(); ++i){
					m_gas_vol_frac_old[ i ] = m_gas_vol_frac[ i ];	
					m_oil_vol_frac_old[ i ] = m_oil_vol_frac[ i ];					
					m_pressure_old[ i ] = m_pressure[ i ];	
					m_mean_velocity_old[ i ] = m_mean_velocity[ i ]; 				
				}
			}
			
			uint_type r = 0;
			real_type norma;
			do
			{
				
				this->compute_Jacobian();
				GMRES_Solve( *m_matrix, *m_variables, *m_source ); 
				this->update_variables();			
				
				//cout << setprecision(10);   				
				norma = itl::two_norm(*m_source);				
				//cout << "\n----norma residuo: " << norma << "-----\n";
				
				/*for( uint_type i = 0; i < number_of_nodes()-1; ++i ){
					cout << setprecision(10);				
					cout << "pressure[ "<< i <<" ] = " << this->m_pressure[ i ] <<"\t";
					if( i == 0)
						cout <<"\t";
					cout << "Gas_vol_frac[ "<< i <<" ] = " << this->m_gas_vol_frac[ i ] <<"\t";
					cout << "Oil_vol_frac[ "<< i <<" ] = " << this->m_oil_vol_frac[ i ] <<"\t";	
					if( i == number_of_nodes() - 1)
						cout <<"\t";				
					cout << "mean_velocity[ "<< i <<" ] = " << this->m_mean_velocity[ i ] << "\n";			

				}*/

				
				++r;
			}while(norma > this->NEWTON_CRIT && r < 100); 
			
			m_current_time += this->dt();
            double transient_norm_pressure      = 0.0;
            double transient_norm_gas_vol_frac  = 0.0;
            double transient_norm_oil_vol_frac  = 0.0;
            double transient_norm_velocity      = 0.0;
            for(unsigned i = 0; i < m_nnodes; ++i){
                transient_norm_pressure     += (m_pressure_old[i] - m_pressure[i])*(m_pressure_old[i] - m_pressure[i]);
                transient_norm_gas_vol_frac += (m_gas_vol_frac_old[i] - m_gas_vol_frac[i])*(m_gas_vol_frac_old[i] - m_gas_vol_frac[i]);	
                transient_norm_oil_vol_frac += (m_oil_vol_frac_old[i] - m_oil_vol_frac[i])*(m_oil_vol_frac_old[i] - m_oil_vol_frac[i]);
                transient_norm_velocity     += (m_mean_velocity_old[i] - m_mean_velocity[i])*(m_mean_velocity_old[i] - m_mean_velocity[i]);	
            }
            double transient_norm = sqrt( transient_norm_pressure    ) 
                + sqrt( transient_norm_gas_vol_frac)
                + sqrt( transient_norm_oil_vol_frac)
                + sqrt( transient_norm_velocity    );


			
            for( uint_type i = 0; i < number_of_nodes(); ++i){
                m_gas_vol_frac_old[ i ] = m_gas_vol_frac[ i ];	
                m_oil_vol_frac_old[ i ] = m_oil_vol_frac[ i ];					
                m_pressure_old[ i ] = m_pressure[ i ];	
                m_mean_velocity_old[ i ] = m_mean_velocity[ i ]; 	
            }

            m_total_production[OilPhase]    += m_oil_vol_frac[0]  *abs(m_oil_velocity[0])  *area()*dt();
            m_total_production[GasPhase]    += m_gas_vol_frac[0]  *abs(m_gas_velocity[0])  *area()*dt();
            m_total_production[WaterPhase]  += m_water_vol_frac[0]*abs(m_water_velocity[0])*area()*dt();
            /*cout << "------> well transient volume production: \n" 
                 << "\t\tPhaseWater: "<< m_total_production[WaterPhase]  << " m^3\n"
                 << "\t\tPhaseOil: "<< m_total_production[OilPhase]      << " m^3\n"
                 << "\t\tPhaseGas: "<< m_total_production[GasPhase]      << " m^3\n";*/
           

            if(transient_norm < 1.0e-3 && TIMESTEP > 0)
            {
                std::cout << "------> well transient TIME: " << m_current_time << " seconds\n";
                std::cout << "----WELL---- >>>TRANSIENT NORM = " << transient_norm << "\n";
                std::cout << "total inflow: \tOIL-> " << m_oil_vol_frac[0]  *abs(m_oil_velocity[0])  *area() << "\n" 
                    << "total inflow: \tGAS-> " << m_gas_vol_frac[0]  *abs(m_gas_velocity[0])  *area() << "\n" 
                    << "total inflow: \tWATER-> " << m_water_vol_frac[0]*abs(m_water_velocity[0])*area() << "\n";                 
                std::ofstream results_file;
                results_file.open( make_filename( "WellData\\results", STEPS, ".dat" ).c_str() );                 
                for( uint_type i = 0; i < number_of_nodes(); ++i )
                {
                    results_file << std::setprecision(10);				
                    results_file << this->m_pressure      [ i ]	<< "\t"
                        << this->m_gas_vol_frac  [ i ]	<< "\t"
                        << this->m_oil_vol_frac  [ i ]	<< "\t"
                        << this->m_water_vol_frac[ i ]	<< "\t"
                        << this->m_mean_velocity [ i ]	<< "\t"
                        << this->m_gas_velocity  [ i ]	<< "\t"
                        << this->m_oil_velocity  [ i ]	<< "\t"
                        << this->m_water_velocity[ i ]	<< "\n";												
                }
                results_file.close();
                ++STEPS;
                break;
            }
		}
        for(uint_type i = 1; i < this->m_pressure.size(); ++i){
            p_pressure[i-1] = this->m_pressure[i];
        }
        // TODO:
        // Calculate total volume of a phase 'p' produced:
        // something like
        // 
        // During transient time:
        // Total_vol_phase_p += m_vol_frac_p * velocity_p * Area * dt();
        // After reaching steady state:
        // TIME = Reservoir_timestep - total_time_well;
        // Total_vol_phase_p += m_vol_frac_p * velocity_p * Area * TIME;
        //real_type m_reservoir_timestep = 8640.0;
       // m_total_production[OilPhase]    += m_oil_vol_frac[0]  *abs(m_oil_velocity[0])  *area()*(m_reservoir_timestep-time);
        //m_total_production[GasPhase]    += m_gas_vol_frac[0]  *abs(m_gas_velocity[0])  *area()*(m_reservoir_timestep-time);
       // m_total_production[WaterPhase]  += m_water_vol_frac[0]*abs(m_water_velocity[0])*area()*(m_reservoir_timestep-time);
       /* cout << "-----------> well TOTAL volume production: \n" 
            << "\t\t---PhaseWater: "<< m_total_production[WaterPhase]  << " m^3\n"
            << "\t\t---PhaseOil: "<< m_total_production[OilPhase]      << " m^3\n"
            << "\t\t---PhaseGas: "<< m_total_production[GasPhase]      << " m^3\n";*/
	}




// Namespace =======================================================================================
} // namespace WellSimulator