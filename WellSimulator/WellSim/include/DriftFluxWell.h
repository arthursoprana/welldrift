#ifndef H_WellSimulator_DRIFTFLUXWELL
#define H_WellSimulator_DRIFTFLUXWELL
#include <Models.hpp>

#include <SharedPointer.h>

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

//#include <WellSolver.h>
#include <GenericWell.h>
#include <string>




// Namespace =======================================================================================
namespace WellSimulator {

    typedef mtl::dense1D< double >  svector_type;
    typedef mtl::matrix< double, mtl::rectangle<>, mtl::dense<>, mtl::row_major>::type 					   dmatrix_type; // MATRIZ DENSA
    typedef mtl::matrix< double, mtl::rectangle<>, mtl::array< mtl::compressed<> >, mtl::row_major >::type smatrix_type; // MATRIZ ESPARSA
    typedef SharedPointer<svector_type> vector_ptr;
    typedef SharedPointer<smatrix_type> matrix_ptr;

	typedef char									string_type;
	typedef double									real_type;
	typedef unsigned int							uint_type;
	//typedef WellSimulator::WellVector				vector_type;
    typedef std::vector<double>				vector_type;
	typedef NodeCoordinates							coord_type;
	typedef std::vector< std::vector<uint_type> >	id_type;
	enum	v_variables{P, alpha_g, alpha_o, v,total_var = 4};



    

	// DriftFluxWell ===================================================================================
	class DriftFluxWell : public GenericWell
	{	
		//------------------------------------------------------------------------- Constructor & Destructor
	public:
		DriftFluxWell();
		DriftFluxWell(
			const uint_type& p_nnodes,
			const real_type& p_radius
			);
		virtual ~DriftFluxWell();

		//--------------------------------------------------------------------- DriftFluxWell main functions
	public:
        void Initialize(
            uint_type   p_number_of_completions, 
            uint_type   p_direction, 
            real_type   p_well_radius, 
            real_type   p_BHPressure, 
            vector_type p_well_length);

		real_type Volume( real_type dS );
		real_type dt();
		void set_dt( real_type p_dt );
		real_type area();
		real_type ksi( real_type p_velocity );

		real_type mean_density(
			const real_type& p_oil_vol_frac,
			const real_type& p_water_vol_frac,
			const real_type& p_gas_vol_frac,
			const real_type& p_pressure
			);
		void set_mean_velocity();
		void set_constant_oil_vol_frac( real_type p_value );
		void set_constant_vol_frac( 
			real_type p_oil_vol_frac, 
			real_type p_gas_vol_frac, 
			real_type p_water_vol_frac 
			);
		void set_constant_pressure( real_type p_pressure );
		void set_constant_velocity( real_type p_velocity );

		void solve();
        void solve(vector_type& p_pressure);

		void GMRES_Solve( smatrix_type &A, svector_type &x, svector_type &b );
		void compute_Jacobian();
		void update_variables();
        void update_variables_for_new_timestep();
		
		real_type liquid_density(
			real_type p_oil_vol_frac,
			real_type p_water_vol_frac,
			real_type p_pressure
			);
		real_type gas_density	( real_type p_pressure );
		real_type oil_density	( real_type p_pressure );
		real_type water_density	( real_type p_pressure );

		real_type mean_velocity( uint_type p_index );
		real_type C_0();
		void set_C_0(real_type p_C_0);
		real_type v_drift_flux( real_type p_gas_vol_frac, real_type p_pressure, real_type p_liquid_density );	
		real_type friction_factor( real_type p_reynolds );

		real_type mod_v_drift_flux(											
			real_type p_mean_velocity, 
			real_type p_gas_vol_frac, 
			real_type p_oil_vol_frac,
			real_type p_water_vol_frac,
			real_type p_pressure  
			); // Gas and Liquid Drift Velocity
		real_type mod_v_drift_flux_ow(											
			real_type p_mean_velocity, 
			real_type p_gas_vol_frac, 
			real_type p_oil_vol_frac,
			real_type p_water_vol_frac,
			real_type p_pressure  
			); // Oil and Water Drift Velocity

		void set_gravity( real_type p_valueX, real_type p_valueY, real_type p_valueZ );
		real_type gravity();
		void set_delta( real_type p_delta_P, real_type p_delta_alpha_g, real_type p_delta_alpha_o, real_type p_delta_v );
		void set_bottom_pressure( real_type p_pressure );		
		void set_boundary_velocity( real_type p_velocity );
		

		real_type dot( vector_type& p_vec1, vector_type& p_vec2 );

		real_type gas_viscosity( real_type p_pressure );
		real_type oil_viscosity( real_type p_pressure );
		real_type water_viscosity( real_type p_pressure );

		uint_type id( uint_type p_node , uint_type p_variable );
		real_type segment_length( coord_type p_coord_i, coord_type p_coord_j );

        real_type GasMassFlowRate() {return gas_density(m_pressure[0])*m_gas_vol_frac[0]*abs(m_gas_velocity[0])*area();}
        real_type OilMassFlowRate() {return oil_density(m_pressure[0])*m_oil_vol_frac[0]*abs(m_oil_velocity[0])*area();}				
      
        real_type DriftCompletionGasFlowRate  () {return gas_density  (m_pressure[0])*m_gas_vol_frac  [0]*abs(m_gas_velocity[0])*area();}
        real_type DriftCompletionOilFlowRate  () {return oil_density  (m_pressure[0])*m_oil_vol_frac  [0]*abs(m_oil_velocity[0])*area();}				
        real_type DriftCompletionWaterFlowRate() {return water_density(m_pressure[0])*m_water_vol_frac[0]*abs(m_water_velocity[0])*area();}

		real_type R_m( 
			real_type	p_pressureW,
			real_type	p_pressureP,
			real_type	p_pressureE,
			real_type	p_gas_vol_fracW,
			real_type	p_gas_vol_fracP,
			real_type	p_gas_vol_fracE,
			real_type	p_oil_vol_fracW,
			real_type	p_oil_vol_fracP,
			real_type	p_oil_vol_fracE,
			real_type	p_velocityW,
			real_type	p_velocityP,
			uint_type	p_node,
			string_type	position  
			);
		real_type R_g( 
			real_type	p_pressureW,
			real_type	p_pressureP,
			real_type	p_pressureE,
			real_type	p_gas_vol_fracW,
			real_type	p_gas_vol_fracP,
			real_type	p_gas_vol_fracE,
			real_type	p_oil_vol_fracW,
			real_type	p_oil_vol_fracP,
			real_type	p_oil_vol_fracE,
			real_type	p_velocityW,
			real_type	p_velocityP,
			uint_type	p_node,
			string_type	position
			);
		real_type R_o( 
			real_type	p_pressureW,
			real_type	p_pressureP,
			real_type	p_pressureE,
			real_type	p_gas_vol_fracW,
			real_type	p_gas_vol_fracP,
			real_type	p_gas_vol_fracE,
			real_type	p_oil_vol_fracW,
			real_type	p_oil_vol_fracP,
			real_type	p_oil_vol_fracE,
			real_type	p_velocityW,
			real_type	p_velocityP,
			uint_type	p_node,
			string_type	position
			);

		real_type R_v(
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
			string_type    position
			);




		void set_with_gas (bool p_choice = false);
        void set_mass_flux(bool p_choice = false){
            m_mass_flux = p_choice;
        }
		void set_newton_criteria(real_type p_tolerance);
		void set_final_timestep(uint_type p_final_timestep );
		void set_heel_pressure(real_type p_pressure){
			m_HEEL_PRESSURE = p_pressure;
		}

        void set_gas_liquid_drift_velocity_model(SharedPointer<IDriftVelocityModel> p_gas_liquid_drift_velocity_model){
            m_gas_liquid_drift_velocity_model = p_gas_liquid_drift_velocity_model;
        }
        void set_oil_water_drift_velocity_model(SharedPointer<IDriftVelocityModel> p_oil_water_drift_velocity_model){
            m_oil_water_drift_velocity_model = p_oil_water_drift_velocity_model;
        }
        void set_gas_liquid_profile_parameter_model(SharedPointer<IProfileParameterModel> p_gas_liquid_profile_parameter_model){
            m_gas_liquid_profile_parameter_model = p_gas_liquid_profile_parameter_model;
        }
        void set_oil_water_profile_parameter_model(SharedPointer<IProfileParameterModel> p_oil_water_profile_parameter_model){
            m_oil_water_profile_parameter_model = p_oil_water_profile_parameter_model;
        }
        void set_gas_oil_interfacial_tension_model(SharedPointer<IInterfacialTensionModel> p_gas_oil_interfacial_tension_model){
            m_gas_oil_interfacial_tension_model = p_gas_oil_interfacial_tension_model;
        }
        void set_gas_water_interfacial_tension_model(SharedPointer<IInterfacialTensionModel> p_gas_water_interfacial_tension_model){
            m_gas_water_interfacial_tension_model = p_gas_water_interfacial_tension_model;
        }
        void set_gas_density_model(SharedPointer<IDensityModel> p_gas_density_model){
            m_gas_density_model = p_gas_density_model;
        }
        void set_oil_density_model(SharedPointer<IDensityModel> p_oil_density_model){
            m_oil_density_model = p_oil_density_model;
        }
        void set_water_density_model(SharedPointer<IDensityModel> p_water_density_model){
            m_water_density_model = p_water_density_model;
        }
        void set_gas_viscosity_model(SharedPointer<IViscosityModel> p_gas_viscosity_model){
            m_gas_viscosity_model = p_gas_viscosity_model;
        }
        void set_oil_viscosity_model(SharedPointer<IViscosityModel> p_oil_viscosity_model){
            m_oil_viscosity_model = p_oil_viscosity_model;
        }
        void set_water_viscosity_model(SharedPointer<IViscosityModel> p_water_viscosity_model){
            m_water_viscosity_model = p_water_viscosity_model;
        }

        real_type get_gas_volume_fraction(uint_type p_index){ return m_gas_vol_frac[p_index]; }
        real_type get_oil_volume_fraction(uint_type p_index){ return m_oil_vol_frac[p_index]; }
        real_type get_water_volume_fraction(uint_type p_index){ return m_water_vol_frac[p_index]; }
        real_type get_mixture_velocity(uint_type p_index){ return m_mean_velocity[p_index]; }
        real_type get_gas_volumetric_flux(uint_type p_index){ return m_gas_vol_frac[p_index]*m_gas_velocity[p_index]*area(); }
        real_type get_oil_volumetric_flux(uint_type p_index){ return m_oil_vol_frac[p_index]*m_oil_velocity[p_index]*area(); }
        real_type get_water_volumetric_flux(uint_type p_index){ return m_water_vol_frac[p_index]*m_water_velocity[p_index]*area(); }

        void set_has_inclination_correction(bool p_has_inclination_correction){
            m_has_inclination_correction = p_has_inclination_correction;
        }

        void set_inclination(real_type p_inclination){
            m_well_inclination = p_inclination;
        }

        real_type get_inclination(){
            return m_well_inclination;
        }

        void set_final_time(real_type p_final_time){
            m_final_time = p_final_time;
        }

        void set_max_delta_t(real_type p_max_delta_t){
            m_max_delta_t = p_max_delta_t;
        }                     

        real_type calculate_new_delta_t_size_converged_solution(real_type delta_t_old);
        real_type calculate_new_delta_t_size_diverged_solution(real_type delta_t_old);
        void restore_initial_guess();         

		//--------------------------------------------------------------------------------------------- Data
	protected:
        // DRIFT MODELS
        SharedPointer<IDriftVelocityModel>       m_gas_liquid_drift_velocity_model;
        SharedPointer<IDriftVelocityModel>       m_oil_water_drift_velocity_model;
        SharedPointer<IProfileParameterModel>    m_gas_liquid_profile_parameter_model;
        SharedPointer<IProfileParameterModel>    m_oil_water_profile_parameter_model;
        SharedPointer<IInterfacialTensionModel>  m_gas_oil_interfacial_tension_model;
        SharedPointer<IInterfacialTensionModel>  m_gas_water_interfacial_tension_model;
        SharedPointer<IDensityModel>    m_gas_density_model;
        SharedPointer<IDensityModel>    m_oil_density_model;
        SharedPointer<IDensityModel>    m_water_density_model;
        SharedPointer<IViscosityModel>    m_gas_viscosity_model;
        SharedPointer<IViscosityModel>    m_oil_viscosity_model;
        SharedPointer<IViscosityModel>    m_water_viscosity_model;
        
        
		vector_type m_oil_velocity;
		vector_type m_water_velocity;
		vector_type m_gas_velocity;
		vector_type m_mean_velocity;
		vector_type m_oil_vol_frac;
		vector_type m_water_vol_frac;
		vector_type m_gas_vol_frac;
		
		vector_type m_gas_vol_frac_old;
		vector_type m_oil_vol_frac_old;
		vector_type m_water_vol_frac_old;
		vector_type m_pressure_old;
		vector_type m_mean_velocity_old;

		bool		m_with_gas;
        bool        m_mass_flux;
        bool        m_convergence_status;
        bool        m_has_inclination_correction;
		real_type	m_profile_parameter_C_0;
		real_type   m_HEEL_PRESSURE;
		real_type   m_oil_API;
		vector_type m_gravity; // gravity vector
		vector_type m_delta;
        vector_type m_total_production; // [ m³ ]
		real_type   m_dt;
        real_type m_well_inclination;
		uint_type   m_FINAL_TIMESTEP;
        real_type   m_final_time;
        real_type   m_max_delta_t;
        real_type   m_current_time;
		real_type NEWTON_CRIT;

		vector_ptr m_variables;
		vector_ptr m_source;
		matrix_ptr m_matrix;

		id_type		 m_id;


	}; // class DriftFluxWell


// Namespace =======================================================================================
} // namespace WellSimulator

#endif // H_WellSimulator_DRIFTFLUXWELL