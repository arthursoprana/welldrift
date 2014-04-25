#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <fstream>   
#include <string>
#include <cmath>
#include <vector>
#include <BoostWrapper/SmartPointer.h>
// Well Includes      
#include <DriftFluxWell.h>
using namespace boost_wrapper;
using namespace WellSimulator;



struct WellInitialData{
    float64 m_well_length;
    float64 m_well_diameter;
    float64 m_well_inclination;
    float64 m_heel_pressure;
    float64 m_toe_mixture_velocity;

    SharedPointer<vector_type> m_oil_inflow;
    SharedPointer<vector_type> m_water_inflow;
    SharedPointer<vector_type> m_gas_inflow;

    float64 m_initial_gas_vol_frac;
    float64 m_initial_oil_vol_frac;
    float64 m_initial_pressure;
    float64 m_initial_mixture_velocity;
    float64 m_final_time;
    float64 m_delta_t;
    float64 m_max_delta_t;
    float64 m_tolerance;
    int     m_number_of_nodes;
    // Fluid data
    float64 m_ref_temperature;
    float64 m_oil_density;
    float64 m_oil_viscosity;

    float64 m_water_density;
    float64 m_water_viscosity;

    float64 m_gas_ref_density;
    float64 m_gas_ref_pressure;
    float64 m_gas_sound_speed;
    float64 m_gas_viscosity;

    // Drift flux correlation parameters

};
using namespace std;
SharedPointer<DriftFluxWell> simulate(SharedPointer<WellInitialData> p_initial_data){
    //typedef WellSimulator::WellVector	vector_type;
    typedef std::vector<float64>     	vector_type;
    typedef NodeCoordinates				coord_type;   


    
    ////// CREATING WELL COORDINATE VECTOR //
    std::vector<coord_type> COORD_VECTOR( p_initial_data->m_number_of_nodes);

    float64 dx = p_initial_data->m_well_length/float64( p_initial_data->m_number_of_nodes - 1 );

    COORD_VECTOR[ 0 ][0] = 0;
    COORD_VECTOR[ 0 ][1] = 0;
    COORD_VECTOR[ 0 ][2] = 0;
    for( unsigned i = 1; i < p_initial_data->m_number_of_nodes; ++i ){
        COORD_VECTOR[ i ][0] = COORD_VECTOR[ i-1 ][0] + dx;
        COORD_VECTOR[ i ][1] = 0;
        COORD_VECTOR[ i ][2] = 0;			
    }


    ////// CREATED....

    SharedPointer<DriftFluxWell> wells( new DriftFluxWell(p_initial_data->m_number_of_nodes , 0.5*p_initial_data->m_well_diameter) );

    wells->set_inclination(p_initial_data->m_well_inclination);
    wells->set_constant_vol_frac( 
                                p_initial_data->m_initial_oil_vol_frac, 
                                p_initial_data->m_initial_gas_vol_frac, 
                                1.0 - (p_initial_data->m_initial_oil_vol_frac + p_initial_data->m_initial_gas_vol_frac) 
                                );	
    wells->set_dt( p_initial_data->m_delta_t );
    wells->set_max_delta_t( p_initial_data->m_max_delta_t );
    wells->set_final_time( p_initial_data->m_final_time );
    wells->set_final_timestep( p_initial_data->m_final_time/p_initial_data->m_delta_t );
    float64 PROFILE_PARAM_C_0 = 1.2;
    wells->set_C_0( PROFILE_PARAM_C_0 );
    // Define Drift-Flux Models
    SharedPointer<IDriftVelocityModel> GasLiquidDriftVelocityModel  (new ShiGasLiquidDriftVelocityModel(0.2,0.4));
    SharedPointer<IDriftVelocityModel> OilWaterDriftVelocityModel   (new ShiOilWaterDriftVelocityModel);

    SharedPointer<IProfileParameterModel> GasLiquidProfileParameterModel(new ShiGasLiquidProfileParameterModel(1.2,0.3,1.0));
    SharedPointer<IProfileParameterModel> OilWaterProfileParameterModel (new ShiOilWaterProfileParameterModel(1.2,0.4,0.7));

    real_type relative_standard_density = p_initial_data->m_oil_density/p_initial_data->m_water_density;
    SharedPointer<IInterfacialTensionModel> GasOilInterfacialTensionModel   (new BeggsGasOilInterfacialTensionModel  (p_initial_data->m_ref_temperature, relative_standard_density));
    SharedPointer<IInterfacialTensionModel> GasWaterInterfacialTensionModel (new BeggsGasWaterInterfacialTensionModel(p_initial_data->m_ref_temperature));     

    SharedPointer<IDensityModel> _GasDensityModel   (new WellCompressibleDensityModel(p_initial_data->m_gas_ref_density, p_initial_data->m_gas_ref_pressure, p_initial_data->m_gas_sound_speed) );
    SharedPointer<IDensityModel> _OilDensityModel   (new ConstantDensityModel  (p_initial_data->m_oil_density)  );
    SharedPointer<IDensityModel> _WaterDensityModel (new ConstantDensityModel  (p_initial_data->m_water_density) );

    SharedPointer<IViscosityModel> _GasViscosityModel   (new PowerViscosityModel (p_initial_data->m_gas_viscosity, 0.0) );
    SharedPointer<IViscosityModel> _OilViscosityModel   (new PowerViscosityModel (p_initial_data->m_oil_viscosity  , 0.0)  );
    SharedPointer<IViscosityModel> _WaterViscosityModel (new PowerViscosityModel (p_initial_data->m_water_viscosity  , 0.0) );

    // -----------------------

    wells->set_gas_liquid_drift_velocity_model   ( GasLiquidDriftVelocityModel   );
    wells->set_oil_water_drift_velocity_model    ( OilWaterDriftVelocityModel    );
    wells->set_gas_liquid_profile_parameter_model( GasLiquidProfileParameterModel);
    wells->set_oil_water_profile_parameter_model ( OilWaterProfileParameterModel );
    wells->set_gas_oil_interfacial_tension_model    ( GasOilInterfacialTensionModel );
    wells->set_gas_water_interfacial_tension_model  ( GasWaterInterfacialTensionModel );
    wells->set_gas_density_model  ( _GasDensityModel   );
    wells->set_oil_density_model  ( _OilDensityModel   );
    wells->set_water_density_model( _WaterDensityModel );
    wells->set_gas_viscosity_model  ( _GasViscosityModel   );
    wells->set_oil_viscosity_model  ( _OilViscosityModel   );
    wells->set_water_viscosity_model( _WaterViscosityModel );

    wells->set_delta( 1.220703125e-4, 1.220703125e-4, 1.220703125e-4, 1.220703125e-4 );
    wells->set_constant_pressure( p_initial_data->m_initial_pressure );
    wells->set_heel_pressure( p_initial_data->m_heel_pressure );
    wells->set_constant_velocity( p_initial_data->m_initial_mixture_velocity );


    wells->set_boundary_velocity( p_initial_data->m_toe_mixture_velocity );
    wells->set_with_gas( true );
    wells->set_mass_flux( false );
    wells->set_newton_criteria( p_initial_data->m_tolerance );

    inflow_vector_type inflow_gas(p_initial_data->m_number_of_nodes, MakeShared<ConstantInflow>(0.0));
    inflow_vector_type inflow_oil(p_initial_data->m_number_of_nodes, MakeShared<ConstantInflow>(0.0));
    inflow_vector_type inflow_water(p_initial_data->m_number_of_nodes, MakeShared<ConstantInflow>(0.0));

    for(int i = 0 ; i < p_initial_data->m_number_of_nodes; ++i){
        real_type Qoil = (*(p_initial_data->m_oil_inflow))[i];
        real_type Qwater = (*(p_initial_data->m_water_inflow))[i];
        real_type Qgas = (*(p_initial_data->m_gas_inflow))[i];           

        inflow_oil[i] = MakeShared<ConstantInflow>(Qoil);
        inflow_water[i] = MakeShared<ConstantInflow>(Qwater);
        inflow_gas[i] = MakeShared<ConstantInflow>(Qgas);             
    }

    wells->initialize_flow(inflow_oil,inflow_water,inflow_gas);

    wells->set_coordinates(COORD_VECTOR);
    wells->set_gravity( 0., 0., 9.8 );      

    wells->solve();	

    return wells;
}

SharedPointer<DriftFluxWell> simulate_provenzano(SharedPointer<WellInitialData> p_initial_data){
    //typedef WellSimulator::WellVector	vector_type;
    typedef std::vector<float64>     	vector_type;
    typedef NodeCoordinates				coord_type;   



    ////// CREATING WELL COORDINATE VECTOR //
    std::vector<coord_type> COORD_VECTOR( p_initial_data->m_number_of_nodes);

    float64 dx = p_initial_data->m_well_length/float64( p_initial_data->m_number_of_nodes - 1 );

    COORD_VECTOR[ 0 ][0] = 0;
    COORD_VECTOR[ 0 ][1] = 0;
    COORD_VECTOR[ 0 ][2] = 0;
    for( unsigned i = 1; i < p_initial_data->m_number_of_nodes; ++i ){
        COORD_VECTOR[ i ][0] = COORD_VECTOR[ i-1 ][0] + dx;
        COORD_VECTOR[ i ][1] = 0;
        COORD_VECTOR[ i ][2] = 0;			
    }


    ////// CREATED....

    SharedPointer<DriftFluxWell> wells( new DriftFluxWell(p_initial_data->m_number_of_nodes , 0.5*p_initial_data->m_well_diameter) );

    wells->set_has_inclination_correction(false);
    wells->set_inclination(p_initial_data->m_well_inclination);
    wells->set_constant_vol_frac( 
        p_initial_data->m_initial_oil_vol_frac, 
        p_initial_data->m_initial_gas_vol_frac, 
        1.0 - (p_initial_data->m_initial_oil_vol_frac + p_initial_data->m_initial_gas_vol_frac) 
        );	
    wells->set_dt( p_initial_data->m_delta_t );
    wells->set_max_delta_t( p_initial_data->m_max_delta_t );
    wells->set_final_time( p_initial_data->m_final_time );
    wells->set_final_timestep( p_initial_data->m_final_time/p_initial_data->m_delta_t );
    float64 PROFILE_PARAM_C_0 = 1.2;
    wells->set_C_0( PROFILE_PARAM_C_0 );
    // Define Drift-Flux Models
    SharedPointer<IDriftVelocityModel> GasLiquidDriftVelocityModel  (new ConstantDriftVelocityModel(-0.5));
    SharedPointer<IDriftVelocityModel> OilWaterDriftVelocityModel   (new ConstantDriftVelocityModel(0.0));

    SharedPointer<IProfileParameterModel> GasLiquidProfileParameterModel(new ConstantProfileParameterModel(1.2));
    SharedPointer<IProfileParameterModel> OilWaterProfileParameterModel (new ConstantProfileParameterModel(1.0));

    real_type relative_standard_density = p_initial_data->m_oil_density/p_initial_data->m_water_density;
    SharedPointer<IInterfacialTensionModel> GasOilInterfacialTensionModel   (new BeggsGasOilInterfacialTensionModel  (p_initial_data->m_ref_temperature, relative_standard_density));
    SharedPointer<IInterfacialTensionModel> GasWaterInterfacialTensionModel (new BeggsGasWaterInterfacialTensionModel(p_initial_data->m_ref_temperature));     

    SharedPointer<IDensityModel> _GasDensityModel   (new WellCompressibleDensityModel(  p_initial_data->m_gas_ref_density, 
                                                                                        p_initial_data->m_gas_ref_pressure, 
                                                                                        p_initial_data->m_gas_sound_speed) 
                                                    );

    SharedPointer<IDensityModel> _OilDensityModel   (new WellCompressibleDensityModel(  1000.0, 
                                                                                        1.0e5, 
                                                                                        1000.0) 
                                                    );

    SharedPointer<IDensityModel> _WaterDensityModel (new ConstantDensityModel  (p_initial_data->m_water_density) );

    SharedPointer<IViscosityModel> _GasViscosityModel   (new PowerViscosityModel (p_initial_data->m_gas_viscosity, 0.0) );
    SharedPointer<IViscosityModel> _OilViscosityModel   (new PowerViscosityModel (p_initial_data->m_oil_viscosity  , 0.0)  );
    SharedPointer<IViscosityModel> _WaterViscosityModel (new PowerViscosityModel (p_initial_data->m_water_viscosity  , 0.0) );

    // -----------------------

    wells->set_gas_liquid_drift_velocity_model   ( GasLiquidDriftVelocityModel   );
    wells->set_oil_water_drift_velocity_model    ( OilWaterDriftVelocityModel    );
    wells->set_gas_liquid_profile_parameter_model( GasLiquidProfileParameterModel);
    wells->set_oil_water_profile_parameter_model ( OilWaterProfileParameterModel );
    wells->set_gas_oil_interfacial_tension_model    ( GasOilInterfacialTensionModel );
    wells->set_gas_water_interfacial_tension_model  ( GasWaterInterfacialTensionModel );
    wells->set_gas_density_model  ( _GasDensityModel   );
    wells->set_oil_density_model  ( _OilDensityModel   );
    wells->set_water_density_model( _WaterDensityModel );
    wells->set_gas_viscosity_model  ( _GasViscosityModel   );
    wells->set_oil_viscosity_model  ( _OilViscosityModel   );
    wells->set_water_viscosity_model( _WaterViscosityModel );

    wells->set_delta( 1.220703125e-4, 1.220703125e-4, 1.220703125e-4, 1.220703125e-4 );
    wells->set_constant_pressure( p_initial_data->m_initial_pressure );
    wells->set_heel_pressure( p_initial_data->m_heel_pressure );
    wells->set_constant_velocity( p_initial_data->m_initial_mixture_velocity );


    wells->set_boundary_velocity( p_initial_data->m_toe_mixture_velocity );
    wells->set_with_gas( true );
    wells->set_mass_flux( true ); // mass inflow flux
    wells->set_newton_criteria( p_initial_data->m_tolerance );

    real_type n_nodes = p_initial_data->m_number_of_nodes;
    real_type length = p_initial_data->m_well_length;
    real_type delta_x = length/n_nodes;

    inflow_vector_type inflow_gas(n_nodes, MakeShared<ConstantInflow>(0.0));
    inflow_vector_type inflow_oil(n_nodes, MakeShared<ConstantInflow>(0.0));
    inflow_vector_type inflow_water(n_nodes, MakeShared<ConstantInflow>(0.0));  
    
    for(int i = 0 ; i < n_nodes; ++i){
        
        real_type pos = i*delta_x + 0.5*delta_x;

        if( pos > 450.0 && pos < 550.0){
            real_type influx_length = 550.0 - 450.0;
            real_type Qoil_per_meter   = 0.0  / influx_length;
            real_type Qwater_per_meter = 0.0  / influx_length;
            real_type Qgas_per_meter   = 0.8 / influx_length;           

            real_type Qoil   = Qoil_per_meter   * delta_x;
            real_type Qwater = Qwater_per_meter * delta_x;
            real_type Qgas   = Qgas_per_meter   * delta_x;           

            inflow_oil[i] = MakeShared<SlugInflow>(Qoil, 0.0, 70.0, 10.0, 20.0);
            inflow_water[i] = MakeShared<ConstantInflow>(Qwater);
            inflow_gas[i] = MakeShared<SlugInflow>(Qgas, 0.0, 70.0, 10.0, 20.0);
        }              
    }

    uint_type toe_idx = n_nodes-1;
    // CASE 1
    //inflow_gas[ toe_idx ] = MakeShared<ConstantInflow>(0.02);
    //inflow_oil[ toe_idx ] = MakeShared<ConstantInflow>(3.0);

    // CASE 2
    inflow_gas[ toe_idx ] = MakeShared<ProvenzanoCase2GasInflow>();
    inflow_oil[ toe_idx ] = MakeShared<ProvenzanoCase2OilInflow>(12.0);

    wells->initialize_flow(inflow_oil,inflow_water,inflow_gas);

    wells->set_coordinates(COORD_VECTOR);
    wells->set_gravity( 0., 0., 9.8 );      

    wells->solve();	

    return wells;
}
#endif // SIMULATOR_H