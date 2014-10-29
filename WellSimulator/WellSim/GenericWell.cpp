#include <GenericWell.h>

#include <iostream>
#include <fstream>
#include <iomanip>

// Namespace =======================================================================================
namespace WellSimulator {

	typedef double						real_type;
	typedef unsigned int				uint_type;
	//typedef WellSimulator::WellVector   vector_type;
    typedef std::vector<double>				vector_type;
	typedef NodeCoordinates				coord_type;

	GenericWell::GenericWell()
	{
	}

	GenericWell::GenericWell(
							 const uint_type& p_nnodes,
							 const real_type& p_radius
							 )
							 : m_coordinates( p_nnodes )
	{	
		this->m_nnodes = p_nnodes;
		this->m_radius = p_radius;
	}

	GenericWell::~GenericWell(){}

	void GenericWell::set_size(const uint_type &p_nnodes){
		this->m_nnodes = p_nnodes;
	}

	void GenericWell::initialize_flow( 
									 inflow_vector_type p_oil_flow_vector,
									 inflow_vector_type p_water_flow_vector,
									 inflow_vector_type p_gas_flow_vector
									 )
	{
       /* for(uint_type i = 1; i < this->m_oil_flow->size(); ++i){
            (*this->m_oil_flow)[i]	= p_oil_flow_vector[i-1];
            (*this->m_water_flow)[i]= p_water_flow_vector[i-1];
            (*this->m_gas_flow)[i]	= p_gas_flow_vector[i-1];
        }*/      
        this->m_oil_flow	= p_oil_flow_vector;
		this->m_water_flow	= p_water_flow_vector;
		this->m_gas_flow	= p_gas_flow_vector;
	}
	void GenericWell::set_radius(const real_type& p_radius)	{
		this->m_radius = p_radius;
	}
	void GenericWell::set_coordinates(const std::vector<coord_type>& p_coord_vector)
	{
		std::ofstream coordinates_file("..\\WellData\\coordinates.txt");					
		for( uint_type i = 0; i < this->number_of_nodes(); ++i )
		{			
			this->m_coordinates[ i ] = p_coord_vector[ i ];

			coordinates_file << std::setprecision(10);				
			coordinates_file << m_coordinates[ i ][ 0 ] 
					 << "\t" << m_coordinates[ i ][ 1 ] 
					 << "\t" << m_coordinates[ i ][ 2 ] << "\n";		
		}
		coordinates_file.close();
	}
	void GenericWell::read_coordinates(std::ifstream& p_infile)
	{
		p_infile.ignore( 50 ,'\n' );
		p_infile.ignore( 50 ,'\n' );
		real_type tempX, tempY, tempZ;
		for( uint_type i = 0; i < this->number_of_nodes(); ++i )
		{
			p_infile.ignore( 12 );
			p_infile >> tempX >> tempY >> tempZ;
			NodeCoordinates tempCoord( tempX, tempY, tempZ );
			this->m_coordinates[ i ] = tempCoord;
		}
	}
	coord_type* GenericWell::coordinates(uint_type p_index)	{
		return &(this->m_coordinates[p_index]);	
	}
	real_type* GenericWell::pressure(uint_type p_index)	{
		return &(this->m_pressure[p_index]);		
	}
	real_type GenericWell::radius() const{
		return this->m_radius;
	}
	uint_type GenericWell::number_of_nodes() const{
		return this->m_nnodes;
	}
	void GenericWell::solve()
	{
	}


// Namespace =======================================================================================
} // namespace WellSimulator
