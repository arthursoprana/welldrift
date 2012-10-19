#ifndef H_WellSimulator_GENERICWELL
#define H_WellSimulator_GENERICWELL

#include <AbstractWell.h>
#include <fstream>
#include <memory>
#include <Typedefs.h>
// Includes Boost
#include <BoostWrapper/SmartPointer.h>
// Namespace =======================================================================================
namespace WellSimulator {
using namespace boost_wrapper;  

// GenericWell ====================================================================================
class GenericWell : public AbstractWell
{
//--------------------------------------------------------------------------------- Type Definitions
public:
	

//------------------------------------------------------------------------- Constructor & Destructor
public:	
	GenericWell();
	GenericWell(
				const uint_type& p_nnodes,
				const real_type& p_radius
				);
	virtual ~GenericWell();
//---------------------------------------------------------- Métodos de interface com o reservatório
public:
	virtual void set_size(const uint_type& p_nnodes);
	virtual void initialize_flow( 
								 SharedPointer<vector_type> p_oil_flow_vector,
								 SharedPointer<vector_type> p_water_flow_vector,
								 SharedPointer<vector_type> p_gas_flow_vector
								 );
	
	virtual void set_radius(const real_type& p_radius);
	virtual void set_coordinates(const std::vector<coord_type>& p_coord_vector);
	virtual void read_coordinates(std::ifstream& p_infile);
	virtual coord_type* coordinates(uint_type p_index); 
	virtual real_type* pressure(uint_type p_index);
	virtual real_type radius() const;
	virtual uint_type number_of_nodes() const;	
	virtual void solve();
//--------------------------------------------------------------------------------------------- Data
protected:
	uint_type				m_nnodes;
	real_type				m_radius;
	std::vector<coord_type> m_coordinates; 
    SharedPointer<vector_type>	m_oil_flow;
	SharedPointer<vector_type>   m_water_flow;
	SharedPointer<vector_type>   m_gas_flow;
	vector_type				m_pressure;
};




// Namespace =======================================================================================
} // namespace WellSimulator

#endif // H_WellSimulator_GENERICWELL





