#ifndef H_WellSimulator_NODECOORDINATES
#define H_WellSimulator_NODECOORDINATES

#include <stdexcept>

// Namespace =======================================================================================
namespace WellSimulator {
	enum Directions {X,Y,Z,N_DIRECTIONS=3};
// Vector ==========================================================================================
class NodeCoordinates
{
//--------------------------------------------------------------------------------- Type Definitions
public:
	typedef double					real_type;
	typedef real_type			    array_type[N_DIRECTIONS];
	typedef unsigned int			uint_type;
	typedef int						direction_type;
	

//------------------------------------------------------------------------- Constructor & Destructor
public:
	NodeCoordinates(){}		

	NodeCoordinates(const real_type& p_x, const real_type& p_y, const real_type& p_z)		
	{
		this->m_values[X] = p_x;
		this->m_values[Y] = p_y;
		this->m_values[Z] = p_z;
	}
	NodeCoordinates(const NodeCoordinates& p_other)		
	{
		for (direction_type d = X; d < N_DIRECTIONS; ++d){
			this->m_values[d] = p_other[d];
		}
	}

	~NodeCoordinates()
	{		
	}
//-------------------------------------------------------- Métodos de acesso e atribuição de valores
public:
	void operator =(const NodeCoordinates& p_other)		
	{
		for (direction_type d = X; d < N_DIRECTIONS; ++d){
			this->m_values[d] = p_other[d];
		}
	}


	real_type operator [](uint_type p_pos) const
	{
		if (p_pos >= N_DIRECTIONS){
			throw std::runtime_error("Exceeded position...");
		}
		return this->m_values[p_pos];
	}

	real_type& operator [](uint_type p_pos)
	{
		if (p_pos >= N_DIRECTIONS){
			throw std::runtime_error("Exceeded position...");
		}
		return this->m_values[p_pos];
	}
	
	real_type& getX(){
		return this->m_values[ X ];
	}
	real_type& getY(){
		return this->m_values[ Y ];
	}
	real_type& getZ(){
		return this->m_values[ Z ];
	}
//--------------------------------------------------------------------------------------------- Data
protected:
	array_type m_values;

};


// Namespace =======================================================================================
} // namespace WellSimulator

#endif // H_WellSimulator_NODECOORDINATES