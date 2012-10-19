#ifndef H_WellSimulator_WELLVECTOR
#define H_WellSimulator_WELLVECTOR

#include <vector>
#include <stdexcept>
#include <cmath>

// Namespace =======================================================================================
namespace WellSimulator {

// WellVector ==========================================================================================
class WellVector
{
//--------------------------------------------------------------------------------- Type Definitions
public:
	typedef double					real_type;
	typedef std::vector<real_type>  array_type;
	typedef unsigned int			uint_type;


//------------------------------------------------------------------------- Constructor & Destructor
public:
	WellVector(){}
	WellVector(uint_type p_size)
		: m_values(p_size)
	{
	}

	WellVector(uint_type p_size, real_type p_value)
		: m_values(p_size, p_value)
	{
	}

	~WellVector()
	{
		this->m_values.clear();
	}
//-------------------------------------------------------- Métodos de acesso e atribuição de valores
public:
	void operator =(WellVector& p_vector)
	{
		this->m_values = p_vector.m_values;
	}

	real_type operator [](uint_type p_pos) const
	{
		if (p_pos >= this->m_values.size()){
			throw std::runtime_error("Exceeded position...");
		}
		return this->m_values[p_pos];
	}

	real_type& operator [](uint_type p_pos)
	{
		if (p_pos >= this->m_values.size()){
			throw std::runtime_error("Exceeded position...");
		}
		return this->m_values[p_pos];
	}

	inline void resize(uint_type p_newsize)
	{
		this->m_values.resize(p_newsize);
	}

	inline void assign(uint_type p_size, real_type p_val)
	{
		this->m_values.assign( p_size, p_val);
	}

	inline uint_type size()
	{
		return this->m_values.size();
	}
	inline real_type norm()
	{
		real_type Norm = 0;
		for( uint_type i = 0; i < this->size(); ++i )
			Norm += m_values[ i ]*m_values[ i ];
		return sqrt(Norm);
	}
//--------------------------------------------------------------------------------------------- Data
protected:
	array_type m_values;

};


// Namespace =======================================================================================
} // namespace WellSimulator

#endif // H_WellSimulator_WELLVECTOR
