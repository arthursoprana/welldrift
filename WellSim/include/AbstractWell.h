#ifndef H_WellSimulator_ABSTRACTWELL
#define H_WellSimulator_ABSTRACTWELL

#include <WellVector.h>
#include <NodeCoordinates.h>
#include <Typedefs.h>
#include <BoostWrapper/SmartPointer.h>

// Namespace =======================================================================================
namespace WellSimulator {

    class IPhaseInflowExpression {
    public:
        virtual void calculate_value_at_time(real_type p_time) = 0;

        real_type get_current_value(){
            return m_value;
        }

    protected:
        real_type m_value;
    };

    class ConstantInflow 
        : public IPhaseInflowExpression
    {
    public:
        ConstantInflow(real_type p_value)  
        {   
            m_value = p_value;
        }

        void calculate_value_at_time(real_type p_time){     
        }
    };

    typedef std::vector< SharedPointer<IPhaseInflowExpression> > inflow_vector_type;

// AbstractWell ===================================================================================
class AbstractWell
{
//--------------------------------------------------------------------------------- Type Definitions
public:
	//typedef double						real_type;
	//typedef unsigned int				uint_type;
	////typedef WellSimulator::WellVector	vector_type;
 //   typedef std::vector<double>				vector_type;
	//typedef NodeCoordinates				coord_type;

//------------------------------------------------------------------------- Constructor & Destructor
public:
	AbstractWell(){}
	virtual ~AbstractWell(){}
//---------------------------------------------------------- Métodos de interface com o reservatório
public:
	virtual void set_size(const uint_type& p_nnodes) = 0;
	virtual void initialize_flow(
								 inflow_vector_type p_oil_flow_vector,
								 inflow_vector_type p_water_flow_vector,
								 inflow_vector_type p_gas_flow_vector
								 ) = 0;
	virtual void set_radius(const real_type& p_radius) = 0;
	virtual void read_coordinates( std::ifstream& p_infile) = 0;
	virtual coord_type* coordinates(uint_type p_index) = 0;
	virtual real_type* pressure(uint_type p_index) = 0;
	virtual real_type radius() const = 0;
	virtual uint_type number_of_nodes() const = 0;
	virtual void solve() = 0;

};

// Namespace =======================================================================================
} // namespace WellSimulator

#endif // H_WELLRES10_ABSTRACTWELL

/*
{
.
.    Criação dos poços, armazenamento em wells (vetor de poços)
.
	for( int i = 0; i < number_f_wells; ++i ){
		calculaFluxoNoPoço( i, flow_vector );
		wells[i].initialize_flow( flow_vector );
	}
	while( residuo maior que tolerancia ){
		zeraSistemaDeEquaçõesDoReservatorio();
		for( int i = 0; i < number_f_wells; ++i ){
			wells[i].solve();
			double radius = wells[i].radius();
			NodeCoordinates* coord;
			double* wellPressure;
			for( int j = 0; j < wells[i].number_of_nodes(); ++j ){
				coord = well[i].coordinates( j );
				wellPressure = well[i].pressure( j );

				Localiza o elemento que contém este nó do poço;
				Calcula vazão utilizando a pressão do reservatório e wellPressure;
				Adiciona a vazao no sistema de equações do reservatório;
			}
		}
		terminaAContrucaoDoSistemaDeEquacoesDoRes();
		resolveOReservatorio();
	}
.
.
.}

*/
