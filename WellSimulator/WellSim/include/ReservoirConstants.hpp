#ifndef HPP_WELLRES10_RESERVOIRCONSTANTS
#define HPP_WELLRES10_RESERVOIRCONSTANTS

// Includes ========================================================================================


// Namespace =======================================================================================
namespace WellSimulator {
typedef double float64;
// Constants =======================================================================================
class ReservoirConstants
{
//---------------------------------------------------------------------------------------- Constants
public:
    static const int number_of_unknowns   = 3;
    static const int number_of_phases     = 3;
    static const int number_of_components = 3;

    static inline float64 convert_mD_To_m2(){
        return 9.869233E-16;
    }

    static inline float64 convert_m2_To_mD(){
        return (1.0/9.869233E-16);
    }

    static inline float64 convert_ft_To_m(){
        return 0.3048;
    }

    static inline float64 convert_m_To_ft(){
        return (1.0/0.3048);
    }

    // Converts from [mD * m] to [cP.m3/(d.kPa)]
    // (86400 s * 1e3 Pa * 0.9869233e-15 m3 * m) / (1e-3 Pa * s * m3)
    static inline float64 conversion_factor() {
        return 8.527017312e-5; // (mD*m) / cP.m3/(d.kPa)
    }

    // Converts from [mD * m] to [cP.m3/(s.Pa)]
    // (s * Pa * 0.9869233e-15 m3 * m) / (1e-3 Pa * s * m3)
    static inline float64 well_conversion_factor() {
        return 9.869233e-13; // (mD*m) / cP.m3/(s.Pa)
    }

    static inline float64 convert_kPa_to_Pa() {
        return 1.0e3;
    }

    static inline float64 convert_Pa_to_kPa() {
        return 1.0e-3;
    }

    static inline float64 convert_cP_to_Pa_s() {
        return 1.0e-3;
    }

    static inline float64 convert_Dynes_per_cm_to_Pa_m() {
        return 1.0e-3;
    }

    static inline float64 convert_Pa_to_psi() {
        return 145.04e-6;
    }

    static inline float64 convert_day_to_s() {
        return 86400;
    }

    static inline float64 convert_s_to_day() {
        return (1.0/86400);
    }

    static inline float64 gravity_SI() {
        return 9.80665;
    }

    static inline float64 gravity(){
        return 9.80665e-3;
    }
};

typedef WellSimulator::ReservoirConstants WellConstants;

} // namespace WellSimulator

#endif // HPP_WELLRES10_RESERVOIRCONSTANTS
