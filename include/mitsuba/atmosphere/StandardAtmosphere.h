#ifndef _STANDARD_ATMOSPHERE_H_
#define _STANDARD_ATMOSPHERE_H_

#include <vector>
#include <mitsuba/atmosphere/Utils.h>

namespace StandardAtmosphere {

	const static std::vector<std::vector<float>> TabulatedValues1 =
	{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86 },
	{ 288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65, 236.15, 229.65, 223.15, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 217.65, 218.65, 219.65, 220.65, 221.65, 222.65, 223.65, 224.65, 225.65, 226.65, 227.65, 228.65, 231.45, 234.25, 237.05, 239.85, 242.65, 245.45, 248.25, 251.05, 253.85, 256.65, 259.45, 262.25, 265.05, 267.85, 270.65, 270.65, 270.65, 270.65, 270.65, 267.85, 265.05, 262.25, 259.45, 256.65, 253.85, 251.05, 248.25, 245.45, 242.65, 239.85, 237.05, 234.25, 231.45, 228.65, 225.85, 223.05, 220.25, 217.45, 214.65, 212.65, 210.65, 208.65, 206.65, 204.65, 202.65, 200.65, 198.65, 196.65, 194.65, 192.65, 190.65, 188.65, 186.946, 186.946 },
	{ 101325, 89874.6, 79495.2, 70108.5, 61640.2, 54019.9, 47181, 41060.7, 35599.8, 30742.5, 26436.3, 22632.1, 19330.4, 16510.4, 14101.8, 12044.6, 10287.5, 8786.68, 7504.84, 6410.01, 5474.89, 4677.89, 3999.79, 3422.43, 2930.49, 2511.02, 2153.09, 1847.46, 1586.29, 1362.96, 1171.87, 1008.23, 868.019, 748.228, 646.122, 558.924, 484.317, 420.367, 365.455, 318.220, 277.522, 242.395, 212.030, 185.738, 162.937, 143.135, 125.910, 110.906, 97.7545, 86.1623, 75.9448, 66.9389, 58.9622, 51.8668, 45.5632, 39.9700, 35.0137, 30.6274, 26.7509, 23.3296, 20.3143, 17.6606, 15.3287, 13.2826, 11.4900, 9.92203, 8.55275, 7.35895, 6.31992, 5.41717, 4.63422, 3.95642, 3.37176, 2.86917, 2.43773, 2.06792, 1.75140, 1.48092, 1.25012, 1.05351, 0.88628, 0.74428, 0.623905, 0.522037, 0.435981, 0.36342, 0.302723 },
	{ 1.225, 1.11164, 1.00649, 0.909122, 0.819129, 0.736116, 0.659697, 0.589501, 0.525168, 0.466348, 0.412707, 0.363918, 0.310828, 0.265483, 0.226753, 0.193674, 0.16542, 0.141288, 0.120676, 0.103071, 0.0880349, 0.0748737, 0.0637273, 0.0542803, 0.0462674, 0.0394658, 0.0336882, 0.0287769, 0.0245988, 0.021042, 0.0180119, 0.0154288, 0.013225, 0.011262, 0.00960889, 0.00821392, 0.00703441, 0.00603513, 0.00518691, 0.00446557, 0.00385101, 0.00332648, 0.00287802, 0.00249393, 0.00216443, 0.00188129, 0.0016376, 0.00142753, 0.00125825, 0.00110904, 0.000977525, 0.000861606, 0.000766867, 0.00068171, 0.000605252, 0.000536684, 0.000475263, 0.000420311, 0.000371207, 0.000327382, 0.000288321, 0.00025355, 0.00022264, 0.0001952, 0.000170875, 0.000149342, 0.000130308, 0.00011351, 9.87069e-05, 0.000085683, 0.000074243, 0.000064211, 0.000055237, 4.74496e-05, 0.000040701, 3.48607e-05, 2.98135e-05, 2.54579e-05, 2.17046e-05, 1.84751e-05, 1.57005e-05, 1.33205e-05, 0.000011282, 9.53899e-06, 8.05098e-06, 6.77222e-06, 5.64114e-06 } };

	/*
	* [1] Dutsch '73. The Ozone distribution in the atmosphere.
	*
	*	It stores the monthly mean values and standard deviation of 
	*	total ozone from 45 years in Arosa, Switzerland(47ºN) in
	*	Dobson(10e-3 cm normal temperature and pressure).
	*	row 1 = month, row 2 = value, row 3 = std deviation
	*/

	const static std::vector<std::vector<float>> ozoneDobson =
	{ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
	{347,370,381,384,372,352,333,317,298,285,290,315}, 
	{43.7,50.5,42.8,36.3,28.6,22.9,19.7,18.4, 20.4,23.1,26.1, 33.3}};
	/*
	* [2] Gorshelev, 2014, High spectral resolution ozone absorption
	*	cross-sections. Part 1: Measurementes, data analysis and comparisions
	*	with previous measurements around 293K.
	*
	*	Absolute absorption cross-section in the UV at 295 +-3K,
	*	cm^2 molecule-1 x 10e-20
	*	Wavelenght, nm, row 1; cross-section, row 2.
	*/
	


//};

//{
	// Tabulated values for the US Standard Atmosphere [1], as a function of height
	// (row 1, in km). It stores the anual mean values of the temperature (row 1,
	// in K), presure (row 2, in kPa) and density (row 3, in kg/m^3), for latitude
	// 45ºN. The values have been obtained from [2]. 
	// ----------------------------------------------------------------------------
	// Known Issues:
	// 1-	It seems that the values obtained are given by a simple model fitting 
	//		the data in [1]. This model seems to be a linear model described in [3].
	// 2-	According to [4], there are some typos in the numbers given. It is 
	//		needed to check that the tabulated parameters here correct these errors.
	// 3-	Row 2 are Pa instead of kPa
	//		
	// [1]	US COESA 1976. Standard Atmosphere, 1976. US Government Printing 
	//		Office.
	//		http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf
	// [2]	1976 Standard Atmosphere Calculator. Digital Dutch. Online, last access
	//		Jul 8, 2014.
	//		http://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
	// [3]	The Standard Atmosphere. Atmoscalculator. Online, last access Jul 8, 2014.
	//		http://www.atmosculator.com/The%20Standard%20Atmosphere.html?
	// [4]	AIAA 2010. Guide to Reference and Standard Atmosphere Models. American 
	//		Institute of Aeronautics and Astronautics.
	//		http://www.spacewx.com/Docs/AIAA_G_003C_2010_9-10.pdf
	

// This class models the 1976 US Standard model of atmosphere [1]. It 
// defines the structural properties of an idealized atmosphere, inclu-
// ding the temperature, pressure, and molecular concentration for a 
// given geopotential altitude. It is an idealized, steady state repre-
// sentation of the earth's atmosphere from the surface to 1000 km, as 
// it is assumed to exist during a period of moderate solar activity. 
// The air is assumed to be dry and to obey the perfect gas law and 
// hydrostatic equation, which, taken together, relate temperature, pre-
// ssure and density with geopotential altitude. It should also be noted 
// that since the standard atmosphere model does not include humidity, 
// and since water has a lower molecular weight than air, its presence 
// produces a lower density.
//
// [1]	US COESA 1976. Standard Atmosphere, 1976. US Government Printing 
//		Office.
//		http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf
//
template <typename Float, typename UInt32, typename Mask>
class StandardAtmosphere {
public:

	// ---------------------------------------------------------------------------- 
	// Returns the temperature (in K) for the input altitude z (in Km).
	Float get_temperature(const Float &z) const {
		return Utils::interpolate<Float, UInt32, Mask>(TabulatedValues1[0], TabulatedValues1[1], z);
	}
	// Returns the pressure (in kPa) for the input altitude z (in Km). 
	Float get_pressure(const Float &z) const {
		return Utils::interpolate<Float, UInt32, Mask>(TabulatedValues1[0], TabulatedValues1[2], z);
	}
	// Returns the density (in kg/m^3) for the input altitude z (in Km). 
	Float get_density(const Float &z) const {
		Float density = Utils::interpolate<Float, UInt32, Mask>(TabulatedValues1[0], TabulatedValues1[3], z);
		//printf("Density: %f \n Z: %f\n", density, z);
		return density;
	}

	// Returns the number of density in m^-3 for the input altitude z (in km).
    static Float get_number_density(const Float &z) {
		const Float T_0 = 273.15; //in K
        const Float P_0 = 101.325; // in kPa

        const Float T = Utils::interpolate<Float, UInt32, Mask>(TabulatedValues1[0], TabulatedValues1[1], z);
        const Float P = Utils::interpolate<Float, UInt32, Mask>(TabulatedValues1[0], TabulatedValues1[2], z) * Float(1e-3);

        const Float amg = P / P_0 * T_0 / T;
		// 1 Amagat =  2.6867805e25 m^-3 = 44.615036 mol/m^3
        const Float finalResult = amg * Float(2.6867805e25);

        return enoki::select(
                z > Float(86),
                Float(0.),
                finalResult
        );
	}

	// Returns the amount of ozone at normal temperature 
	// and pressure for the input month in Dobson
	// [1] Ramanathan and R.N. Kulkarni. Height distibution of atmospheric ozone. 1953.
	static Float get_robson_ozone(const Float &z, int month) {
		// 1 Dobson = 2.96 * 10e26 molecules/km^2
		Float proportion = enoki::select(
                z >= Float(0) && z <= Float(9),
                Float(9 / 210.f),
                enoki::select(
                        z > Float(9) && z <= Float(18),
                        Float(14 / 210.f),
                        enoki::select(
                                z > Float(18) && z <= Float(27),
                                Float(111 / 210.f),
                                enoki::select(
                                        z > Float(27) && z <= Float(36),
                                        Float(64 / 210.f),
                                        enoki::select(
                                                z > Float(36) && z < Float(54),
                                                Float(12 / 210.f),
                                                Float(0.)
                                        )
                                )
                        )
                )
        );

		// 1 Dobson = 2.96 * 10e26 molecules/km^2
		return proportion * Float(ozoneDobson[1][month]) * Float(2.96e26);
	}

	// ----------------------------------------------------------------------------
	// Returns the temperature (in K) for the input altitude z (in Km), date (in 
	// days, from 1 to 365), and latitude (in degrees, from 0 to 90).
	Float get_temperature(const Float &z, const int date, const Float latitude) const {
        return get_temperature(z);
    }
	
	// Returns the pressure (in kPa) for the input altitude z (in Km), date (in 
	// days, from 1 to 365), and latitude (in degrees, from 0 to 90).
	Float get_pressure(const Float &z, const int date, const Float latitude) const {
        return get_pressure(z);
    }
	
	// Returns the density (in kg/m^3) for the input altitude z (in Km), date (in 
	// days, from 1 to 365), and latitude (in degrees, from 0 to 90).	
	Float get_density(const Float &z, const int date, const Float latitude) const {
        return get_density(z);
    }

	/*std::vector<float> get_height() {
		return TabulatedValues1[0];
	}*/
	
};
}; //StandardAtmosphere




#endif //_STANDARD_ATMOSPHERE_H_ 
