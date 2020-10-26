#ifndef _RAYLEIGH_SCATTERING_H_
#define _RAYLEIGH_SCATTERING_H_

#include <vector>
#include <mitsuba/atmosphere/Utils.h>

namespace RayleighScattering {

	// Tabulated values presented in [1] (Table 1), including the King correction factor (F_k),
	// the depolarization factor (\rho_n) and the \gamma term used in the phase function, for
	// wavelengths from 200 to 1000 nm.
	//
	// [1] A. Bucholtz 1995. Rayleigh-scattering calculations for the terrestrial atmosphere. Applied
	//	   Optics, Vol. 34(15).
	//	   http://augerlal.lal.in2p3.fr/pmwiki/uploads/Bucholtz.pdf
	const static std::vector<std::vector<float>> LambdaDependentValues =
		{{200, 205, 210, 215, 220, 225, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000},
		{1.080, 1.077, 1.074, 1.072, 1.070, 1.068, 1.066, 1.064, 1.062, 1.060, 1.059, 1.057, 1.056, 1.055, 1.055, 1.054, 1.053, 1.053, 1.052, 1.052, 1.052, 1.051, 1.051, 1.051, 1.050, 1.049, 1.049, 1.048, 1.048, 1.048, 1.048, 1.047, 1.047, 1.047, 1.047, 1.047},
		{0.04545, 0.04384, 0.04221, 0.04113, 0.04004, 0.03895, 0.03785, 0.03675, 0.03565, 0.03455, 0.03400, 0.03289, 0.03233, 0.03178, 0.03178, 0.03122, 0.03066, 0.03066, 0.03010, 0.03010, 0.03010, 0.02955, 0.02955, 0.02955, 0.02899, 0.02842, 0.02842, 0.02786, 0.02786, 0.02786, 0.02786, 0.02730, 0.02730, 0.02730, 0.02730, 0.02730},
		{0.02326, 0.02241, 0.02156, 0.02100 ,0.02043, 0.01986, 0.01930, 0.01872, 0.01815, 0.01758, 0.01729, 0.01672, 0.01643, 0.01614, 0.01614, 0.01586, 0.01557, 0.01557, 0.01528, 0.01528, 0.01528, 0.01499, 0.01499, 0.01499, 0.01471, 0.01442, 0.01442, 0.01413, 0.01413, 0.01413, 0.01413, 0.01384, 0.01384, 0.01384, 0.01384, 0.01384}};

	// Returns the gamma value (Eq. 13 [1]) for the Rayleigh phase function including
	// the depolarization term, as defined in [2]. The values are obtained from the 
	// tabulated values in [1] (Table 1).
	//
	// [1] A. Bucholtz 1995. Rayleigh-scattering calculations for the terrestrial atmosphere. Applied
	//	   Optics, Vol. 34(15).
	//	   http://augerlal.lal.in2p3.fr/pmwiki/uploads/Bucholtz.pdf
	// [2] S. Chandrasekhar 1960. Radiative Transfer. Dover.
    /*template <typename Float, typename Spectrum>
	const inline Spectrum gamma() {
		return Spectrum(std::vector<Float>(LambdaDependentValues<Float>[0], LambdaDependentValues<Float>[0]+36),&(LambdaDependentValues<Float>[3][0]));
	};*/

	//need to be done
    template <typename Float, typename UInt32, typename Mask>
	inline Float gamma(const Float &wl) {
		return Utils::interpolate<Float, UInt32, Mask>(LambdaDependentValues[0], LambdaDependentValues[3], wl);
	};


	// Tabulated scattering cross section values for wavelengths from 200 to 4000 nm [1]. Row 2 stores
	// the values measured in [1] (Table 2), while Row 1 stores the values obtained using the fitting 
	// described in [1], Eq. 8 with the values stored in [1], Table 3.
	//
	// [1] A. Bucholtz 1995. Rayleigh-scattering calculations for the terrestrial atmosphere. Applied
	//	   Optics, Vol. 34(15).
	//	   http://augerlal.lal.in2p3.fr/pmwiki/uploads/Bucholtz.pdf
    const static std::vector<std::vector<float>> ScatteringCrossSection =
		{{ 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 900, 1000, 1100, 1200, 1300, 1500, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500, 4000},
		{3.59637e-25, 2.90136e-25, 2.36416e-25, 1.94405e-25, 1.61196e-25, 1.34686e-25, 1.13332e-25, 9.59864e-26, 8.17885e-26, 7.00833e-26, 6.03686e-26, 5.22556e-26, 4.54407e-26, 3.96848e-26, 3.47985e-26, 3.06302e-26, 2.70583e-26, 2.39843e-26, 2.1328e-26, 1.90238e-26, 1.70178e-26, 1.52652e-26, 1.3729e-26, 1.23783e-26, 1.11871e-26, 1.01335e-26, 9.1991e-27, 8.36827e-27, 7.62767e-27, 6.96591e-27, 6.37324e-27, 6.48194e-27, 5.98211e-27, 5.52926e-27, 5.11822e-27, 4.74446e-27, 4.40401e-27, 4.09337e-27, 3.80949e-27, 3.54966e-27, 3.31148e-27, 3.09283e-27, 2.89182e-27, 2.70679e-27, 2.53624e-27, 2.37883e-27, 2.23338e-27, 2.09881e-27, 1.97417e-27, 1.85858e-27, 1.75129e-27, 1.65158e-27, 1.55883e-27, 1.47245e-27, 1.39195e-27, 1.31684e-27, 1.24669e-27, 1.18113e-27, 1.1198e-27, 1.06237e-27, 1.00856e-27, 6.1988e-28, 4.01061e-28, 2.7049e-28, 1.88793e-28, 1.35621e-28, 7.50751e-29, 7.50751e-29, 5.74996e-29, 4.47567e-29, 3.53404e-29, 2.8264e-29, 2.28652e-29, 1.54211e-29, 1.07634e-29, 7.73202e-30, 5.69228e-30, 4.28017e-30, 2.26358e-30, 1.30358e-30}, 
		{3.612e-25, 2.836e-25, 2.269e-25, 1.841e-25, 1.515e-25, 1.259e-25, 1.056e-25, 8.939e-26, 7.614e-26, 6.534e-26, 5.642e-26, 4.903e-26, 4.279e-26, 3.752e-26, 3.307e-26, 2.924e-26, 2.598e-26, 2.317e-26, 2.071e-26, 1.858e-26, 1.673e-26, 1.509e-26, 1.366e-26, 1.239e-26, 1.127e-26, 1.027e-26, 9.378e-27, 8.583e-27, 7.871e-27, 7.232e-27, 6.656e-27, 6.138e-27, 5.669e-27, 5.245e-27, 4.86e-27, 4.509e-27, 4.189e-27, 3.897e-27, 3.63e-27, 3.385e-27, 3.161e-27, 2.956e-27, 2.767e-27, 2.593e-27, 2.432e-27, 2.284e-27, 2.147e-27, 2.02e-27, 1.903e-27, 1.793e-27, 1.692e-27, 1.597e-27, 1.51e-27, 1.428e-27, 1.351e-27, 1.28e-27, 1.213e-27, 1.15e-27, 1.092e-27, 1.037e-27, 9.854e-28, 6.129e-28, 4.01e-28, 2.734e-28, 1.927e-28, 1.398e-28, 1.038e-28, 7.872e-29, 6.077e-29, 4.766e-29, 3.79e-29, 3.052e-29, 2.485e-29, 1.697e-29, 1.197e-29, 8.691e-30, 6.46e-30, 4.901e-30, 2.645e-30, 1.55e-30}};
	//


	/*
	 * 1e-20
	 Absolute absorption cross-sections in the UV at 295 ± 3 K, cm2 molecule−1 × 10−20
	 */
    const static std::vector<std::vector<float>> ozoneUVcrossSection =
        { { 244,248,253,257,289,296,302 },
		{ 946e-20,105.1e-20,1120e-20,1107e-20,151e-20,61.1e-20,29.6e-20} };

	/*
	*	Absoluted absorption cross-section in the near 380 nm at
	*	295+-3K, cm^2 molecue^-1 x 10e-23 gave by [2].
	*	Wavelenght, nm, row 1; cross-section, row 2.
	*/
    const static std::vector<std::vector<float>> ozoneNearCrossSection =
        { { 365, 405,455 },
		{ 4.9e-23,1.46e-23,20.6e-23 } };

	/*
	*	Absoluted absorption cross-section in the visible at 295 +-3K,
	*	cm^2 molecule^-1 x 10e-21 gave by [2].
	*	Wavelength, nm, row 1; cross-section, row 2.
	*/
    const static std::vector<std::vector<float>> ozoneVisibleCrossSection =
        { { 543,576,594,604,611,632 },
		{ 3.08e-21,4.70e-21,4.63e-21, 5.10e-21, 4.54e-21,3.36e-21 } };

	/*
	*	Absolute absorption cross-section in the near infra red at 295+-3K,
	*	cm^2 molecule^-1 x 10e-22 gave by [2].
	*	Wavelenght, nm, row 1; cross-section, row 2.
	*/
    const static std::vector<std::vector<float>> ozoneNIRCrossSection =
        { { 748,755,760,765,770,779,802,817,853,877,889,898,933,944,991,1046 },
		{ 4.38e-22,3.22e-22,2.77e-22,2.53e-22,2.49e-22,3.15e-22,
		1.45e-22,2.20e-22,1.46e-22,0.377e-22, 0.510e-22,0.638e-22,
		0.162e-22,0.424e-22,0.407e-22,0.0773e-22 } };

	// Returns the Rayleigh scattering cross section, from the values obtained in [1]. This values are 
	// computed for standard air (15ºC, 1013.25mb), and are given, for each wavelength, in cm^2. 
	// Wavelengths are given in nm. 
	//
	// [1] A. Bucholtz 1995. Rayleigh-scattering calculations for the terrestrial atmosphere. Applied
	//	   Optics, Vol. 34(15).
	//	   http://augerlal.lal.in2p3.fr/pmwiki/uploads/Bucholtz.pdf
    float get_cross_section() {
        return Utils::max(ScatteringCrossSection[2]);
    }

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
	Spectrum get_cross_section(const Wavelength &wl) {
		// interpolates with one given value for wavelenght
        Spectrum s(0.);

        for (size_t i = 0; i < wl.Size; i++)
            s[i] = Utils::interpolate<Float, UInt32, Mask>(ScatteringCrossSection[0], ScatteringCrossSection[2], wl[i]);

        return s;
	}

    float get_ozone_cross_section() {
        const float ozone_cross_sections[4] = {Utils::max(ozoneUVcrossSection[1]),
                                         Utils::max(ozoneNearCrossSection[1]),
                                         Utils::max(ozoneVisibleCrossSection[1]),
                                         Utils::max(ozoneNIRCrossSection[1])};
        float result = ozone_cross_sections[0];
        for (int i = 1; i < 4; i++) {
            if (ozone_cross_sections[i] > result)
                result = ozone_cross_sections[i];
        }

        return result;
    }

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
	Spectrum get_ozone_cross_section(const Wavelength &wl) {
        Spectrum s(0.);

        for (size_t i = 0; i < wl.Size; i++) {
            s[i] = enoki::select(
                    wl[i] >= Float(244) && wl[i] <= Float(302),
                    Utils::interpolate<Float, UInt32, Mask>(ozoneUVcrossSection[0], ozoneUVcrossSection[1], wl[i]),
                    enoki::select(
                            wl[i] > Float(302) && wl[i] < Float(365),
                            Utils::fast_interpolate<Float>(ozoneUVcrossSection[0][6], ozoneNearCrossSection[0][0], wl[i],
                                                           ozoneUVcrossSection[1][6], ozoneNearCrossSection[1][0]),
                            enoki::select(
                                    wl[i] >= Float(365) && wl[i] <= Float(455),
                                    Utils::interpolate<Float, UInt32, Mask>(ozoneNearCrossSection[0], ozoneNearCrossSection[1], wl[i]),
                                    enoki::select(
                                            wl[i] > Float(455) && wl[i] < Float(543),
                                            Utils::fast_interpolate<Float>(ozoneNearCrossSection[0][2], ozoneVisibleCrossSection[0][0], wl[i],
                                                                           ozoneNearCrossSection[1][2], ozoneVisibleCrossSection[1][0]),
                                            enoki::select(
                                                    wl[i] >= Float(543) && wl[i] <= Float(632),
                                                    Utils::interpolate<Float, UInt32, Mask>(ozoneVisibleCrossSection[0], ozoneVisibleCrossSection[1], wl[i]),
                                                    enoki::select(
                                                            wl[i] > Float(632) && wl[i] < Float(748),
                                                            Utils::fast_interpolate<Float>(ozoneVisibleCrossSection[0][5], ozoneNIRCrossSection[0][0], wl[i],
                                                                                           ozoneVisibleCrossSection[1][5], ozoneNIRCrossSection[1][0]),
                                                            enoki::select(
                                                                    wl[i] >= Float(748) && wl[i] <= Float(1046),
                                                                    Utils::interpolate<Float, UInt32, Mask>(ozoneNIRCrossSection[0], ozoneNIRCrossSection[1], wl[i]),
                                                                    Float(0.)
                                                            )
                                                    )
                                            )
                                    )
                            )
                    )
            );
        }

        return s;
	}

	

}; //RayleighScattering
#endif //_RAYLEIGH_SCATTERING_H_
