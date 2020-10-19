#ifndef _GLOBAL_ATMOSPHERIC_AEROSOL_H_
#define _GLOBAL_ATMOSPHERIC_AEROSOL_H_

// This class models the Global Aerosol Model described in [1], that
// describes the distribution of aerosol in the atmosphere.
//
// [1]	AIAA 1999. Guide to Global Aerosol Models (GAM). American 
//		Institute of Aeronautics and Astronautics.
//		http://www.spacewx.com/Docs/AIAA-656-598.pdf


template <typename Float, typename Spectrum, typename Wavelength>
class GlobalAerosolModel {
public:

    virtual Spectrum get_scattering() const = 0;
    virtual Spectrum get_scattering(const Wavelength &wl) const = 0;
    virtual Spectrum get_absorption() const = 0;
    virtual Spectrum get_absorption(const Wavelength &wl) const = 0;
	virtual Float get_density(Float z) const = 0;
}; //GlobalAerosolModel

#endif //_ATMOSPHERIC_AEROSOL_H_