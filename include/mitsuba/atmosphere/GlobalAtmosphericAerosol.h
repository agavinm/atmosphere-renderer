#ifndef _GLOBAL_ATMOSPHERIC_AEROSOL_H_
#define _GLOBAL_ATMOSPHERIC_AEROSOL_H_

#include <array>
#include <mitsuba/atmosphere/Utils.h>

// This class models the Global Aerosol Model described in [1], that
// describes the distribution of aerosol in the atmosphere.
//
// [1]	AIAA 1999. Guide to Global Aerosol Models (GAM). American 
//		Institute of Aeronautics and Astronautics.
//		http://www.spacewx.com/Docs/AIAA-656-598.pdf

#define UNIT_CONVERSION 1e-6f // cross section values are in m^2

template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
class GlobalAerosolModel {
public:
    std::array<float, 100001> m_absorption{}, m_scattering{};
    float m_max_absorption, m_max_scattering;

protected:
    explicit GlobalAerosolModel(const std::array<std::array<float, 1001>, 3> &tabulatedValues) {
        m_max_absorption = tabulatedValues[1][0] * UNIT_CONVERSION;
        m_max_scattering = tabulatedValues[2][0] * UNIT_CONVERSION;
        for (size_t wl = 0; wl <= tabulatedValues[0][0] * 100; wl++) {
            m_absorption[wl] = tabulatedValues[1][0] * UNIT_CONVERSION;
            m_scattering[wl] = tabulatedValues[2][0] * UNIT_CONVERSION;
        }

        for (size_t i = 1; i < tabulatedValues[0].size(); i++) {
            if (tabulatedValues[1][i] * UNIT_CONVERSION > m_max_absorption)
                m_max_absorption = tabulatedValues[1][i] * UNIT_CONVERSION;
            if (tabulatedValues[2][i] * UNIT_CONVERSION > m_max_scattering)
                m_max_scattering = tabulatedValues[2][i] * UNIT_CONVERSION;

            for (size_t wl = tabulatedValues[0][i - 1] * 100.f + 1; wl < tabulatedValues[0][i] * 100.f; wl++) {
                m_absorption[wl] = Utils::fast_interpolate<float>(tabulatedValues[0][i - 1] * 100.f,
                                                                  tabulatedValues[0][i] * 100.f,
                                                                  wl,
                                                                  tabulatedValues[1][i - 1] * UNIT_CONVERSION,
                                                                  tabulatedValues[1][i] * UNIT_CONVERSION);
                m_scattering[wl] = Utils::fast_interpolate<float>(tabulatedValues[0][i - 1] * 100.f,
                                                                  tabulatedValues[0][i] * 100.f,
                                                                  wl,
                                                                  tabulatedValues[2][i - 1] * UNIT_CONVERSION,
                                                                  tabulatedValues[2][i] * UNIT_CONVERSION);
            }
            m_absorption[size_t(tabulatedValues[0][i] * 100.f)] = tabulatedValues[1][i] * UNIT_CONVERSION;
            m_scattering[size_t(tabulatedValues[0][i] * 100.f)] = tabulatedValues[2][i] * UNIT_CONVERSION;
        }
    }

public:
    [[nodiscard]] float get_absorption() const {
        return m_max_absorption;
    }

    Spectrum get_absorption(const Wavelength &wl) const {
        Spectrum s(0.);

        for (size_t i = 0; i < wl.Size; i++)
            s[i] = Utils::get_<Float, UInt32, Mask>(wl[i] * Float(100.f), m_absorption);

        return s;
    }

    [[nodiscard]] float get_scattering() const {
        return m_max_scattering;
    }

    Spectrum get_scattering(const Wavelength &wl) const {
        Spectrum s(0.);

        for (size_t i = 0; i < wl.Size; i++)
            s[i] = Utils::get_<Float, UInt32, Mask>(wl[i] * Float(100.f), m_scattering);

        return s;
    }

	virtual Float get_density(const Float &z) const = 0;

    [[nodiscard]] virtual float get_density_float(const float &z) const = 0;
}; //GlobalAerosolModel

#endif //_ATMOSPHERIC_AEROSOL_H_