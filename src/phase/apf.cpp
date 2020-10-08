#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/atmosphere/RayleighScattering.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-apf:

Atmosphere phase function (:monosp:`apf`)
-----------------------------------------------

This phase function simulates completely the Earth atmosphere.
It does not have any parameters.

*/
template <typename Float, typename Spectrum>
class AtmospherePhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(PhaseFunction, m_flags)
    MTS_IMPORT_TYPES(PhaseFunctionContext)

    AtmospherePhaseFunction(const Properties &props) : Base(props), m_constant(Float(3) / (Float(16) * M_PI)) {
        m_flags = +PhaseFunctionFlags::Anisotropic; // TODO ??
    }

    std::pair<Vector3f, Float> sample(const PhaseFunctionContext & /* ctx */,
                                      const MediumInteraction3f &mi, const Point2f &sample,
                                      Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        auto wo  = warp::square_to_uniform_sphere(sample);
        auto pdf = warp::square_to_uniform_sphere_pdf(wo);
        return std::make_pair(wo, pdf);
    } // TODO ??

    Float eval(const PhaseFunctionContext & /* ctx */, const MediumInteraction3f &mi,
               const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);

        const auto gamma = RayleighScattering::gamma<Float>(mi.wavelengths[0]); // TODO: 0 is first wavelength, select the correct wavelength

        return (m_constant / (Float(1) + Float(2) * gamma)) *
               (Float(1) + Float(3) * gamma + (Float(1) - gamma) * pow(dot(wo, mi.wi), Float(2)));
    } // TODO: Only is F^m

    std::string to_string() const override { return "AtmospherePhaseFunction[]"; }

    MTS_DECLARE_CLASS()
private:
    const Float m_constant;
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmospherePhaseFunction, PhaseFunction)
MTS_EXPORT_PLUGIN(AtmospherePhaseFunction, "Atmosphere phase function")
NAMESPACE_END(mitsuba)