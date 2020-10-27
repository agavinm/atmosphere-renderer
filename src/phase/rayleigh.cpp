#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-rayleigh:

Rayleigh phase function (:monosp:`rayleigh`)
-----------------------------------------------

This plugin implements the Rayleigh phase function model.
It does not have any parameters.

*/
template <typename Float, typename Spectrum>
class RayleighPhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(PhaseFunction, m_flags)
    MTS_IMPORT_TYPES(PhaseFunctionContext)

    RayleighPhaseFunction(const Properties &props) : Base(props), m_constant(ScalarFloat(3) / (ScalarFloat(16) * ScalarFloat(M_PI))) {
        m_flags = +PhaseFunctionFlags::Anisotropic;
    }

    std::pair<Vector3f, Float> sample(const PhaseFunctionContext & /* ctx */,
                                      const MediumInteraction3f &, const Point2f &sample,
                                      Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        auto wo  = warp::square_to_uniform_sphere(sample);
        auto pdf = warp::square_to_uniform_sphere_pdf(wo);
        return std::make_pair(wo, pdf);
    }

    Spectrum eval(const PhaseFunctionContext & /* ctx */, const MediumInteraction3f &mi,
               const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);

        const auto cosine = dot(wo, mi.wi);
        return m_constant * (Float(1) + cosine * cosine);
    }

    std::string to_string() const override { return "RayleighPhaseFunction[]"; }

    MTS_DECLARE_CLASS()
private:
    const ScalarFloat m_constant;
};

MTS_IMPLEMENT_CLASS_VARIANT(RayleighPhaseFunction, PhaseFunction)
MTS_EXPORT_PLUGIN(RayleighPhaseFunction, "Rayleigh phase function")
NAMESPACE_END(mitsuba)