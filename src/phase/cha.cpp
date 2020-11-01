#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/atmosphere/RayleighScattering.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-cha:

Chandrasekhar phase function (:monosp:`cha`)
-----------------------------------------------

This plugin implements the phase function model proposed by Chandrasekhar.
It does not have any parameters.

*/
template <typename Float, typename Spectrum>
class ChandrasekharPhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(PhaseFunction, m_flags)
    MTS_IMPORT_TYPES(PhaseFunctionContext)

    ChandrasekharPhaseFunction(const Properties &props) : Base(props), m_constant(ScalarFloat(3) / (ScalarFloat(16) * math::Pi<ScalarFloat>)) {
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
        Spectrum result = 0.;

        for (size_t i = 0; i < result.Size; i++) {
            const auto gamma = m_rayleighScattering.gamma<Float, UInt32, Mask>(mi.wavelengths[i]);

            const auto cosine = dot(wo, mi.wi);
            result[i] = (m_constant / (Float(1) + Float(2) * gamma)) *
                    (Float(1) + Float(3) * gamma + (Float(1) - gamma) * (cosine * cosine));
        }

        return result;
    }

    std::string to_string() const override { return "ChandrasekharPhaseFunction[]"; }

    MTS_DECLARE_CLASS()
private:
    const ScalarFloat m_constant;
    const RayleighScattering::RayleighScattering m_rayleighScattering;
};

MTS_IMPLEMENT_CLASS_VARIANT(ChandrasekharPhaseFunction, PhaseFunction)
MTS_EXPORT_PLUGIN(ChandrasekharPhaseFunction, "Chandrasekhar phase function")
NAMESPACE_END(mitsuba)