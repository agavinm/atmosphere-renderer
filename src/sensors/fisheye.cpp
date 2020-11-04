#include <mitsuba/render/sensor.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/bbox.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _sensor-fisheye:

Angular fisheye camera (:monosp:`fisheye`)
--------------------------------------------------

.. pluginparameters::

 * - to_world
   - |transform|
   - Specifies an optional camera-to-world transformation.
     (Default: none (i.e. camera space = world space))
 * - aperture
   - |float|
   - Denotes the camera's aperture in degrees---must be between 0 and 360,
     excluding 0.

This plugin implements an angular fisheye camera model, based on http://paulbourke.net/dome/fisheye/

The exact camera position and orientation is most easily expressed using the
:monosp:`lookat` tag, i.e.:

.. code-block:: xml

    <sensor type="fisheye">
        <transform name="to_world">
            <!-- Move and rotate the camera so that looks from (1, 1, 1) to (1, 2, 1)
                and the direction (0, 0, 1) points "up" in the output image -->
            <lookat origin="1, 1, 1" target="1, 2, 1" up="0, 0, 1"/>
        </transform>
    </sensor>

 */

template <typename Float, typename Spectrum>
class FisheyeCamera final : public ProjectiveCamera<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(ProjectiveCamera, m_world_transform, m_needs_sample_3,
                    m_film, m_sampler, m_resolution, m_shutter_open,
                    m_shutter_open_time, m_near_clip, m_far_clip)
    MTS_IMPORT_TYPES()

    // =============================================================
    //! @{ \name Constructors
    // =============================================================

    FisheyeCamera(const Properties &props) : Base(props) {
        ScalarVector2i size = m_film->size();

        if (m_world_transform->has_scale())
            Throw("Scale factors in the camera-to-world transformation are not allowed!");

        m_aperture = props.float_("aperture");
        if (m_aperture < 1 || m_aperture > 360)
            Throw("Invalid aperture angle: Valid angle is from 1 to 360");
        m_aperture = enoki::deg_to_rad(m_aperture);
        m_aperture_div = m_aperture / ScalarFloat(2.f);
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Sampling methods (Sensor interface)
    // =============================================================

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &position_sample,
                                          const Point2f & /*aperture_sample*/,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] = sample_wavelength<Float, Spectrum>(wavelength_sample);
        Ray3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;

        const Float x = position_sample.x() * Float(2.f) - Float(1.f),
                y = position_sample.y() * Float(2.f) - Float(1.f);

        const Float r = enoki::sqrt(x * x + y * y),
                phi = enoki::select(
                        r == Float(0.f),
                        Float(0.f),
                        enoki::select(
                                x < Float(0.f),
                                math::Pi<Float> - enoki::asin(y / r),
                                enoki::asin(y / r)
                                )
                        );

        const Float theta = r * m_aperture_div;
        const Float sinTheta = enoki::sin(theta);

        Vector3f d(sinTheta * enoki::cos(phi), sinTheta * enoki::sin(phi), enoki::cos(theta));

        auto trafo = m_world_transform->eval(ray.time, active);
        ray.o = trafo.translation();
        ray.d = trafo * d;
        ray.update();

        wav_weight = enoki::select(
                r > Float(1.f),
                Spectrum(0.f),
                wav_weight);

        return std::make_pair(ray, wav_weight);
    }

    ScalarBoundingBox3f bbox() const override {
        return m_world_transform->translation_bounds();
    }

    //! @}
    // =============================================================

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        Base::parameters_changed(keys);
    }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "FisheyeCamera[" << std::endl
            //<< "  x_fov = " << m_x_fov << "," << std::endl
            << "  near_clip = " << m_near_clip << "," << std::endl
            << "  far_clip = " << m_far_clip << "," << std::endl
            << "  film = " << indent(m_film) << "," << std::endl
            << "  sampler = " << indent(m_sampler) << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  shutter_open = " << m_shutter_open << "," << std::endl
            << "  shutter_open_time = " << m_shutter_open_time << "," << std::endl
            << "  world_transform = " << indent(m_world_transform) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ScalarFloat m_aperture, m_aperture_div;
};

MTS_IMPLEMENT_CLASS_VARIANT(FisheyeCamera, ProjectiveCamera)
MTS_EXPORT_PLUGIN(FisheyeCamera, "Fisheye Camera");
NAMESPACE_END(mitsuba)
