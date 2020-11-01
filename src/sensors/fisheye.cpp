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
 * - fov
   - |float|
   - Denotes the camera's field of view in degrees---must be between 0 and 180,
     excluding the extremes. Alternatively, it is also possible to specify a
     field of view using the :monosp:`focal_length` parameter.
 * - focal_length
   - |string|
   - Denotes the camera's focal length specified using *35mm* film
     equivalent units. Alternatively, it is also possible to specify a field of
     view using the :monosp:`fov` parameter. See the main description for further
     details. (Default: :monosp:`50mm`)
 * - fov_axis
   - |string|
   - When the parameter :monosp:`fov` is given (and only then), this parameter further specifies
     the image axis, to which it applies.

     1. :monosp:`x`: :monosp:`fov` maps to the :monosp:`x`-axis in screen space.
     2. :monosp:`y`: :monosp:`fov` maps to the :monosp:`y`-axis in screen space.
     3. :monosp:`diagonal`: :monosp:`fov` maps to the screen diagonal.
     4. :monosp:`smaller`: :monosp:`fov` maps to the smaller dimension
        (e.g. :monosp:`x` when :monosp:`width` < :monosp:`height`)
     5. :monosp:`larger`: :monosp:`fov` maps to the larger dimension
        (e.g. :monosp:`y` when :monosp:`width` < :monosp:`height`)

     The default is :monosp:`x`.
 * - near_clip, far_clip
   - |float|
   - Distance to the near/far clip planes. (Default: :monosp:`near_clip=1e-2` (i.e. :monosp:`0.01`)
     and :monosp:`far_clip=1e4` (i.e. :monosp:`10000`))

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/sensor_perspective.jpg
   :caption: The material test ball viewed through a perspective pinhole camera. (:monosp:`fov=28`)
.. subfigure:: ../../resources/data/docs/images/render/sensor_perspective_large_fov.jpg
   :caption: The material test ball viewed through a perspective pinhole camera. (:monosp:`fov=40`)
.. subfigend::
   :label: fig-perspective

This plugin implements a simple idealizied perspective camera model, which
has an infinitely small aperture. This creates an infinite depth of field,
i.e. no optical blurring occurs.

By default, the camera's field of view is specified using a 35mm film
equivalent focal length, which is first converted into a diagonal field
of view and subsequently applied to the camera. This assumes that
the film's aspect ratio matches that of 35mm film (1.5:1), though the
parameter still behaves intuitively when this is not the case.
Alternatively, it is also possible to specify a field of view in degrees
along a given axis (see the :monosp:`fov` and :monosp:`fov_axis` parameters).

The exact camera position and orientation is most easily expressed using the
:monosp:`lookat` tag, i.e.:

.. code-block:: xml

    <sensor type="perspective">
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

        const Float r = enoki::sqrt(position_sample.x() * position_sample.x() + position_sample.y() * position_sample.y()),
                phi = enoki::select(
                        r == Float(0.f),
                        Float(0.f),
                        enoki::select(
                                position_sample.x() < Float(0.f),
                                math::Pi<Float> - enoki::asin(position_sample.y() / r),
                                enoki::asin(position_sample.y() / r)
                                )
                        );

        const Float theta = r * m_aperture_div;
        const Float sinTheta = enoki::sin(theta);

        Vector3f d(sinTheta * enoki::cos(phi), sinTheta * enoki::sin(phi), enoki::cos(theta)); // 1
        //Vector3f d(-sinTheta * enoki::cos(phi), -sinTheta * enoki::sin(phi), enoki::cos(theta)); // 2

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
