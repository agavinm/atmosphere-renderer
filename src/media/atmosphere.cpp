#include <enoki/stl.h>

#include <memory>
#include <mitsuba/atmosphere/Aerosol/BackgroundAerosol.h>
/*#include <mitsuba/atmosphere/Aerosol/DesertDustAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeCleanAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeMineralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarAntarticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarArticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RemoteContinentalAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RuralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/UrbanAerosol.h>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>*/
#include <mitsuba/atmosphere/RayleighScattering.h>
#include <mitsuba/atmosphere/StandardAtmosphere.h>
#include <mitsuba/core/class.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>
#include <string>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class AtmosphereMedium final : public Medium<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction)
    MTS_IMPORT_TYPES(Volume)

    AtmosphereMedium(const Properties &props) : Base(props), m_lat(45.), m_up(0.,1.,0.), m_date(100),
                                                m_earth_radius(6356.766), m_earth_albedo(0.7), m_earth_emission(0.) {
        m_D = props.int_("D", 3);
        if(m_D>3 || m_D<2)
            throw("Invalid dimension D for 'AtmosphereMedium'");

        /*std::string aerosolModel = props.string("aerosol_model", "BackgroundAerosol");
        if (aerosolModel == "DesertDustAerosol")
            m_AerosolModel = std::make_shared<desertDustAerosol::DesertDustAerosol<Float, Spectrum>>();
        else if (aerosolModel == "MaritimeCleanAerosol")
            m_AerosolModel = std::make_shared<maritimeCleanAerosol::MaritimeCleanAerosol<Float, Spectrum>>();
        else if (aerosolModel == "MaritimeMineralAerosol")
            m_AerosolModel = std::make_shared<maritimeMineralAerosol::MaritimeMineralAerosol<Float, Spectrum>>();
        else if (aerosolModel == "PolarAntarticAerosol")
            m_AerosolModel = std::make_shared<polarAntarticAerosol::PolarAntarticAerosol<Float, Spectrum>>();
        else if (aerosolModel == "PolarArticAerosol")
            m_AerosolModel = std::make_shared<polarArticAerosol::PolarArticAerosol<Float, Spectrum>>();
        else if (aerosolModel == "RemoteContinentalAerosol")
            m_AerosolModel = std::make_shared<remoteContinentalAerosol::RemoteContinentalAerosol<Float, Spectrum>>();
        else if (aerosolModel == "RuralAerosol")
            m_AerosolModel = std::make_shared<ruralAerosol::RuralAerosol<Float, Spectrum>>();
        else if (aerosolModel == "UrbanAerosol")
            m_AerosolModel = std::make_shared<urbanAerosol::UrbanAerosol<Float, Spectrum>>();
        else // BackgroundAerosol
            m_AerosolModel = std::make_shared<backgroundAerosol::BackgroundAerosol<Float, Spectrum>>();*/


        m_is_homogeneous = false;
        m_sigmat = props.volume<Volume>("sigma_t", 1.f);

        m_scale = props.float_("scale", 1.0f);
        m_has_spectral_extinction = props.bool_("has_spectral_extinction", true);

        //m_max_density = m_scale * m_sigmat->max();
        m_aabb        = m_sigmat->bbox();


        // init
        init();

        // Precompute values
        //precompute_values(wl); // TODO wl ??
    }

    void init() {
        m_lat_rad = m_lat / 180.*M_PI;
        m_uw = Vector3f(cos(m_lat_rad), sin(m_lat_rad), 0); // TODO Vector2f??
        m_hw = Vector3f(-sin(m_lat_rad), cos(m_lat_rad), 0);
    }

    /*void precompute_values(const int wl) {
        //m_max_extinction = get_extinction(Vector3(0, 65.0995790, 0));

        Vector3f p(0, 0, 0);
        Spectrum extinction;
        m_max_extinction = get_extinction(p,wl);
        Float max_height;

        for (int i = 0; i < 8700; i++) {

            extinction = get_extinction(p,wl);

            if (extinction.get_max() > m_max_extinction.get_max()) {
                m_max_extinction = extinction;
                max_height = p[1];
            }
            p = p + Vector3f(0, 0.01, 0);
        }

    }*/

    UnpolarizedSpectrum
    get_combined_extinction(const MediumInteraction3f & /* mi */,
                            Mask active) const override {
        // TODO: This could be a spectral quantity (at least in RGB mode)
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_max_density;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        auto sigmat = m_scale * m_sigmat->eval(mi, active);
        auto sigmas = sigmat * get_albedo(mi.p, mi.wavelengths[0]); // mi.wavelengths is Color<Float, 1>
        auto sigman = get_combined_extinction(mi, active) - sigmat;
        return { sigmas, sigman, sigmat };
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        return m_aabb.ray_intersect(ray);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("scale", m_scale);
        //callback->put_object("albedo", m_albedo.get());
        callback->put_object("sigma_t", m_sigmat.get());
        Base::traverse(callback);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "AtmosphereMedium[" << std::endl
            //<< "  albedo  = " << string::indent(m_albedo) << std::endl
            << "  sigma_t = " << string::indent(m_sigmat) << std::endl
            << "  scale   = " << string::indent(m_scale) << std::endl
            << "]";
        return oss.str();
    }

    Float get_height(const Vector3f &p) const
    {

        Float earth_radious = 6356.766;
        Vector3f earth_center(0, -earth_radious, 0);
        Float l = (earth_center - p).Size; // TODO Size correcto??
        return (l - earth_radious);

    }

    // Computes the latitude (in degrees) of a point p (in Km) given in local coordinates.
    Float get_latitude(const Vector3f &p) const
    {
        Vector3f pw;
        if (m_D == 3)
            pw = (dot(p, m_hw), dot(p, m_uw), p[2]);
        else if (m_D == 2)
            pw = Vector3f(dot(p, m_hw), dot(p, m_uw), 0); // TODO Vector2f ??

        pw += m_uw; //Check

        return acosf(dot(Vector3f(1., 0., 0), normalize(pw)))*M_1_PI*180.; // TODO Vector2f ??
    }

    bool is_out_of_medium(const Vector3f &p) const {

        Float height = get_height(p);

        if (!(height >= 0 && height <= 86))
            return true;

        return false;
    }

    Spectrum get_extinction(const Vector3f &p, const int wl) const {
        if (is_out_of_medium(p))
            return Spectrum(0.);
        return get_absorption(p, wl) + get_scattering(p, wl);
    }

    Spectrum get_albedo(const Vector3f &p, const int wl) const {
        return get_scattering(p, wl) / get_extinction(p, wl);
    }

    Spectrum get_absorption(const Vector3f &p, const int wl) const {
        //return Spectrum(0.);
        return get_ozone_absorption(p, wl);// +get_aerosol_absorption(p, wl);
    }

    Spectrum get_scattering(const Vector3f &p, const int wl) const {
        return get_rayleigh_scattering(p, wl);
        //return  get_rayleigh_scattering(p, wl) +get_aerosol_scattering(p, wl);
    }

    Spectrum get_rayleigh_scattering(const Vector3f &p, const int wl) const {
        Float h = get_height(p);
        Float lat = get_latitude(p);

        Spectrum cross_section(0.);

        Float density = m_standardAtmosphere.get_number_density(h);

        cross_section = RayleighScattering::get_cross_section<Float, Spectrum>(wl);

        Spectrum finalResult = density * cross_section;

        return finalResult*1e-1;
    }

    Spectrum get_ozone_absorption(const Vector3f &p, const int wl, const int month = 1) const {
        Float h = get_height(p);
        Spectrum cross_section(0.);
        cross_section = RayleighScattering::get_ozone_cross_section<Float, Spectrum>(wl);
        cross_section *= 1e-10;
        Float density = StandardAtmosphere::StandardAtmosphere<Float>::get_robson_ozone(h, 1);//  get_robson_ozone(h, 1);
        Spectrum result = density*cross_section;
        return result;

        //return Spectrum(0.);
    }

    MTS_DECLARE_CLASS()
private:
    ref<Volume> m_sigmat;
    ScalarFloat m_scale;

    ScalarBoundingBox3f m_aabb;
    ScalarFloat m_max_density;


    int m_D;
    // Molecular description
    StandardAtmosphere::StandardAtmosphere<Float> m_standardAtmosphere;
    // Aerosol description
    //std::shared_ptr<GlobalAerosolModel<Float, Spectrum>> m_AerosolModel; // TODO add

    // Sun description (probably useless)
    Float m_sun_phi, m_sun_thita;
    // So that light direction is -[sin(thita)cos(phi), cos(thita), sin(thita)sin(phi)]
    Vector3f m_sun_dir;

    // Local description
    Float m_lat;
    // In degrees, from 0 to 90. It is assumed that the atmosphere is symmetrical wrt the Ecuator.
    Float m_lat_rad;

    Vector3f m_up; // Up-vector in local coordinates.
    Vector3f m_uw, m_hw; // Up and tangent vector in local coordinates.

    int m_date; // In days, where 1 is Jan 1st, and 365 is Dec 31st (no leap year)

    // Earth description
    Float m_earth_radius; // In kilometers
    Spectrum m_earth_albedo;
    Spectrum m_earth_emission;

    // Aerosol description
    // Float m_T; // Turbicity
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmosphereMedium, Medium)
MTS_EXPORT_PLUGIN(AtmosphereMedium, "Atmosphere Medium")
NAMESPACE_END(mitsuba)
