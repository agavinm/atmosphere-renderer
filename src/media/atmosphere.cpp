#include <enoki/stl.h>

#include <memory>
#include <mitsuba/atmosphere/Aerosol/BackgroundAerosol.h>
#include <mitsuba/atmosphere/Aerosol/DesertDustAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeCleanAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeMineralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarAntarticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarArticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RemoteContinentalAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RuralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/UrbanAerosol.h>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>
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

#define LOG_MODE LogLevel::Debug

// Statistics
static long long m_stat_nan = 0, m_stat_inf = 0, m_stat_normal = 0;

template <typename Float, typename Spectrum>
class AtmosphereMedium final : public Medium<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction)
    MTS_IMPORT_TYPES(Volume)

    AtmosphereMedium(const Properties &props) : Base(props), m_lat(45.), m_up(0.,1.,0.), m_date(100),
                                                m_earth_albedo(0.7), m_earth_emission(0.) {
        m_D = props.int_("D", 3);
        if(m_D>3 || m_D<2)
            throw("Invalid dimension D for 'AtmosphereMedium'");

        std::string aerosolModel = props.string("aerosol_model", "");
        if (aerosolModel == "BackgroundAerosol")
            m_AerosolModel = std::make_shared<backgroundAerosol::BackgroundAerosol<Float, Spectrum>>();
        else if (aerosolModel == "DesertDustAerosol")
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
        else
            m_AerosolModel = nullptr;

        m_is_homogeneous = false;
        //m_albedo = props.volume<Volume>("albedo", 0.75f);
        m_sigmat = props.volume<Volume>("sigma_t", 1.f);
        Log(LOG_MODE, "Sigma_t bbox: \"%s\"", m_sigmat.get());

        m_scale = props.float_("scale", 1.0f);
        m_has_spectral_extinction = props.bool_("has_spectral_extinction", true);

        //m_max_density = m_scale * m_sigmat->max();
        //m_aabb        = m_sigmat->bbox();

        m_earth_radius = props.float_("earth_radius_km");
        Log(LOG_MODE, "Initialized Earth with \"%s\" radius.", std::to_string(m_earth_radius));

        m_earth_scale = Float(6356.766) / m_earth_radius;
        Log(LOG_MODE, "Initialized Earth scale as \"%s\".", std::to_string(m_earth_scale));

        Point3f earth_center = props.point3f("earth_center");
        worldToLocal = Transform4f(Matrix4f(Point4f(1, 0, 0, 0),
                                   Point4f(0, 1, 0, 0),
                                   Point4f(0, 0, 1, 0),
                                   Point4f(earth_center.x(), earth_center.y(), earth_center.z(), 1))).inverse();
        Log(LOG_MODE, "Initialized worldToLocal transform as \"%s\"", worldToLocal);
        
        const Float box_radius = m_earth_radius + 160.0f / m_earth_scale;
        m_aabb = ScalarBoundingBox3f(Point3f(earth_center.x() - box_radius, earth_center.y() - box_radius, earth_center.z() - box_radius),
                                     Point3f(earth_center.x() + box_radius, earth_center.y() + box_radius, earth_center.z() + box_radius));
        std::ostringstream oss;
        oss << m_aabb;
        Log(LOG_MODE, "Initialized Bounding Box as \"%s\"", oss.str());

        // init
        init();

        // Precompute values
        //precompute_values(wl); // TODO wl ??
    }

    void init() {
        m_lat_rad = m_lat / 180.*M_PI;
        m_uw = Vector3f(cos(m_lat_rad), sin(m_lat_rad), 0);
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
    get_combined_extinction(const MediumInteraction3f &mi,
                            Mask active) const override {
        // TODO: This could be a spectral quantity (at least in RGB mode)
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_sigmat->eval(mi) * m_scale;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        auto sigmat = m_scale * m_sigmat->eval(mi, active);
        Log(LOG_MODE, "World vector p \"%s\"", mi.p);
        auto sigmas = sigmat * get_albedo(worldToLocal.transform_affine(mi.p), mi.wavelengths[0]); // TODO: 0 is first wavelength, select the correct wavelength
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

    /**
     * Computes the geopotential height of a point p given in local coordinates.
     * @param p
     * @return
     */
    Float get_height(const Vector3f &p) const {
        Float l = norm(p);
        Float aux = (l - m_earth_radius) * m_earth_scale;
        Log(LOG_MODE, "Vector p \"%s\"", p);
        Log(LOG_MODE, "Pre-Height of ray collision from the earth center \"%s\".", std::to_string(l));
        Log(LOG_MODE, "Height of ray collision from the earth center \"%s\" km.", std::to_string(aux));
        return aux;
        //return (l - m_earth_radius);
    }

    /**
     * Computes the latitude (in degrees) of a point p given in local coordinates.
     * @param p
     * @return
     */
    Float get_latitude(const Vector3f &p) const {
        Vector3f pw;
        if (m_D == 3)
            pw = (dot(p, m_hw), dot(p, m_uw), p[2]);
        else if (m_D == 2)
            pw = Vector3f(dot(p, m_hw), dot(p, m_uw), 0);

        pw += m_uw; //TODO Check

        return acos(dot(Vector3f(1, 0, 0), normalize(pw))) * Float(M_1_PI) * Float(180);
    }

    /**
     * If this distance between center of the Earth and position of the ray is
     * less than earth_radius+athmosphere_size, ray is inside of medium.
     * Otherwise, not.
     * @param p
     * @return
     */
    bool is_out_of_medium(const Vector3f &p) const {
        Float height = get_height(p);

        if (!(height >= Float(0) && height <= Float(86)))
            return true;

        return false;
    }

    Spectrum get_extinction(const Vector3f &p, const int wl) const {
        if (is_out_of_medium(p))
            return Spectrum(0.);
        return get_absorption(p, wl) + get_scattering(p, wl);
    }

    Spectrum get_albedo(const Vector3f &p, const int wl) const {
        Spectrum scattering = get_scattering(p, wl);
        Log(LOG_MODE, "Scattering: \"%s\"", scattering);
        Spectrum extinction = get_extinction(p, wl);
        Log(LOG_MODE, "Extinction: \"%s\"", extinction);
        Spectrum albedo = scattering / extinction;
        for (auto &a : albedo) { // TODO: Why is nan or inf ??
            if (isnan(a)) {
                a = 0;
                m_stat_nan++;
            }
            else if (isinf(a)) {
                a = 1;
                m_stat_inf++;
            }
            else {
                m_stat_normal++;
            }
        }
        Log(LOG_MODE, "Albedo: \"%s\"", albedo);
        return albedo;
    }

    Spectrum get_absorption(const Vector3f &p, const int wl) const {
        if (m_AerosolModel == nullptr)
            return get_ozone_absorption(p, wl);
        else
            return get_ozone_absorption(p, wl) + get_aerosol_absorption(p, wl);
    }

    Spectrum get_scattering(const Vector3f &p, const int wl) const {
        if (m_AerosolModel == nullptr)
            return get_rayleigh_scattering(p, wl);
        else
            return get_rayleigh_scattering(p, wl) + get_aerosol_scattering(p, wl);
    }

    Spectrum get_rayleigh_scattering(const Vector3f &p, const int wl) const {
        Float h = get_height(p);
        //Log(LOG_MODE, "Height of ray collision from the earth center \"%s\" km.", std::to_string(h));
        Float lat = get_latitude(p);
        //Log(LOG_MODE, "Latitude of ray collision: \"%s\"", std::to_string(lat));

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


    // ----------------------------------------------------------------------------
    // Aerosols scattering and absorption functions
    Spectrum get_aerosol_absorption(const Vector3f &p, const int wl) const {
        Float h = get_height(p);
        Spectrum cross_section(0.);

        Float density = m_AerosolModel->get_density(h);

        cross_section = m_AerosolModel->get_absorption(wl);
        Spectrum final_result = cross_section * density;

        return final_result;
    }

    Spectrum get_aerosol_scattering(const Vector3f &p, const int wl) const {
        Float h = get_height(p);
        Spectrum cross_section(0.);

        Float density = m_AerosolModel->get_density(h);

        cross_section = m_AerosolModel->get_scattering(wl);
        Spectrum final_result = cross_section * density;

        return final_result;
    }

    ~AtmosphereMedium() {
        Float total = Float(100) / (m_stat_nan + m_stat_inf + m_stat_normal);
        Log(Info, "Number of NaN in get_albedo(): \"%s\" (\"%s\"\"%\")", std::to_string(m_stat_nan), std::to_string(m_stat_nan * total));
        Log(Info, "Number of Inf in get_albedo(): \"%s\" (\"%s\"\"%\")", std::to_string(m_stat_inf), std::to_string(m_stat_inf * total));
        Log(Info, "Number of Normal in get_albedo(): \"%s\" (\"%s\"\"%\")", std::to_string(m_stat_normal), std::to_string(m_stat_normal * total));
    }

    MTS_DECLARE_CLASS()
private:
    ref<Volume> m_sigmat;
    ScalarFloat m_scale;

    ScalarBoundingBox3f m_aabb;
    //ScalarFloat m_max_density;


    int m_D;
    // Molecular description
    StandardAtmosphere::StandardAtmosphere<Float> m_standardAtmosphere;
    // Aerosol description
    std::shared_ptr<GlobalAerosolModel<Float, Spectrum>> m_AerosolModel;

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
    Float m_earth_scale; // From 0 to 1 (where 1 is real radius = m_earth_radius)
    Spectrum m_earth_albedo;
    Spectrum m_earth_emission;

    Transform4f worldToLocal;

    // Aerosol description
    // Float m_T; // Turbicity
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmosphereMedium, Medium)
MTS_EXPORT_PLUGIN(AtmosphereMedium, "Atmosphere Medium")
NAMESPACE_END(mitsuba)
