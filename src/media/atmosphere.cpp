#include <enoki/stl.h>

#include <memory>
/*#include <mitsuba/atmosphere/Aerosol/BackgroundAerosol.h>
#include <mitsuba/atmosphere/Aerosol/DesertDustAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeCleanAerosol.h>
#include <mitsuba/atmosphere/Aerosol/MaritimeMineralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarAntarticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/PolarArticAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RemoteContinentalAerosol.h>
#include <mitsuba/atmosphere/Aerosol/RuralAerosol.h>
#include <mitsuba/atmosphere/Aerosol/UrbanAerosol.h>*/
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
#include <mitsuba/render/phase.h>
#include <string>
#include <random>

NAMESPACE_BEGIN(mitsuba)

#define LOG_MODE LogLevel::Debug

// Statistics
//static long long m_stat_nan = 0, m_stat_inf = 0, m_stat_normal = 0;

// Random generator
static thread_local std::random_device rd;
static thread_local std::mt19937 mt(rd());
static thread_local std::uniform_real_distribution<float> random_dist(0.0f, 1.0f);

template <typename Float, typename Spectrum>
class AtmosphereMedium final : public Medium<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction)
    MTS_IMPORT_TYPES(PhaseFunction, Volume)

    AtmosphereMedium(const Properties &props) : Base(props, true), m_lat(45.), m_up(0.,1.,0.), m_date(100),
                                                m_earth_albedo(0.7), m_earth_emission(0.) {
        m_D = props.int_("D", 3);
        if(m_D>3 || m_D<2)
            throw("Invalid dimension D for 'AtmosphereMedium'");

        // Phase functions
        m_molecular_phase_function = PluginManager::instance()->create_object<PhaseFunction>(Properties("cha"));
        Properties hg("hg");
        hg.set_float("g", 0.76f);
        m_aerosol_phase_function = PluginManager::instance()->create_object<PhaseFunction>(hg);

        // Aerosol model
        /*std::string aerosolModel = props.string("aerosol_model", "");
        if (aerosolModel == "BackgroundAerosol")
            m_AerosolModel = std::make_shared<backgroundAerosol::BackgroundAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "DesertDustAerosol")
            m_AerosolModel = std::make_shared<desertDustAerosol::DesertDustAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "MaritimeCleanAerosol")
            m_AerosolModel = std::make_shared<maritimeCleanAerosol::MaritimeCleanAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "MaritimeMineralAerosol")
            m_AerosolModel = std::make_shared<maritimeMineralAerosol::MaritimeMineralAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "PolarAntarticAerosol")
            m_AerosolModel = std::make_shared<polarAntarticAerosol::PolarAntarticAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "PolarArticAerosol")
            m_AerosolModel = std::make_shared<polarArticAerosol::PolarArticAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "RemoteContinentalAerosol")
            m_AerosolModel = std::make_shared<remoteContinentalAerosol::RemoteContinentalAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "RuralAerosol")
            m_AerosolModel = std::make_shared<ruralAerosol::RuralAerosol<Float, Spectrum, Wavelength>>();
        else if (aerosolModel == "UrbanAerosol")
            m_AerosolModel = std::make_shared<urbanAerosol::UrbanAerosol<Float, Spectrum, Wavelength>>();
        else*/
            m_AerosolModel = nullptr;

        m_is_homogeneous = false;
        //m_albedo = props.volume<Volume>("albedo", 0.75f);
        //m_sigmat = props.volume<Volume>("sigma_t", 1.f);
        //Log(LOG_MODE, "Sigma_t bbox: \"%s\"", m_sigmat.get());

        //m_scale = props.float_("scale", 1.0f);
        m_has_spectral_extinction = props.bool_("has_spectral_extinction", true);

        //m_max_density = m_scale * m_sigmat->max();
        //m_aabb        = m_sigmat->bbox();

        m_earth_radius = props.float_("earth_radius_km");
        Log(LOG_MODE, "Initialized Earth with \"%s\" radius.", m_earth_radius);

        m_earth_scale = Float(6356.766) / m_earth_radius;
        Log(LOG_MODE, "Initialized Earth scale as \"%s\".", m_earth_scale);

        Point3f earth_center = props.point3f("earth_center");
        worldToLocal = Transform4f(Matrix4f(Point4f(1, 0, 0, 0),
                                   Point4f(0, 1, 0, 0),
                                   Point4f(0, 0, 1, 0),
                                   Point4f(earth_center.x(), earth_center.y(), earth_center.z(), 1))).inverse();
        Log(LOG_MODE, "Initialized worldToLocal transform as \"%s\"", worldToLocal);
        
        const Float box_radius = m_earth_radius + Float(160.) / m_earth_scale;
        m_aabb = BoundingBox3f(Point3f(earth_center.x() - box_radius, earth_center.y() - box_radius, earth_center.z() - box_radius),
                                     Point3f(earth_center.x() + box_radius, earth_center.y() + box_radius, earth_center.z() + box_radius));
        std::ostringstream oss;
        oss << m_aabb;
        Log(LOG_MODE, "Initialized Bounding Box as \"%s\"", oss.str());

        // init
        m_lat_rad = m_lat / Float(180.*M_PI);
        m_uw = Vector3f(cos(m_lat_rad), sin(m_lat_rad), 0);
        m_hw = Vector3f(-sin(m_lat_rad), cos(m_lat_rad), 0);

        // Precompute values
        Vector3f p(0, m_earth_radius, 0);
        Spectrum extinction = get_extinction(p, -1);
        m_max_extinction = extinction[0];

        for (int i = 0; i < 8699; i++) {
            p = p + Vector3f(0, Float(0.01) / m_earth_scale, 0);

            extinction = get_extinction(p, -1);

            if (all(extinction[0] > m_max_extinction)) {
                m_max_extinction = extinction[0];
            }
        }

        Log(Info, "Max extinction = \"%s\"", m_max_extinction);

        /*const Mask msk = Float(0.4) < Float(1);
        Log(Info, "0.4 < 1; msk = \"%s\", count(msk) = \"%s\", count(!msk) = \"%s\"", msk, count(msk), count(!msk));

        std::vector<float> array = {100, 20, 10, 30, 40, 50};
        auto x = Float(25);

        UInt32 j = binary_search(
                0,
                array.size() - 2,
                [&](UInt32 index) {
                    return !(gather<Float>(array.data(), index) <= x && gather<Float>(array.data(), index + 1) >= x);
                } // It search until return false
        );

        Log(Info, "array = \"%s\"", array);
        Log(Info, "x = \"%s\"", x);
        Log(Info, "j = \"%s\"", j); // array[j] <= x && array[j+1] >= x*/

        /*const std::vector<float> keys = {1, 2, 3, 4, 5, 6}, values = {0.1, 0.2, 0.3, 0.74, 0.5, 0.6};
        const auto x = Utils::interpolate<Float, UInt32>(keys, values, Float(3.5));
        Log(Info, "x = \"%s\"", x);*/

        /*const auto keys = RayleighScattering::ScatteringCrossSection[0];
        const Float x = 453;
        const auto idx1 = Utils::get_low(keys, x);
        const auto idx2 = Utils::get_low2<Float, UInt32, Mask>(keys, x);
        Log(Info, "idx1 = \"%s\"", idx1);
        Log(Info, "idx2 = \"%s\"", idx2);*/
    }

    UnpolarizedSpectrum
    get_combined_extinction(const MediumInteraction3f & /* mi */,
                            Mask active) const override {
        // TODO: This could be a spectral quantity (at least in RGB mode)
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_max_extinction;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        /*for (const auto &a : mi.p) {
            if (isinf(a) || isnan(a)) // TODO Why ??
                Log(Info, "World vector p \"%s\"", mi.p);
                //return { Spectrum(0.), Spectrum(0.), Spectrum(0.) };
        }*/
        Log(LOG_MODE, "World vector p \"%s\"", mi.p);

        // After:
        auto p = worldToLocal.transform_affine(mi.p);
        const Mask out_of_medium = is_out_of_medium(p);
        UnpolarizedSpectrum sigmas = enoki::select(
                out_of_medium,
                Spectrum(0.),
                get_scattering(p, mi.wavelengths)
        );
        UnpolarizedSpectrum sigmat = enoki::select(
                out_of_medium,
                Spectrum(0.),
                get_absorption(p, mi.wavelengths) + sigmas
        );
        //UnpolarizedSpectrum sigman = 0.f;

        // Before:
        //auto sigmat = m_scale * m_sigmat->eval(mi, active); // TODO: Use correct sigma_t
        //auto sigmas = sigmat * get_albedo(worldToLocal.transform_affine(mi.p), mi.wavelengths); // TODO: 0 is first wavelength, select the correct wavelength
        auto sigman = get_combined_extinction(mi, active) - sigmat;

        /*Log(Info, "p \"%s\"", p);
        Log(Info, "Sigma_t \"%s\"", sigmat);
        Log(Info, "Sigma_s \"%s\"", sigmas);
        Log(Info, "Sigma_n \"%s\"", sigman);*/
        /*for (const auto &a : sigmat) {
            if (isnan(a) || isinf(a))
                Log(Info, "Sigma_t \"%s\"", sigmat);
        }*/

        return { sigmas, sigman, sigmat };
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        return m_aabb.ray_intersect(ray);
    }

    void traverse(TraversalCallback *callback) override {
        //callback->put_parameter("scale", m_scale);
        //callback->put_object("albedo", m_albedo.get());
        //callback->put_object("sigma_t", m_sigmat.get());
        Base::traverse(callback);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "AtmosphereMedium[" << std::endl
            //<< "  albedo  = " << string::indent(m_albedo) << std::endl
            //<< "  sigma_t = " << string::indent(m_sigmat) << std::endl
            //<< "  scale   = " << string::indent(m_scale) << std::endl
            << "]";
        return oss.str();
    }

    const PhaseFunction *phase_function(const MediumInteraction3f &mi) const override {
        if (m_AerosolModel == nullptr)
            return m_molecular_phase_function.get();

        // Get scattering
        auto p = worldToLocal.transform_affine(mi.p);
        auto rayleigh_scattering = get_rayleigh_scattering(p, mi.wavelengths); // Molecular scattering (sigmas^m)
        auto aerosol_scattering = get_aerosol_scattering(p, mi.wavelengths); // Aerosol scattering (sigmas^a)

        // Get scattering probability
        Float rayleigh_scattering_prob = 0.;//, aerosol_scattering_prob = 0.;
        Float size = 0;
        for (size_t i = 0; i < rayleigh_scattering.Size; i++) {
            Float total = rayleigh_scattering[i] + aerosol_scattering[i];
            const Mask msk = total != Float(0.);
            size += enoki::select(
                    msk,
                    Float(1),
                    Float(0)
            );
            rayleigh_scattering_prob += enoki::select(
                    msk,
                    rayleigh_scattering[i] / total,
                    Float(0)
            );
            //rayleigh_scattering_prob += rayleigh_scattering[i] / total;
            //aerosol_scattering_prob += aerosol_scattering[i] / total;
        }
        rayleigh_scattering_prob /= size;
        //aerosol_scattering_prob /= size;

        Log(Info, "size = \"%s\"; rayleigh_scattering_prob = \"%s\"; hmean = \"%s\"", size, rayleigh_scattering_prob, enoki::hmean(rayleigh_scattering_prob));

        //Log(Info, "rayleigh_scattering = \"%s\"; aerosol_scattering = \"%s\"", rayleigh_scattering, aerosol_scattering);
        //Log(Info, "rayleigh_scattering_prob = \"%s\"; aerosol_scattering_prob = \"%s\"", rayleigh_scattering_prob, aerosol_scattering_prob);

        // Russian roulette
        if (random_dist(mt) < enoki::hmean(rayleigh_scattering_prob))
            return m_molecular_phase_function.get(); // P(molecular) = [0, rayleigh_scattering_prob) //TODO: Use all wl
        else
            return m_aerosol_phase_function.get(); // P(aerosol) = [rayleigh_scattering_prob,  1) //TODO: Use all wl
    } // TODO: Check

    /**
     * Computes the geopotential height of a point p given in local coordinates.
     * @param p
     * @return
     */
    Float get_height(const Vector3f &p) const {
        Float l = norm(p);
        Float height = (l - m_earth_radius) * m_earth_scale;
        Log(LOG_MODE, "Vector p \"%s\"", p);
        Log(LOG_MODE, "Pre-Height of ray collision from the earth center \"%s\".", l);
        Log(LOG_MODE, "Height of ray collision from the earth center \"%s\" km.", height);
        return height;
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
    Mask is_out_of_medium(const Vector3f &p) const {
        Float height = get_height(p);
        return !(height >= Float(0) && height <= Float(86));
    }

    Spectrum get_extinction(const Vector3f &p, const Wavelength &wl) const {
        return enoki::select(
                is_out_of_medium(p),
                Spectrum(0.),
                get_absorption(p, wl) + get_scattering(p, wl)
                );
    }

    /*Spectrum get_albedo(const Vector3f &p, const Wavelength &wl) const {
        if (is_out_of_medium(p))
            return Spectrum(0.); // TODO: Is it correct?

        Spectrum scattering = get_scattering(p, wl);
        Spectrum extinction = get_extinction(p, wl);
        Spectrum albedo = scattering / extinction;
        for (auto &a : albedo) { // TODO: Why is nan or inf ??
            if (isnan(a)) { // Because Extinction == Scattering == 0
                a = 0;
                m_stat_nan++;
                Log(Info, "[albedo = NaN] Scattering: \"%s\"", scattering);
                Log(Info, "[albedo = NaN] Extinction: \"%s\"", extinction);
            }
            else if (isinf(a)) { // Because Extinction == 0 and Scattering != 0
                a = 1;
                m_stat_inf++;
                Log(Info, "[albedo = Inf] Scattering: \"%s\"", scattering);
                Log(Info, "[albedo = Inf] Extinction: \"%s\"", extinction);
            }
            else {
                m_stat_normal++;
            }
        }
        Log(LOG_MODE, "Albedo: \"%s\"", albedo);
        return albedo;

        return get_scattering(p, wl) / get_extinction(p, wl);
    }*/

    Spectrum get_absorption(const Vector3f &p, const Wavelength &wl) const {
        if (m_AerosolModel == nullptr)
            return get_ozone_absorption(p, wl);
        else
            return get_ozone_absorption(p, wl) + get_aerosol_absorption(p, wl);
    }

    Spectrum get_scattering(const Vector3f &p, const Wavelength &wl) const {
        if (m_AerosolModel == nullptr)
            return get_rayleigh_scattering(p, wl);
        else
            return get_rayleigh_scattering(p, wl) + get_aerosol_scattering(p, wl);
    }

    Spectrum get_rayleigh_scattering(const Vector3f &p, const Wavelength &wl) const {
        Float h = get_height(p);
        //Log(LOG_MODE, "Height of ray collision from the earth center \"%s\" km.", h);
        //Float lat = get_latitude(p);
        //Log(LOG_MODE, "Latitude of ray collision: \"%s\"", lat);

        Spectrum cross_section(0.);

        Float density = m_standardAtmosphere.get_number_density(h);

        if (wl == -1)
            cross_section = RayleighScattering::get_cross_section<Spectrum>();
        else
            cross_section = RayleighScattering::get_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl);

        Spectrum finalResult = density * cross_section;

        //Log(Info, "get_rayleigh_scattering = \"%s\"", finalResult);

        return finalResult*1e-1;
    }

    Spectrum get_ozone_absorption(const Vector3f &p, const Wavelength &wl, const int month = 1) const {
        Float h = get_height(p);

        Spectrum cross_section(0.);

        if (wl == -1)
            cross_section = RayleighScattering::get_ozone_cross_section<Spectrum>();
        else
            cross_section = RayleighScattering::get_ozone_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl);

        cross_section *= 1e-10;

        Float density = StandardAtmosphere::StandardAtmosphere<Float, UInt32, Mask>::get_robson_ozone(h, 1);//  get_robson_ozone(h, 1);

        Spectrum result = cross_section * density;

        //Log(Info, "h = \"%s\"", h);
        //Log(Info, "cross_section = \"%s\"", cross_section);
        //Log(Info, "cross_section2 = \"%s\"", cross_section);
        //Log(Info, "density = \"%s\"", density);
        //Log(Info, "result = \"%s\"", result);
        return result;

        //return Spectrum(0.);
    }


    // ----------------------------------------------------------------------------
    // Aerosols scattering and absorption functions
    Spectrum get_aerosol_absorption(const Vector3f &p, const Wavelength &wl) const {
        Float h = get_height(p);
        Spectrum cross_section(0.);

        Float density = m_AerosolModel->get_density(h);

        if (wl == -1)
            cross_section = m_AerosolModel->get_absorption();
        else
            cross_section = m_AerosolModel->get_absorption(wl);

        Spectrum final_result = cross_section * density;

        //Log(Info, "get_aerosol_absorption = \"%s\"", final_result);

        return final_result;
    }

    Spectrum get_aerosol_scattering(const Vector3f &p, const Wavelength &wl) const {
        Float h = get_height(p);
        Spectrum cross_section(0.);

        Float density = m_AerosolModel->get_density(h);

        if (wl == -1)
            cross_section = m_AerosolModel->get_scattering();
        else
            cross_section = m_AerosolModel->get_scattering(wl);

        Spectrum final_result = cross_section * density;

        //Log(Info, "get_aerosol_scattering = \"%s\"", final_result);

        return final_result;
    }

    /*~AtmosphereMedium() {
        Float total = Float(100) / (m_stat_nan + m_stat_inf + m_stat_normal);
        Log(Info, "Number of NaN in get_albedo(): \"%s\" (\"%s\"\"%\")", m_stat_nan, m_stat_nan * total);
        Log(Info, "Number of Inf in get_albedo(): \"%s\" (\"%s\"\"%\")", m_stat_inf, m_stat_inf * total);
        Log(Info, "Number of Normal in get_albedo(): \"%s\" (\"%s\"\"%\")", m_stat_normal), m_stat_normal * total);
    }*/

    MTS_DECLARE_CLASS()
private:
    //ref<Volume> m_sigmat;
    //ScalarFloat m_scale;
    ref<PhaseFunction> m_molecular_phase_function, m_aerosol_phase_function;

    BoundingBox3f m_aabb;
    //ScalarFloat m_max_density;


    int m_D;
    // Molecular description
    StandardAtmosphere::StandardAtmosphere<Float, UInt32, Mask> m_standardAtmosphere;
    // Aerosol description
    std::shared_ptr<GlobalAerosolModel<Float, Spectrum, Wavelength>> m_AerosolModel;

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

    // Total atmosphere values
    Float m_max_extinction;
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmosphereMedium, Medium)
MTS_EXPORT_PLUGIN(AtmosphereMedium, "Atmosphere Medium")
NAMESPACE_END(mitsuba)
