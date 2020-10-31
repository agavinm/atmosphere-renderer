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
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/phase.h>
#include <random>

NAMESPACE_BEGIN(mitsuba)

#define LOG_MODE LogLevel::Debug

// Statistics
//static long long m_stat_nan = 0, m_stat_inf = 0, m_stat_normal = 0;

// Random
static thread_local std::random_device rd;
template<typename ScalarFloat>
static thread_local PCG32<ScalarFloat> m_random_generator(rd());

template <typename Float, typename Spectrum>
class AtmosphereMedium final : public Medium<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Medium, m_phase_function, m_is_homogeneous, m_has_spectral_extinction)
    MTS_IMPORT_TYPES(PhaseFunction, Volume)

    explicit AtmosphereMedium(const Properties &props) : Base(props, true), m_lat(45.), m_up(0.,1.,0.) {
        /*m_D = props.int_("D", 3);
        if(m_D>3 || m_D<2)
            Throw("Invalid dimension D for 'AtmosphereMedium'");*/

        // Phase functions
        std::string molecular_phase_function = props.string("molecular_phase_function", "cha");
        if (molecular_phase_function != "rayleigh")
            m_phase_function = PluginManager::instance()->create_object<PhaseFunction>(Properties("cha"));
        else
            m_phase_function = PluginManager::instance()->create_object<PhaseFunction>(Properties("rayleigh"));

        Properties hg("hg");
        hg.set_float("g", 0.76f);
        m_aerosol_phase_function = PluginManager::instance()->create_object<PhaseFunction>(hg);

        m_month = props.int_("month", 1);
        if(m_month < 1 || m_month > 12)
            Throw("Invalid month number: Valid month is from 1 to 12");
        m_month--;

        // Aerosol model
        std::string aerosolModel = props.string("aerosol_model", "");
        if (aerosolModel == "BackgroundAerosol")
            m_AerosolModel = std::make_shared<backgroundAerosol::BackgroundAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "DesertDustAerosol")
            m_AerosolModel = std::make_shared<desertDustAerosol::DesertDustAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "MaritimeCleanAerosol")
            m_AerosolModel = std::make_shared<maritimeCleanAerosol::MaritimeCleanAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "MaritimeMineralAerosol")
            m_AerosolModel = std::make_shared<maritimeMineralAerosol::MaritimeMineralAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "PolarAntarticAerosol")
            m_AerosolModel = std::make_shared<polarAntarticAerosol::PolarAntarticAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "PolarArticAerosol")
            m_AerosolModel = std::make_shared<polarArticAerosol::PolarArticAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "RemoteContinentalAerosol")
            m_AerosolModel = std::make_shared<remoteContinentalAerosol::RemoteContinentalAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "RuralAerosol")
            m_AerosolModel = std::make_shared<ruralAerosol::RuralAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else if (aerosolModel == "UrbanAerosol")
            m_AerosolModel = std::make_shared<urbanAerosol::UrbanAerosol<Float, UInt32, Mask, Spectrum, Wavelength>>();
        else
            m_AerosolModel = nullptr;

        m_is_homogeneous = false;
        //m_albedo = props.volume<Volume>("albedo", 0.75f);
        //m_sigmat = props.volume<Volume>("sigma_t", 1.f);
        //Log(LOG_MODE, "Sigma_t bbox: \"%s\"", m_sigmat.get());

        //m_scale = props.float_("scale", 1.0f);
        m_has_spectral_extinction = props.bool_("has_spectral_extinction", true);

        //m_max_density = m_scale * m_sigmat->max();
        //m_aabb        = m_sigmat->bbox();

        m_earth_radius = props.float_("earth_radius");
        Log(LOG_MODE, "Initialized Earth with \"%s\" radius.", m_earth_radius);

        m_earth_scale = ScalarFloat(6356.766) / m_earth_radius;
        Log(LOG_MODE, "Initialized Earth scale as \"%s\".", m_earth_scale);

        ScalarPoint3f earth_center = props.point3f("earth_center");
        worldToLocal = ScalarTransform4f(ScalarMatrix4f(ScalarPoint4f(1, 0, 0, 0),
                                                        ScalarPoint4f(0, 1, 0, 0),
                                                        ScalarPoint4f(0, 0, 1, 0),
                                                        ScalarPoint4f(earth_center.x(), earth_center.y(), earth_center.z(), 1))).inverse();
        Log(LOG_MODE, "Initialized worldToLocal transform as \"%s\"", worldToLocal);
        
        const ScalarFloat box_radius = m_earth_radius + ScalarFloat(160.) / m_earth_scale;
        m_aabb = ScalarBoundingBox3f(ScalarPoint3f(earth_center.x() - box_radius, earth_center.y() - box_radius, earth_center.z() - box_radius),
                                     ScalarPoint3f(earth_center.x() + box_radius, earth_center.y() + box_radius, earth_center.z() + box_radius));
        std::ostringstream oss;
        oss << m_aabb;
        Log(LOG_MODE, "Initialized Bounding Box as \"%s\"", oss.str());

        // init
        m_lat_rad = m_lat / ScalarFloat(180.*M_PI);
        m_uw = ScalarVector3f(cos(m_lat_rad), sin(m_lat_rad), 0);
        m_hw = ScalarVector3f(-sin(m_lat_rad), cos(m_lat_rad), 0);

        // Precompute values
        ScalarVector3f p(0, m_earth_radius, 0);
        ScalarFloat extinction = get_max_extinction(p);
        m_max_extinction = extinction;

        for (int i = 0; i < 8699; i++) {
            p = p + ScalarVector3f(0, ScalarFloat(0.01) / m_earth_scale, 0);

            extinction = get_max_extinction(p);

            if (enoki::all(extinction > m_max_extinction)) {
                m_max_extinction = extinction;
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

        /*const auto keys = m_rayleighScattering.ScatteringCrossSection[0];
        const Float x = 453;
        const auto idx1 = Utils::get_low(keys, x);
        const auto idx2 = Utils::get_low2<Float, UInt32, Mask>(keys, x);
        Log(Info, "idx1 = \"%s\"", idx1);
        Log(Info, "idx2 = \"%s\"", idx2);*/

        /*PCG32<ScalarFloat> m_random_generator;
        Log(Info, "random1 = \"%s\"", m_random_generator.next_float32());
        Log(Info, "random2 = \"%s\"", m_random_generator.next_float32());
        Log(Info, "random3 = \"%s\"", m_random_generator.next_float32());*/

        //Log(Info, "random = \"%s\"", m_random_generator<ScalarFloat>.next_float32());
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
        const auto p = worldToLocal.transform_affine(mi.p);
        const Mask out_of_medium = is_out_of_medium(p);
        auto sigmas = UnpolarizedSpectrum(0.);
        auto sigmat = UnpolarizedSpectrum(0.);
        if (!enoki::all(out_of_medium)) {
            sigmas = enoki::select(
                    out_of_medium,
                    Spectrum(0.),
                    get_scattering(p, mi.wavelengths)
            );
            sigmat = enoki::select(
                    out_of_medium,
                    Spectrum(0.),
                    get_absorption(p, mi.wavelengths) + sigmas
            );
        }
        //UnpolarizedSpectrum sigman = 0.f;

        // Before:
        //auto sigmat = m_scale * m_sigmat->eval(mi, active); // TODO: Use correct sigma_t
        //auto sigmas = sigmat * get_albedo(worldToLocal.transform_affine(mi.p), mi.wavelengths); // TODO: 0 is first wavelength, select the correct wavelength
        const auto sigman = get_combined_extinction(mi, active) - sigmat;

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
            return m_phase_function.get();

        const auto p = worldToLocal.transform_affine(mi.p);
        const Mask out_of_medium = is_out_of_medium(p);
        if (enoki::all(out_of_medium))
            return m_phase_function.get();

        // Get scattering
        const Spectrum rayleigh_scattering = enoki::select(
                out_of_medium,
                Spectrum(0.),
                get_rayleigh_scattering(p, mi.wavelengths)
                ); // Molecular scattering (sigmas^m)
        const Spectrum aerosol_scattering = enoki::select(
                out_of_medium,
                Spectrum(0.),
                get_aerosol_scattering(p, mi.wavelengths)
                ); // Aerosol scattering (sigmas^a)

        // Get scattering probability
        const ScalarFloat total = Utils::get_first<Float, ScalarFloat>(rayleigh_scattering[0]) + Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]);
        if (total == ScalarFloat(0))
            return m_phase_function.get();
        const ScalarFloat rayleigh_scattering_prob = Utils::get_first<Float, ScalarFloat>(rayleigh_scattering[0]) / total;
        //const ScalarFloat aerosol_scattering_prob = Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]) / total;

        //Log(Info, "rayleigh_scattering = \"%s\"; aerosol_scattering = \"%s\"", Utils::get_first<Float, ScalarFloat>(rayleigh_scattering[0]), Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]));
        //Log(Info, "rayleigh_scattering_prob = \"%s\"; aerosol_scattering_prob = \"%s\"", rayleigh_scattering_prob, aerosol_scattering_prob);

        // Russian roulette

        //Log(Info, "rayleigh_scattering_prob = \"%s\"; random = \"%s\"", rayleigh_scattering_prob, m_random_generator<ScalarFloat>.next_float32());
        if (m_random_generator<ScalarFloat>.next_float32() < rayleigh_scattering_prob)
            return m_phase_function.get(); // P(molecular) = [0, rayleigh_scattering_prob)
        else
            return m_aerosol_phase_function.get(); // P(aerosol) = [rayleigh_scattering_prob,  1) //TODO: Use all wl
    }

    /**
     * Computes the geopotential height of a point p given in local coordinates.
     * @param p
     * @return
     */
    Float get_height(const Vector3f &p) const {
        const Float l = norm(p);
        const Float height = (l - m_earth_radius) * m_earth_scale;
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
        //if (m_D == 3)
            pw = (dot(p, m_hw), dot(p, m_uw), p[2]);
        /*else if (m_D == 2)
            pw = Vector3f(dot(p, m_hw), dot(p, m_uw), 0);*/

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
        const Float height = get_height(p);
        return !(height >= Float(0) && height <= Float(86));
    }

    ScalarFloat get_max_extinction(const ScalarVector3f &p) const {
        if (enoki::all(is_out_of_medium(p)))
            return ScalarFloat(0.);

        // get_height
        const ScalarFloat h = (norm(p) - m_earth_radius) * m_earth_scale;

        const StandardAtmosphere::StandardAtmosphere aux_atmosphere;
        auto finalResult = ScalarFloat(0.);

        // get_ozone_absorption
        finalResult += m_rayleighScattering.get_ozone_cross_section() * 1e-10f * aux_atmosphere.get_robson_ozone<float, ScalarUInt32, ScalarMask>(h, m_month);

        // get_rayleigh_scattering
        finalResult += aux_atmosphere.get_number_density<float, ScalarUInt32, ScalarMask>(h) * m_rayleighScattering.get_cross_section() * 1e-1f;

        if (m_AerosolModel == nullptr)
            return finalResult;

        // get_aerosol_absorption
        finalResult += m_AerosolModel->get_density_float(h) * m_AerosolModel->get_absorption();

        // get_aerosol_scattering
        finalResult += m_AerosolModel->get_density_float(h) * m_AerosolModel->get_scattering();

        return finalResult;
    }

    Spectrum get_absorption(const Vector3f &p, const Wavelength &wl) const {
        if (m_AerosolModel == nullptr)
            return get_ozone_absorption(p, wl);
        else
            return get_ozone_absorption(p, wl) + get_aerosol_absorption(p, wl);
        /*else {
            auto ozone_absorption = get_ozone_absorption(p, wl), aerosol_absorption = get_aerosol_absorption(p, wl);
            Log(Info, "ozone_absorption = \"%s\"; aerosol_absorption = \"%s\"", Utils::get_first<Float, ScalarFloat>(ozone_absorption[0]), Utils::get_first<Float, ScalarFloat>(aerosol_absorption[0]));
            return ozone_absorption + aerosol_absorption;
        }*/
    }

    Spectrum get_scattering(const Vector3f &p, const Wavelength &wl) const {
        if (m_AerosolModel == nullptr)
            return get_rayleigh_scattering(p, wl);
        else
            return get_rayleigh_scattering(p, wl) + get_aerosol_scattering(p, wl);
        /*else {
            auto rayleigh_scattering = get_rayleigh_scattering(p, wl), aerosol_scattering = get_aerosol_scattering(p, wl);
            Log(Info, "rayleigh_scattering = \"%s\"; aerosol_scattering = \"%s\"", Utils::get_first<Float, ScalarFloat>(rayleigh_scattering[0]), Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]));
            return rayleigh_scattering + aerosol_scattering;
        }*/
    }

    Spectrum get_rayleigh_scattering(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);
        //Log(LOG_MODE, "Height of ray collision from the earth center \"%s\" km.", h);
        //Float lat = get_latitude(p);
        //Log(LOG_MODE, "Latitude of ray collision: \"%s\"", lat);

        const Spectrum cross_section = m_rayleighScattering.get_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl);

        const auto density = m_standardAtmosphere.get_number_density<Float, UInt32, Mask>(h);

        const Spectrum finalResult = density * cross_section;

        //Log(Info, "get_rayleigh_scattering = \"%s\"", finalResult);

        return finalResult * 1e-1;
    }

    Spectrum get_ozone_absorption(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);

        const Spectrum cross_section = m_rayleighScattering.get_ozone_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl) * Float(1e-10);

        const auto density = m_standardAtmosphere.get_robson_ozone<Float, UInt32, Mask>(h, m_month);//  get_robson_ozone(h, 1);

        const Spectrum result = cross_section * density;

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
        const Float h = get_height(p);

        const Spectrum cross_section = m_AerosolModel->get_absorption(wl);

        const Float density = m_AerosolModel->get_density(h);

        const Spectrum final_result = cross_section * density;

        //Log(Info, "get_aerosol_absorption = \"%s\"", final_result);

        return final_result;
    }

    Spectrum get_aerosol_scattering(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);

        const Spectrum cross_section = m_AerosolModel->get_scattering(wl);

        const Float density = m_AerosolModel->get_density(h);

        const Spectrum final_result = cross_section * density;

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
    ref<PhaseFunction> m_aerosol_phase_function;

    ScalarBoundingBox3f m_aabb;
    //ScalarFloat m_max_density;


    //int m_D;
    // Molecular description
    const StandardAtmosphere::StandardAtmosphere m_standardAtmosphere;
    const RayleighScattering::RayleighScattering m_rayleighScattering;
    // Aerosol description
    std::shared_ptr<GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>> m_AerosolModel{};

    // Sun description (probably useless)
    //Float m_sun_phi, m_sun_thita;
    // So that light direction is -[sin(thita)cos(phi), cos(thita), sin(thita)sin(phi)]
    //Vector3f m_sun_dir;

    // Local description
    ScalarFloat m_lat;
    // In degrees, from 0 to 90. It is assumed that the atmosphere is symmetrical wrt the Ecuator.
    ScalarFloat m_lat_rad;

    ScalarVector3f m_up; // Up-vector in local coordinates.
    ScalarVector3f m_uw, m_hw; // Up and tangent vector in local coordinates.

    //int m_date; // In days, where 1 is Jan 1st, and 365 is Dec 31st (no leap year)
    int m_month;

    // Earth description
    ScalarFloat m_earth_radius; // In kilometers
    ScalarFloat m_earth_scale; // From 0 to 1 (where 1 is real radius = m_earth_radius)
    //Spectrum m_earth_albedo;
    //Spectrum m_earth_emission;

    ScalarTransform4f worldToLocal;

    // Aerosol description
    // Float m_T; // Turbicity

    // Total atmosphere values
    ScalarFloat m_max_extinction;
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmosphereMedium, Medium)
MTS_EXPORT_PLUGIN(AtmosphereMedium, "Atmosphere Medium")
NAMESPACE_END(mitsuba)
