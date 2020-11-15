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

//#define PRINT_DATA
//#define DEBUG_INFO

#ifdef PRINT_DATA
#include <iostream>
#include <fstream>
#endif

NAMESPACE_BEGIN(mitsuba)

#ifdef DEBUG_INFO
static float min_z = 86.f, max_z = 0.f;
#endif

// Random
static thread_local std::random_device rd;
template<typename ScalarFloat>
static thread_local PCG32<ScalarFloat> m_random_generator(rd());

template <typename Float, typename Spectrum>
class AtmosphereMedium final : public Medium<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Medium, m_phase_function, m_is_homogeneous, m_has_spectral_extinction)
    MTS_IMPORT_TYPES(PhaseFunction, Volume)

    explicit AtmosphereMedium(const Properties &props) : Base(props, true) {
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

        m_turbidity = props.float_("turbidity", 1.f);
        if(m_turbidity < 0.f)
            Throw("Invalid turbidity scale parameter");

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
        m_has_spectral_extinction = props.bool_("has_spectral_extinction", true);

        m_earth_radius = props.float_("earth_radius");

        m_earth_scale = ScalarFloat(6356.766) / m_earth_radius;

        ScalarPoint3f earth_center = props.point3f("earth_center");
        worldToLocal = ScalarTransform4f(ScalarMatrix4f(ScalarPoint4f(1, 0, 0, 0),
                                                        ScalarPoint4f(0, 1, 0, 0),
                                                        ScalarPoint4f(0, 0, 1, 0),
                                                        ScalarPoint4f(earth_center.x(), earth_center.y(), earth_center.z(), 1))).inverse();
        
        const ScalarFloat box_radius = m_earth_radius + ScalarFloat(87.f) / m_earth_scale;
        m_aabb = ScalarBoundingBox3f(ScalarPoint3f(earth_center.x() - box_radius, earth_center.y() - box_radius, earth_center.z() - box_radius),
                                     ScalarPoint3f(earth_center.x() + box_radius, earth_center.y() + box_radius, earth_center.z() + box_radius));

        // Precompute values
        ScalarVector3f p(0, m_earth_radius, 0);
        ScalarFloat extinction = get_max_extinction(p);
        m_max_extinction = extinction;

#ifdef PRINT_DATA
        std::ofstream file;
        file.open("dataZ.txt");
        file << "%%z\trayleigh_scattering\taerosol_scattering\tozone_absorption\taerosol_absorption\taerosol_density\taerosol_cross_section_scattering\taerosol_cross_section_absorption" << std::endl;
#endif

        for (int i = 0; i < 8699; i++) {
#ifdef PRINT_DATA
            const Wavelength wl = Wavelength(500.f);
            const auto s_m = get_rayleigh_scattering(p, wl);
            const auto s_a = get_aerosol_scattering(p, wl);
            const auto a_m = get_ozone_absorption(p, wl);
            const auto a_a = get_aerosol_absorption(p, wl);
            file << Utils::get_first<Float, ScalarFloat>(get_height(p)) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(s_m[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(s_a[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(a_m[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(a_a[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>((s_m + s_a + a_m + a_a)[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(m_AerosolModel->get_density(get_height(p)) * m_turbidity) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(m_AerosolModel->get_scattering(wl)[0]) << '\t' <<
                    Utils::get_first<Float, ScalarFloat>(m_AerosolModel->get_absorption(wl)[0]) << std::endl;
#endif

            p = p + ScalarVector3f(0, ScalarFloat(0.01) / m_earth_scale, 0);

            extinction = get_max_extinction(p);

            if (enoki::all(extinction > m_max_extinction)) {
                m_max_extinction = extinction;
            }
        }

#ifdef DEBUG_INFO
        Log(Info, "Initialized Earth with \"%s\" radius.", m_earth_radius);
        Log(Info, "Initialized Earth scale as \"%s\".", m_earth_scale);
        Log(Info, "Initialized worldToLocal transform as \"%s\"", worldToLocal);
        std::ostringstream oss;
        oss << m_aabb;
        Log(Info, "Initialized Bounding Box as \"%s\"", oss.str());
#endif
        Log(Info, "Max extinction = \"%s\"", m_max_extinction);

#ifdef PRINT_DATA
        file.close();
        file.open("dataWl.txt");
        file << "%%wl\taerosol_cross_section_scattering\taerosol_cross_section_absorption\taerosol_cross_section_extinction" << std::endl;
        for (int i = 1000; i < 10000; i++) {
            const Wavelength wl = Wavelength(float(i) / 10.f);
            const auto aerosol_cross_section_scattering = Utils::get_first<Float, ScalarFloat>(m_AerosolModel->get_scattering(wl)[0]);
            const auto aerosol_cross_section_absorption = Utils::get_first<Float, ScalarFloat>(m_AerosolModel->get_absorption(wl)[0]);
            file << Utils::get_first<Float, ScalarFloat>(wl[0]) << '\t' <<
                    aerosol_cross_section_scattering << '\t' <<
                    aerosol_cross_section_absorption << '\t' <<
                    aerosol_cross_section_scattering + aerosol_cross_section_absorption << std::endl;
        }
        file.close();

        exit(0);
#endif
    }

    UnpolarizedSpectrum
    get_combined_extinction(const MediumInteraction3f & /* mi */,
                            Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_max_extinction;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

#ifdef DEBUG_INFO
        Log(Info, "World vector p \"%s\"", mi.p);
#endif

        const auto p = worldToLocal.transform_affine(mi.p);
#ifdef DEBUG_INFO
        const Float z = get_height(p);
        const float hmaxz = enoki::hmax(z), hminz = enoki::hmin(z);
        if (hmaxz > max_z)
            max_z = hmaxz;
        if (hminz < min_z)
            min_z = hminz;
#endif
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
        const auto sigman = get_combined_extinction(mi, active) - sigmat;

        return { sigmas, sigman, sigmat };
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        return m_aabb.ray_intersect(ray);
    }

    void traverse(TraversalCallback *callback) override {
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

#ifdef DEBUG_INFO
        const ScalarFloat aerosol_scattering_prob = Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]) / total;
        Log(Info, "rayleigh_scattering = \"%s\"; aerosol_scattering = \"%s\"; rayleigh_scattering_prob = \"%s\"; aerosol_scattering_prob = \"%s\"", Utils::get_first<Float, ScalarFloat>(rayleigh_scattering[0]), Utils::get_first<Float, ScalarFloat>(aerosol_scattering[0]), rayleigh_scattering_prob, aerosol_scattering_prob);
#endif

        // Russian roulette
        if (m_random_generator<ScalarFloat>.next_float32() < rayleigh_scattering_prob)
            return m_phase_function.get(); // P(molecular) = [0, rayleigh_scattering_prob)
        else
            return m_aerosol_phase_function.get(); // P(aerosol) = [rayleigh_scattering_prob,  1)
    }

    /**
     * Computes the geopotential height of a point p given in local coordinates.
     * @param p
     * @return
     */
    Float get_height(const Vector3f &p) const {
        const Float l = norm(p);
        const Float height = (l - m_earth_radius) * m_earth_scale;
#ifdef DEBUG_INFO
        Log(Info, "Vector p \"%s\"", p);
        Log(Info, "Pre-Height of ray collision from the earth center \"%s\".", l);
        Log(Info, "Height of ray collision from the earth center \"%s\" km.", height);
#endif
        return height;
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
        finalResult += m_AerosolModel->get_density_float(h) * m_turbidity * m_AerosolModel->get_absorption();

        // get_aerosol_scattering
        finalResult += m_AerosolModel->get_density_float(h) * m_turbidity * m_AerosolModel->get_scattering();

        return finalResult;
    }

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
        const Float h = get_height(p);

        const Spectrum cross_section = m_rayleighScattering.get_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl);

        const auto density = m_standardAtmosphere.get_number_density<Float, UInt32, Mask>(h);

        const Spectrum finalResult = density * cross_section;

        return finalResult * Float(1e-1f);
    }

    Spectrum get_ozone_absorption(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);

        const Spectrum cross_section = m_rayleighScattering.get_ozone_cross_section<Float, UInt32, Mask, Spectrum, Wavelength>(wl) * Float(1e-10f);

        const auto density = m_standardAtmosphere.get_robson_ozone<Float, UInt32, Mask>(h, m_month);

        return cross_section * density;
    }


    // ----------------------------------------------------------------------------
    // Aerosols scattering and absorption functions
    Spectrum get_aerosol_absorption(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);

        const Spectrum cross_section = m_AerosolModel->get_absorption(wl);

        const Float density = m_AerosolModel->get_density(h) * m_turbidity;

        return cross_section * density;
    }

    Spectrum get_aerosol_scattering(const Vector3f &p, const Wavelength &wl) const {
        const Float h = get_height(p);

        const Spectrum cross_section = m_AerosolModel->get_scattering(wl);

        const Float density = m_AerosolModel->get_density(h) * m_turbidity;

        return cross_section * density;
    }

#ifdef DEBUG_INFO
    ~AtmosphereMedium() {
        Log(Info, "Min z = \"%s\"; Max z = \"%s\"", min_z, max_z);
    }
#endif

    MTS_DECLARE_CLASS()
private:
    ref<PhaseFunction> m_aerosol_phase_function;

    ScalarBoundingBox3f m_aabb;

    ScalarTransform4f worldToLocal;

    // Molecular description
    const StandardAtmosphere::StandardAtmosphere m_standardAtmosphere;
    const RayleighScattering::RayleighScattering m_rayleighScattering;

    // Aerosol description
    std::shared_ptr<GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>> m_AerosolModel{};
    int m_month;

    // Earth description
    ScalarFloat m_earth_radius; // In kilometers
    ScalarFloat m_earth_scale; // From 0 to 1 (where 1 is real radius = m_earth_radius)
    ScalarFloat m_turbidity;

    // Total atmosphere values
    ScalarFloat m_max_extinction;
};

MTS_IMPLEMENT_CLASS_VARIANT(AtmosphereMedium, Medium)
MTS_EXPORT_PLUGIN(AtmosphereMedium, "Atmosphere Medium")
NAMESPACE_END(mitsuba)
