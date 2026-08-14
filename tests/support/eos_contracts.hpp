#pragma once

#include "derivative_oracle.hpp"
#include "eos_test_state.hpp"
#include "property_catalog.hpp"
#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/numbers.hpp"

#include <algorithm>
#include <boost/ut.hpp>
#include <cmath>
#include <format>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace fugacity_test {

inline constexpr tolerance structural_tolerance{.abs = 1e-8, .rel = 1e-9};
inline constexpr tolerance identity_tolerance{.abs = 1e-7, .rel = 2e-8};
inline constexpr tolerance derivative_tolerance{.abs = 2e-5, .rel = 2e-5};

template<class Contribution>
void check_contribution_contracts(const Contribution& contribution, const eos_test_state& state)
{
    using namespace boost::ut;
    std::vector<double> rho(state.x.size());
    std::ranges::transform(state.x, rho.begin(), [&](double xi) { return state.c * xi; });
    std::vector<double> partial(state.x.size());
    contribution.calc_partial_helmholtz(rho.data(), state.T, partial.data());
    const double density = contribution.calc_helmholtz_density(rho.data(), state.T);
    const double partial_sum = std::accumulate(partial.begin(), partial.end(), 0.0);
    const double molar = contribution.calc_helmholtz(state.c, state.x.data(), state.T);

    check_close("Psi == sum(Psi_i)", density, partial_sum, structural_tolerance, state);
    check_close("Psi == c*a", density, state.c * molar, structural_tolerance, state);
    expect(std::isfinite(density) && std::isfinite(molar) &&
           std::ranges::all_of(partial, [](double value) { return std::isfinite(value); }))
        << std::format("contribution outputs must be finite: {}", describe(state));
}

template<class EoSPair> void check_property_wrapper_contracts(const EoSPair& eos, const eos_test_state& state)
{
    namespace fug = fugacity;
    auto x = state.x;
    std::vector<double> pointer_gradient(eos.size());
    std::vector<double> wrapper_gradient(eos.size());

#define FUGACITY_CHECK_WRAPPER(NAME, VALUE, DT, DC, DX)                                                              \
    check_close(#VALUE " pointer/container", fug::VALUE(eos, state.c, x.data(), state.T),                               \
                fug::VALUE(eos, state.c, x, state.T), structural_tolerance, state);                                     \
    check_close(#DT " pointer/container", fug::DT(eos, state.c, x.data(), state.T), fug::DT(eos, state.c, x, state.T),   \
                structural_tolerance, state);                                                                          \
    check_close(#DC " pointer/container", fug::DC(eos, state.c, x.data(), state.T), fug::DC(eos, state.c, x, state.T),   \
                structural_tolerance, state);                                                                          \
    fug::DX(eos, state.c, x.data(), state.T, pointer_gradient.data());                                                  \
    fug::DX(eos, state.c, x, state.T, wrapper_gradient);                                                                \
    for (std::size_t component = 0; component < eos.size(); ++component) {                                             \
        check_close(#DX " pointer/container", pointer_gradient[component], wrapper_gradient[component],                \
                    structural_tolerance, state);                                                                      \
    }
    FUGACITY_TEST_PROPERTY_CATALOG(FUGACITY_CHECK_WRAPPER)
#undef FUGACITY_CHECK_WRAPPER

#define FUGACITY_CHECK_MOLAR_MASS_WRAPPER(NAME, VALUE, DT, DC, DX)                                                   \
    check_close(#VALUE " pointer/container", fug::VALUE(eos, state.c, x.data(), state.T, state.effective_molar_mass),   \
                fug::VALUE(eos, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);         \
    check_close(#DT " pointer/container", fug::DT(eos, state.c, x.data(), state.T, state.effective_molar_mass),         \
                fug::DT(eos, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);            \
    check_close(#DC " pointer/container", fug::DC(eos, state.c, x.data(), state.T, state.effective_molar_mass),         \
                fug::DC(eos, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);            \
    fug::DX(eos, state.c, x.data(), state.T, state.effective_molar_mass, pointer_gradient.data());                      \
    fug::DX(eos, state.c, x, state.T, state.effective_molar_mass, wrapper_gradient);                                    \
    for (std::size_t component = 0; component < eos.size(); ++component) {                                             \
        check_close(#DX " pointer/container", pointer_gradient[component], wrapper_gradient[component],                \
                    structural_tolerance, state);                                                                      \
    }
    FUGACITY_TEST_MOLAR_MASS_PROPERTY_CATALOG(FUGACITY_CHECK_MOLAR_MASS_WRAPPER)
#undef FUGACITY_CHECK_MOLAR_MASS_WRAPPER
}

template<class EoSPair> void check_public_preconditions(const EoSPair& eos, const eos_test_state& state)
{
    using namespace boost::ut;
    namespace fug = fugacity;
    auto x = state.x;
    const double bad_temperature = -1.0;
    auto must_reject_temperature = [&](std::string_view name, auto&& operation) {
        expect(throws<std::domain_error>(operation)) << name << " must reject T <= 0";
    };

#define FUGACITY_CHECK_TEMPERATURE(NAME, VALUE)                                                                      \
    must_reject_temperature(#VALUE, [&] { (void)fug::VALUE(eos, state.c, x, bad_temperature); })
    FUGACITY_CHECK_TEMPERATURE(pressure, calc_pressure);
    FUGACITY_CHECK_TEMPERATURE(internal_energy, calc_internal_energy);
    FUGACITY_CHECK_TEMPERATURE(enthalpy, calc_enthalpy);
    FUGACITY_CHECK_TEMPERATURE(entropy, calc_entropy);
    FUGACITY_CHECK_TEMPERATURE(gibbs, calc_gibbs);
    FUGACITY_CHECK_TEMPERATURE(dp_dc, calc_dp_dc);
    FUGACITY_CHECK_TEMPERATURE(dp_dT, calc_dp_dT);
    FUGACITY_CHECK_TEMPERATURE(cv, calc_cv);
    FUGACITY_CHECK_TEMPERATURE(cp, calc_cp);
#undef FUGACITY_CHECK_TEMPERATURE

    must_reject_temperature("calc_sound_speed_squared", [&] {
        (void)fug::calc_sound_speed_squared(eos, state.c, x, bad_temperature, state.effective_molar_mass);
    });

#ifndef NDEBUG
    std::vector<double> bad_x(eos.size() + 1, 1.0 / static_cast<double>(eos.size() + 1));
    std::vector<double> bad_gradient(eos.size() + 1);
    auto must_reject_size = [&](std::string_view name, auto&& operation) {
        expect(throws<std::logic_error>(operation)) << name << " must reject mismatched composition size";
    };
#define FUGACITY_CHECK_SIZE(NAME, VALUE, DT, DC, DX)                                                                 \
    must_reject_size(#VALUE, [&] { (void)fug::VALUE(eos, state.c, bad_x, state.T); });                                  \
    must_reject_size(#DT, [&] { (void)fug::DT(eos, state.c, bad_x, state.T); });                                        \
    must_reject_size(#DC, [&] { (void)fug::DC(eos, state.c, bad_x, state.T); });                                        \
    must_reject_size(#DX, [&] { fug::DX(eos, state.c, bad_x, state.T, bad_gradient); });
    FUGACITY_TEST_PROPERTY_CATALOG(FUGACITY_CHECK_SIZE)
#undef FUGACITY_CHECK_SIZE
    must_reject_size("calc_sound_speed_squared", [&] {
        (void)fug::calc_sound_speed_squared(eos, state.c, bad_x, state.T, state.effective_molar_mass);
    });
    must_reject_size("calc_sound_speed_squared_dT", [&] {
        (void)fug::calc_sound_speed_squared_dT(eos, state.c, bad_x, state.T, state.effective_molar_mass);
    });
    must_reject_size("calc_sound_speed_squared_dc", [&] {
        (void)fug::calc_sound_speed_squared_dc(eos, state.c, bad_x, state.T, state.effective_molar_mass);
    });
    must_reject_size("calc_sound_speed_squared_dx", [&] {
        fug::calc_sound_speed_squared_dx(eos, state.c, bad_x, state.T, state.effective_molar_mass, bad_gradient);
    });
#endif
}

template<class EoSPair> void check_complete_eos_identities(const EoSPair& eos, const eos_test_state& state)
{
    namespace fug = fugacity;
    auto x = state.x;
    const double a = fug::calc_helmholtz(eos, state.c, x, state.T);
    const double a_sum = eos.ideal().calc_helmholtz(state.c, x.data(), state.T) +
                         eos.residual().calc_helmholtz(state.c, x.data(), state.T);
    const double pressure = fug::calc_pressure(eos, state.c, x, state.T);
    const double internal_energy = fug::calc_internal_energy(eos, state.c, x, state.T);
    const double enthalpy = fug::calc_enthalpy(eos, state.c, x, state.T);
    const double entropy = fug::calc_entropy(eos, state.c, x, state.T);
    const double gibbs = fug::calc_gibbs(eos, state.c, x, state.T);

    check_close("total a == ideal a + residual a", a, a_sum, identity_tolerance, state);
    check_close("u == a + T*s", internal_energy, a + state.T * entropy, identity_tolerance, state);
    check_close("h == u + p/c", enthalpy, internal_energy + pressure / state.c, identity_tolerance, state);
    check_close("g == a + p/c", gibbs, a + pressure / state.c, identity_tolerance, state);
    check_close("g == h - T*s", gibbs, enthalpy - state.T * entropy, identity_tolerance, state);
    check_close("mixed pressure partials", fug::calc_dp_dc_dT(eos, state.c, x, state.T),
                fug::calc_dp_dT_dc(eos, state.c, x, state.T), derivative_tolerance, state);

    std::vector<double> rho(eos.size());
    std::ranges::transform(x, rho.begin(), [&](double xi) { return state.c * xi; });
    std::vector<double> chemical_potential(eos.size());
    fug::calc_chemical_potential(eos, std::span<const double>{rho}, state.T, std::span<double>{chemical_potential});
    const double psi = state.c * a;
    const double euler_pressure = std::inner_product(rho.begin(), rho.end(), chemical_potential.begin(), -psi);
    check_close("Euler pressure", pressure, euler_pressure, identity_tolerance, state);

    std::vector<double> log_phi(eos.size());
    std::vector<double> fugacity(eos.size());
    fug::calc_log_fugacity_coeff(eos, state.c, std::span<const double>{x}, state.T, std::span<const double>{rho},
                                std::span<double>{log_phi});
    fug::calc_fugacity(eos, std::span<const double>{rho}, state.T, std::span<double>{fugacity});
    for (std::size_t component = 0; component < eos.size(); ++component) {
        check_close("fugacity/log(phi) consistency", fugacity[component],
                    x[component] * pressure * std::exp(log_phi[component]), identity_tolerance, state);
    }
}

template<class EoSPair>
void check_public_first_derivatives(const EoSPair& eos, const eos_test_state& state, const eos_valid_domain& domain)
{
    namespace fug = fugacity;
    auto x = state.x;

#define FUGACITY_CHECK_DERIVATIVES(NAME, VALUE, DT, DC, DX)                                                          \
    {                                                                                                                  \
        const auto reference_T =                                                                                       \
            adaptive_first_derivative<double>([&](double value) { return fug::VALUE(eos, state.c, x, value); },         \
                                              state.T, 1.0, domain.T_min, domain.T_max);                               \
        check_close(#DT, fug::DT(eos, state.c, x, state.T), reference_T.value, derivative_tolerance, state,             \
                    reference_T.step, reference_T.error);                                                              \
        const auto reference_c =                                                                                       \
            adaptive_first_derivative<double>([&](double value) { return fug::VALUE(eos, value, x, state.T); },         \
                                              state.c, 1.0, domain.c_min, domain.c_max);                               \
        check_close(#DC, fug::DC(eos, state.c, x, state.T), reference_c.value, derivative_tolerance, state,             \
                    reference_c.step, reference_c.error);                                                              \
        std::vector<double> gradient(eos.size());                                                                      \
        fug::DX(eos, state.c, x, state.T, gradient);                                                                    \
        for (std::size_t component = 0; component + 1 < eos.size(); ++component) {                                     \
            const std::size_t dependent = eos.size() - 1;                                                              \
            const double lower = std::max(domain.minimum_mole_fraction - x[component], x[dependent] - 1.0);            \
            const double upper = std::min(1.0 - x[component], x[dependent] - domain.minimum_mole_fraction);            \
            const auto reference_x = adaptive_first_derivative<double>(                                                \
                [&](double delta) {                                                                                    \
                    auto perturbed = x;                                                                                \
                    perturbed[component] += delta;                                                                     \
                    perturbed[dependent] -= delta;                                                                     \
                    return fug::VALUE(eos, state.c, perturbed, state.T);                                                \
                },                                                                                                     \
                0.0, 1.0, lower, upper);                                                                               \
            check_close(std::format("{}[{}]-{}[{}] (simplex tangent)", #DX, component, #DX, dependent),                \
                        gradient[component] - gradient[dependent], reference_x.value, derivative_tolerance, state,     \
                        reference_x.step, reference_x.error);                                                          \
        }                                                                                                              \
    }
    FUGACITY_TEST_PROPERTY_CATALOG(FUGACITY_CHECK_DERIVATIVES)
#undef FUGACITY_CHECK_DERIVATIVES

#define FUGACITY_CHECK_MOLAR_MASS_DERIVATIVES(NAME, VALUE, DT, DC, DX)                                               \
    {                                                                                                                  \
        const auto reference_T = adaptive_first_derivative<double>(                                                    \
            [&](double value) { return fug::VALUE(eos, state.c, x, value, state.effective_molar_mass); }, state.T, 1.0, \
            domain.T_min, domain.T_max);                                                                               \
        check_close(#DT, fug::DT(eos, state.c, x, state.T, state.effective_molar_mass), reference_T.value,              \
                    derivative_tolerance, state, reference_T.step, reference_T.error);                                 \
        const auto reference_c = adaptive_first_derivative<double>(                                                    \
            [&](double value) { return fug::VALUE(eos, value, x, state.T, state.effective_molar_mass); }, state.c, 1.0, \
            domain.c_min, domain.c_max);                                                                               \
        check_close(#DC, fug::DC(eos, state.c, x, state.T, state.effective_molar_mass), reference_c.value,              \
                    derivative_tolerance, state, reference_c.step, reference_c.error);                                 \
        std::vector<double> gradient(eos.size());                                                                      \
        fug::DX(eos, state.c, x, state.T, state.effective_molar_mass, gradient);                                        \
        for (std::size_t component = 0; component + 1 < eos.size(); ++component) {                                     \
            const std::size_t dependent = eos.size() - 1;                                                              \
            const double lower = std::max(domain.minimum_mole_fraction - x[component], x[dependent] - 1.0);            \
            const double upper = std::min(1.0 - x[component], x[dependent] - domain.minimum_mole_fraction);            \
            const auto reference_x = adaptive_first_derivative<double>(                                                \
                [&](double delta) {                                                                                    \
                    auto perturbed = x;                                                                                \
                    perturbed[component] += delta;                                                                     \
                    perturbed[dependent] -= delta;                                                                     \
                    return fug::VALUE(eos, state.c, perturbed, state.T, state.effective_molar_mass);                    \
                },                                                                                                     \
                0.0, 1.0, lower, upper);                                                                               \
            check_close(std::format("{}[{}]-{}[{}] (simplex tangent)", #DX, component, #DX, dependent),                \
                        gradient[component] - gradient[dependent], reference_x.value, derivative_tolerance, state,     \
                        reference_x.step, reference_x.error);                                                          \
        }                                                                                                              \
    }
    FUGACITY_TEST_MOLAR_MASS_PROPERTY_CATALOG(FUGACITY_CHECK_MOLAR_MASS_DERIVATIVES)
#undef FUGACITY_CHECK_MOLAR_MASS_DERIVATIVES
}

template<class EoSPair>
void check_ideal_gas_contracts(const EoSPair& eos, const eos_test_state& state, const eos_valid_domain& domain)
{
    namespace fug = fugacity;
    const double R = fug::ideal_gas_constant<double>;
    auto x = state.x;
    const double pressure = fug::calc_pressure(eos, state.c, x, state.T);
    check_close("ideal p == cRT", pressure, state.c * R * state.T, identity_tolerance, state);
    check_close("ideal Z == 1", pressure / (state.c * R * state.T), 1.0, identity_tolerance, state);
    check_close("ideal dp/dc == RT", fug::calc_dp_dc(eos, state.c, x, state.T), R * state.T, identity_tolerance, state);
    check_close("ideal dp/dT == cR", fug::calc_dp_dT(eos, state.c, x, state.T), state.c * R, identity_tolerance, state);
    check_close("ideal cp-cv == R", fug::calc_cp(eos, state.c, x, state.T) - fug::calc_cv(eos, state.c, x, state.T), R,
                identity_tolerance, state);

    const double other_c = state.c == domain.c_min ? (state.c + domain.c_max) / 2.0 : (state.c + domain.c_min) / 2.0;
    check_close("ideal h independent of c", fug::calc_enthalpy(eos, state.c, x, state.T),
                fug::calc_enthalpy(eos, other_c, x, state.T), identity_tolerance, state);
    check_close("ideal u independent of c", fug::calc_internal_energy(eos, state.c, x, state.T),
                fug::calc_internal_energy(eos, other_c, x, state.T), identity_tolerance, state);
    check_close("ideal cp independent of c", fug::calc_cp(eos, state.c, x, state.T),
                fug::calc_cp(eos, other_c, x, state.T), identity_tolerance, state);
    check_close("ideal cv independent of c", fug::calc_cv(eos, state.c, x, state.T),
                fug::calc_cv(eos, other_c, x, state.T), identity_tolerance, state);

    std::vector<double> rho(eos.size());
    std::ranges::transform(x, rho.begin(), [&](double xi) { return state.c * xi; });
    std::vector<double> log_phi(eos.size());
    std::vector<double> fugacity(eos.size());
    fug::calc_log_fugacity_coeff(eos, state.c, std::span<const double>{x}, state.T, std::span<const double>{rho},
                                std::span<double>{log_phi});
    fug::calc_fugacity(eos, std::span<const double>{rho}, state.T, std::span<double>{fugacity});
    for (std::size_t component = 0; component < eos.size(); ++component) {
        check_close("ideal log(phi) == 0", log_phi[component], 0.0, identity_tolerance, state);
        check_close("ideal fugacity == rho_i RT", fugacity[component], rho[component] * R * state.T, identity_tolerance,
                    state);
    }
}

struct residual_contract_options {
    double dilute_concentration{};
    tolerance dilute_tolerance{.abs = 1e-5, .rel = 1e-5};
};

template<class EoSPair>
void check_residual_dilute_limit(const EoSPair& eos, eos_test_state state, residual_contract_options options)
{
    namespace fug = fugacity;
    const double R = fug::ideal_gas_constant<double>;
    state.c = options.dilute_concentration;
    state.label += "/dilute-limit";
    auto x = state.x;
    const double residual_helmholtz = eos.residual().calc_helmholtz(state.c, x.data(), state.T);
    const double pressure = fug::calc_pressure(eos, state.c, x, state.T);
    check_close("dilute residual a -> 0", residual_helmholtz, 0.0, options.dilute_tolerance, state);
    check_close("dilute Z -> 1", pressure / (state.c * R * state.T), 1.0, options.dilute_tolerance, state);

    std::vector<double> rho(eos.size());
    std::ranges::transform(x, rho.begin(), [&](double xi) { return state.c * xi; });
    std::vector<double> log_phi(eos.size());
    fug::calc_log_fugacity_coeff(eos, state.c, std::span<const double>{x}, state.T, std::span<const double>{rho},
                                std::span<double>{log_phi});
    for (double value : log_phi) {
        check_close("dilute log(phi) -> 0", value, 0.0, options.dilute_tolerance, state);
    }
}

template<class StaticEoS, class DynamicEoS>
void check_static_dynamic_equivalence(const StaticEoS& fixed, const DynamicEoS& dynamic, const eos_test_state& state)
{
    namespace fug = fugacity;
    if (fixed.size() != dynamic.size()) {
        throw std::invalid_argument("static/dynamic fixtures have different component counts");
    }
    auto x = state.x;
    std::vector<double> rho(fixed.size());
    std::ranges::transform(x, rho.begin(), [&](double xi) { return state.c * xi; });
    auto check_contribution = [&](std::string_view name, const auto& fixed_contribution,
                                  const auto& dynamic_contribution) {
        check_close(std::format("{} molar static/dynamic", name),
                    fixed_contribution.calc_helmholtz(state.c, x.data(), state.T),
                    dynamic_contribution.calc_helmholtz(state.c, x.data(), state.T), structural_tolerance, state);
        check_close(std::format("{} density static/dynamic", name),
                    fixed_contribution.calc_helmholtz_density(rho.data(), state.T),
                    dynamic_contribution.calc_helmholtz_density(rho.data(), state.T), structural_tolerance, state);
        std::vector<double> fixed_partial(fixed.size());
        std::vector<double> dynamic_partial(dynamic.size());
        fixed_contribution.calc_partial_helmholtz(rho.data(), state.T, fixed_partial.data());
        dynamic_contribution.calc_partial_helmholtz(rho.data(), state.T, dynamic_partial.data());
        for (std::size_t component = 0; component < fixed.size(); ++component) {
            check_close(std::format("{} partial[{}] static/dynamic", name, component), fixed_partial[component],
                        dynamic_partial[component], structural_tolerance, state);
        }
    };
    check_contribution("ideal", fixed.ideal(), dynamic.ideal());
    check_contribution("residual", fixed.residual(), dynamic.residual());

    std::vector<double> fixed_gradient(fixed.size());
    std::vector<double> dynamic_gradient(dynamic.size());
#define FUGACITY_CHECK_STATIC_DYNAMIC(NAME, VALUE, DT, DC, DX)                                                       \
    check_close(#VALUE " static/dynamic", fug::VALUE(fixed, state.c, x, state.T),                                       \
                fug::VALUE(dynamic, state.c, x, state.T), structural_tolerance, state);                                 \
    check_close(#DT " static/dynamic", fug::DT(fixed, state.c, x, state.T), fug::DT(dynamic, state.c, x, state.T),       \
                structural_tolerance, state);                                                                          \
    check_close(#DC " static/dynamic", fug::DC(fixed, state.c, x, state.T), fug::DC(dynamic, state.c, x, state.T),       \
                structural_tolerance, state);                                                                          \
    fug::DX(fixed, state.c, x, state.T, fixed_gradient);                                                                \
    fug::DX(dynamic, state.c, x, state.T, dynamic_gradient);                                                            \
    for (std::size_t component = 0; component < fixed.size(); ++component) {                                           \
        check_close(#DX " static/dynamic", fixed_gradient[component], dynamic_gradient[component],                     \
                    structural_tolerance, state);                                                                      \
    }
    FUGACITY_TEST_PROPERTY_CATALOG(FUGACITY_CHECK_STATIC_DYNAMIC)
#undef FUGACITY_CHECK_STATIC_DYNAMIC

#define FUGACITY_CHECK_MOLAR_MASS_STATIC_DYNAMIC(NAME, VALUE, DT, DC, DX)                                            \
    check_close(#VALUE " static/dynamic", fug::VALUE(fixed, state.c, x, state.T, state.effective_molar_mass),           \
                fug::VALUE(dynamic, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);     \
    check_close(#DT " static/dynamic", fug::DT(fixed, state.c, x, state.T, state.effective_molar_mass),                 \
                fug::DT(dynamic, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);        \
    check_close(#DC " static/dynamic", fug::DC(fixed, state.c, x, state.T, state.effective_molar_mass),                 \
                fug::DC(dynamic, state.c, x, state.T, state.effective_molar_mass), structural_tolerance, state);        \
    fug::DX(fixed, state.c, x, state.T, state.effective_molar_mass, fixed_gradient);                                    \
    fug::DX(dynamic, state.c, x, state.T, state.effective_molar_mass, dynamic_gradient);                                \
    for (std::size_t component = 0; component < fixed.size(); ++component) {                                           \
        check_close(#DX " static/dynamic", fixed_gradient[component], dynamic_gradient[component],                     \
                    structural_tolerance, state);                                                                      \
    }
    FUGACITY_TEST_MOLAR_MASS_PROPERTY_CATALOG(FUGACITY_CHECK_MOLAR_MASS_STATIC_DYNAMIC)
#undef FUGACITY_CHECK_MOLAR_MASS_STATIC_DYNAMIC
}

} // namespace fugacity_test
