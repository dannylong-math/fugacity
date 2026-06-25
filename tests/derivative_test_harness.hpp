#pragma once
//
// Reusable, templated derivative-consistency harness.
//
// Given *any* EoS pair (an instance of synthesize::EoS<Ideal, Residual>) and a
// thermodynamic state, this checks every derivative-based quantity in
// core_calculations.hpp against an independent finite-difference reference.
//
// The references are built from textbook thermodynamic identities applied to
// the model's own fundamental Helmholtz functions, evaluated in `long double`
// with a 4th-order central finite-difference stencil:
//
//     f'(x) = (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h)
//
// Because the references do NOT reuse the library's internal reduced-derivative
// (alpha / lambda) machinery, a passing run validates both the Enzyme autodiff
// AND the closed-form thermodynamic formulas wired on top of it.
//
// To add derivative tests for a newly implemented EoS, simply call
// `run_derivative_consistency_tests(eos, c, x, T)` from a test case.
//
// Units (see eos_test_models.hpp / the library headers):
//   c   [mol/m^3], x [-], T [K], pressure [Pa], energies [J/mol],
//   entropy / heat capacities [J/mol/K], chemical potential [J/mol].
//
#include "synthesize/core/core_calculations.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <format>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace synthesize_test {

// 4th-order central first derivative of a unary callable, in long double.
[[nodiscard]] long double central_diff(auto f, long double x, long double h)
{
    return (-f(x + (2 * h)) + (8 * f(x + h)) - (8 * f(x - h)) + f(x - (2 * h))) / (12 * h);
}

// Relative-error expectation with an informative message on failure.
inline void check_rel(std::string_view name, double actual, double expected, double reltol)
{
    using namespace boost::ut;
    const double denom = std::max(1.0, std::abs(expected));
    const double rel = std::abs(actual - expected) / denom;
    expect(rel <= reltol) << std::format("{}: actual={:.12g} expected={:.12g} rel_err={:.3e} (tol={:.1e})", name,
                                         actual, expected, rel, reltol);
}

// `rtol` is intentionally looser than bitwise: the pointer core and the
// container wrapper evaluate identical source, but the compiler may contract
// floating-point (FMA) differently across the two inlined call sites, so the
// results can disagree in the last digit or two. 1e-9 still catches any real
// delegation/logic mistake while tolerating that last-bit rounding.
template<std::size_t N, class Model, std::floating_point Number = double>
void run_free_function_consistency_tests(const Model& model, Number rtol = 1e-9)
{
    using namespace boost::ut;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<Number> T_dist(1., 5000.);
    std::uniform_real_distribution<Number> rho_dist(0.0001, 100.);

    std::vector<Number> T_vals(50);
    std::vector<std::array<Number, N>> rho_vecs(50);
    std::ranges::generate(T_vals, [&]() { return T_dist(gen); });
    std::ranges::generate(rho_vecs, [&]() {
        std::array<Number, N> rho;
        for (std::size_t i = 0; i < rho.size(); ++i) {
            rho.at(i) = rho_dist(gen);
        }
        return rho;
    });
    "Temperature values"_test = [&](const Number T) {
        "Rho vec"_test = [&](const auto rho_i) {
            const Number c = std::reduce(rho_i.begin(), rho_i.end(), Number{0});
            std::array<Number, N> x{};
            for (std::size_t i = 0; i < N; ++i) {
                x.at(i) = rho_i.at(i) / c;
            }
            // Scratch gradient buffers for the reverse-mode `_dx` checks. Each
            // `calc_*_dx` zeroes its output first, so reuse across calls is safe.
            std::array<Number, N> grad_ptr{};
            std::array<Number, N> grad_vec{};

            // Every scalar free function has a pointer-based core and a
            // container-based wrapper that must agree exactly. `consistency`
            // invokes the same overload set with `x.data()` (pointer) and `x`
            // (container) and checks the two results match. `std::forward`
            // preserves the value category so the prvalue pointer selects the
            // pointer overload while the lvalue container selects the wrapper.
            auto consistency = [&](std::string_view name, auto fn) { check_rel(name, fn(x.data()), fn(x), rtol); };
            // Same idea for the reverse-mode `_dx` functions, which fill a
            // gradient buffer instead of returning a value.
            auto consistency_grad = [&](std::string_view name, auto ptr_fill, auto vec_fill) {
                ptr_fill(grad_ptr);
                vec_fill(grad_vec);
                for (std::size_t i = 0; i < N; ++i) {
                    check_rel(name, grad_ptr.at(i), grad_vec.at(i), rtol);
                }
            };

#define SYNTHESIZE_CHECK_SCALAR(FN)                                                                                    \
    consistency(#FN, [&](auto&& xx) { return FN(model, c, std::forward<decltype(xx)>(xx), T); })
#define SYNTHESIZE_CHECK_GRAD(FN)                                                                                      \
    consistency_grad(                                                                                                  \
        #FN, [&](auto& g) { FN(model, c, x.data(), T, g.data()); }, [&](auto& g) { FN(model, c, x, T, g); })

            SYNTHESIZE_CHECK_SCALAR(calc_helmholtz);
            SYNTHESIZE_CHECK_SCALAR(calc_helmholtz_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_helmholtz_dc);
            SYNTHESIZE_CHECK_GRAD(calc_helmholtz_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_pressure);
            SYNTHESIZE_CHECK_SCALAR(calc_pressure_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_pressure_dc);
            SYNTHESIZE_CHECK_GRAD(calc_pressure_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_internal_energy);
            SYNTHESIZE_CHECK_SCALAR(calc_internal_energy_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_internal_energy_dc);
            SYNTHESIZE_CHECK_GRAD(calc_internal_energy_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_enthalpy);
            SYNTHESIZE_CHECK_SCALAR(calc_enthalpy_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_enthalpy_dc);
            SYNTHESIZE_CHECK_GRAD(calc_enthalpy_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_entropy);
            SYNTHESIZE_CHECK_SCALAR(calc_entropy_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_entropy_dc);
            SYNTHESIZE_CHECK_GRAD(calc_entropy_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_gibbs);
            SYNTHESIZE_CHECK_SCALAR(calc_gibbs_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_gibbs_dc);
            SYNTHESIZE_CHECK_GRAD(calc_gibbs_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_dp_dc);
            SYNTHESIZE_CHECK_SCALAR(calc_dp_dc_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_dp_dc_dc);
            SYNTHESIZE_CHECK_GRAD(calc_dp_dc_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_dp_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_dp_dT_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_dp_dT_dc);
            SYNTHESIZE_CHECK_GRAD(calc_dp_dT_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_cv);
            SYNTHESIZE_CHECK_SCALAR(calc_cv_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_cv_dc);
            SYNTHESIZE_CHECK_GRAD(calc_cv_dx);

            SYNTHESIZE_CHECK_SCALAR(calc_cp);
            SYNTHESIZE_CHECK_SCALAR(calc_cp_dT);
            SYNTHESIZE_CHECK_SCALAR(calc_cp_dc);
            SYNTHESIZE_CHECK_GRAD(calc_cp_dx);

#undef SYNTHESIZE_CHECK_SCALAR
#undef SYNTHESIZE_CHECK_GRAD

            // Speed-of-sound family carries an extra effective-molar-mass
            // argument, so it does not fit the macros above.
            const Number molar_mass = Number{0.02862};
            consistency("calc_sound_speed_squared", [&](auto&& xx) {
                return calc_sound_speed_squared(model, c, std::forward<decltype(xx)>(xx), T, molar_mass);
            });
            consistency("calc_sound_speed_squared_dT", [&](auto&& xx) {
                return calc_sound_speed_squared_dT(model, c, std::forward<decltype(xx)>(xx), T, molar_mass);
            });
            consistency("calc_sound_speed_squared_dc", [&](auto&& xx) {
                return calc_sound_speed_squared_dc(model, c, std::forward<decltype(xx)>(xx), T, molar_mass);
            });
            consistency_grad(
                "calc_sound_speed_squared_dx",
                [&](auto& g) { calc_sound_speed_squared_dx(model, c, x.data(), T, molar_mass, g.data()); },
                [&](auto& g) { calc_sound_speed_squared_dx(model, c, x, T, molar_mass, g); });
        } | rho_vecs;
    } | T_vals;
}

// ===========================================================================
// Precondition (input-validation) checks for the free functions.
//
// Every container-based free function validates that the mole-fraction container
// matches the component count, and every temperature-dependent function rejects
// a non-positive temperature. This helper drives both failure paths generically
// for the whole property catalogue, so a single call exercises:
//
//   * size mismatch  -> std::logic_error  via SYNTHESIZE_ASSERT. This is a
//     debug-only check (elided under NDEBUG), so those expectations are compiled
//     only when NDEBUG is not defined.
//   * non-positive T -> std::domain_error. This is an always-on check, exercised
//     on both the container wrapper and the pointer-core overload of each
//     temperature-dependent function.
//
//   eos : an EoS<Ideal, Residual>
//   c   : a valid molar concentration           [mol/m^3]
//   x   : valid mole fractions (size N)          [-]
//   T   : a valid (positive) temperature         [K]
// ===========================================================================
template<std::size_t N, class EoSPair>
void run_precondition_tests(const EoSPair& eos, double c, std::array<double, N> x, [[maybe_unused]] double T,
                            double effective_molar_mass = 0.02)
{
    using namespace boost::ut;
    namespace ge = synthesize;

    const double bad_T = -1.0; // non-positive temperature
    const double mm = effective_molar_mass;
    std::array<double, N> xv = x; // valid-size mole fractions
    std::array<double, N> grad{}; // valid-size gradient scratch

    // ---- non-positive temperature: always-on std::domain_error --------------
    auto temp_throws = [&](std::string_view name, auto fn) {
        expect(throws<std::domain_error>([&] { fn(); })) << name << "must reject T <= 0";
    };
    // Each temperature-checked property validates T in both its container wrapper
    // and its pointer-core overload; drive both.
#define SYNTHESIZE_CHECK_T(FN)                                                                                         \
    temp_throws(#FN " (wrapper)", [&] { (void)FN(eos, c, xv, bad_T); });                                               \
    temp_throws(#FN " (core)", [&] { (void)FN(eos, c, xv.data(), bad_T); })

    SYNTHESIZE_CHECK_T(ge::calc_pressure);
    SYNTHESIZE_CHECK_T(ge::calc_internal_energy);
    SYNTHESIZE_CHECK_T(ge::calc_enthalpy);
    SYNTHESIZE_CHECK_T(ge::calc_entropy);
    SYNTHESIZE_CHECK_T(ge::calc_gibbs);
    SYNTHESIZE_CHECK_T(ge::calc_dp_dc);
    SYNTHESIZE_CHECK_T(ge::calc_dp_dT);
    SYNTHESIZE_CHECK_T(ge::calc_cv);
    SYNTHESIZE_CHECK_T(ge::calc_cp);
#undef SYNTHESIZE_CHECK_T

    // Speed of sound carries an extra molar-mass argument.
    temp_throws("calc_sound_speed_squared (wrapper)",
                [&] { (void)ge::calc_sound_speed_squared(eos, c, xv, bad_T, mm); });
    temp_throws("calc_sound_speed_squared (core)",
                [&] { (void)ge::calc_sound_speed_squared(eos, c, xv.data(), bad_T, mm); });

    // Pressure's first-derivative cores validate T as well.
    temp_throws("calc_pressure_dT (core)", [&] { (void)ge::calc_pressure_dT(eos, c, xv.data(), bad_T); });
    temp_throws("calc_pressure_dc (core)", [&] { (void)ge::calc_pressure_dc(eos, c, xv.data(), bad_T); });
    temp_throws("calc_pressure_dx (core)", [&] { ge::calc_pressure_dx(eos, c, xv.data(), bad_T, grad.data()); });

    // ---- mismatched mole-fraction size: std::logic_error (debug only) -------
#ifndef NDEBUG
    // Deliberately wrong-sized (N+1) containers for the size-mismatch checks.
    // Using a fixed-size std::array (not a runtime-sized vector) keeps the
    // `x.size() == eos.size()` comparison a compile-time-false constant, so the
    // failing call does not introduce a spurious runtime branch whose passing
    // side would never be taken.
    std::array<double, N + 1> x_bad{};
    std::array<double, N + 1> grad_bad{};
    auto size_throws = [&](std::string_view name, auto fn) {
        expect(throws<std::logic_error>([&] { fn(); })) << name << "must reject mismatched x size";
    };
#define SYNTHESIZE_CHECK_SIZE_SCALAR(FN) size_throws(#FN, [&] { (void)FN(eos, c, x_bad, T); })
#define SYNTHESIZE_CHECK_SIZE_GRAD(FN) size_throws(#FN, [&] { FN(eos, c, x_bad, T, grad_bad); })

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_helmholtz);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_helmholtz_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_helmholtz_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_helmholtz_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_pressure);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_pressure_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_pressure_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_pressure_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_internal_energy);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_internal_energy_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_internal_energy_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_internal_energy_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_enthalpy);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_enthalpy_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_enthalpy_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_enthalpy_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_entropy);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_entropy_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_entropy_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_entropy_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_gibbs);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_gibbs_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_gibbs_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_gibbs_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dc);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dc_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dc_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_dp_dc_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dT_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_dp_dT_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_dp_dT_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cv);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cv_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cv_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_cv_dx);

    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cp);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cp_dT);
    SYNTHESIZE_CHECK_SIZE_SCALAR(ge::calc_cp_dc);
    SYNTHESIZE_CHECK_SIZE_GRAD(ge::calc_cp_dx);
#undef SYNTHESIZE_CHECK_SIZE_SCALAR
#undef SYNTHESIZE_CHECK_SIZE_GRAD

    // Speed of sound carries an extra molar-mass argument.
    size_throws("calc_sound_speed_squared", [&] { (void)ge::calc_sound_speed_squared(eos, c, x_bad, T, mm); });
    size_throws("calc_sound_speed_squared_dT", [&] { (void)ge::calc_sound_speed_squared_dT(eos, c, x_bad, T, mm); });
    size_throws("calc_sound_speed_squared_dc", [&] { (void)ge::calc_sound_speed_squared_dc(eos, c, x_bad, T, mm); });
    size_throws("calc_sound_speed_squared_dx",
                [&] { ge::calc_sound_speed_squared_dx(eos, c, x_bad, T, mm, grad_bad); });
#endif
}

// ===========================================================================
// Reusable single-call checks usable by ANY model's test suite. Each is a
// templated one-liner; loop over a few states in the caller for coverage.
// ===========================================================================

// Structural Helmholtz consistency for a single contribution model (ideal or
// residual). Applies to EVERY model: it only uses the EquationOfState interface.
//   (1) Psi(rho,T) == sum_i partial_i(rho,T)              (decomposition)
//   (2) Psi(rho,T) == c * a(c,x,T),  c = sum rho, x_i = rho_i/c   (density<->molar)
// No physical reference values are needed; this catches density/molar and
// per-component bookkeeping mistakes regardless of the model's physics.
template<std::size_t N, class Model>
void check_helmholtz_consistency(const Model& model, std::array<double, N> rho, double T, double rtol = 1e-11)
{
    double c = 0.0;
    for (const double r : rho) {
        c += r;
    }
    std::array<double, N> x{};
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = rho[i] / c;
    }

    const double psi = model.calc_helmholtz_density(rho.data(), T);

    std::array<double, N> partial{};
    model.calc_partial_helmholtz(rho.data(), T, partial.data());
    double psi_sum = 0.0;
    for (const double v : partial) {
        psi_sum += v;
    }
    check_rel("Psi == sum_i Psi_i", psi, psi_sum, rtol);
    check_rel("Psi == c * a", psi, c * model.calc_helmholtz(c, x.data(), T), rtol);
}

// Ideal-gas pressure: any ideal model paired with a zero residual must give the
// ideal-gas law p = c R T. Pass the EoS pair (ideal model + NoResidual).
template<std::size_t N, class EoSPair>
void check_ideal_gas_pressure(const EoSPair& eos, double c, std::array<double, N> x, double T, double rtol = 1e-12)
{
    const double R = synthesize::ideal_gas_constant<double>;
    std::array<double, N> xarr = x;
    const std::span<const double, N> xs{xarr};
    check_rel("p == c R T", synthesize::calc_pressure(eos, c, xs, T), c * R * T, rtol);
}

// Euler relation (a thermodynamic identity holding for EVERY EoS):
//     p = sum_i rho_i mu_i - Psi
// with mu_i the chemical potentials (reverse-mode autodiff of the total Psi).
// Cross-checks the chemical potentials against calc_pressure; for mixtures this
// is the guard that the partial-molar autodiff is self-consistent with pressure.
template<std::size_t N, class EoSPair>
void check_euler_pressure(const EoSPair& eos, std::array<double, N> rho, double T, double rtol = 1e-8)
{
    namespace ge = synthesize;
    double c = 0.0;
    for (const double r : rho) {
        c += r;
    }
    std::array<double, N> x{};
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = rho[i] / c;
    }

    std::array<double, N> rhoarr = rho;
    const std::span<const double, N> rhos{rhoarr};
    std::array<double, N> mu{};
    ge::calc_chemical_potential(eos, rhos, T, std::span<double, N>{mu});

    const double psi =
        eos.ideal().calc_helmholtz_density(rho.data(), T) + eos.residual().calc_helmholtz_density(rho.data(), T);
    double p_euler = -psi;
    for (std::size_t i = 0; i < N; ++i) {
        p_euler += rho[i] * mu[i];
    }

    const std::span<const double, N> xs{x};
    check_rel("Euler p == calc_pressure", p_euler, ge::calc_pressure(eos, c, xs, T), rtol);
}

// ---------------------------------------------------------------------------
// Main harness.
//   eos : an EoS<Ideal, Residual>
//   c   : molar concentration                 [mol/m^3]
//   x   : mole fractions (must sum to 1)       [-]
//   T   : temperature                          [K]
//   effective_molar_mass : used only for the speed-of-sound check [kg/mol]
// ---------------------------------------------------------------------------
template<std::size_t N, class EoSPair>
void run_derivative_consistency_tests(const EoSPair& eos, double c, std::array<double, N> x, double T,
                                      double effective_molar_mass = 0.02)
{
    namespace ge = synthesize;
    using ld = long double;

    const ld R = ge::ideal_gas_constant<ld>;
    const ld c0 = c;
    const ld T0 = T;
    std::array<ld, N> xl{};
    for (std::size_t i = 0; i < N; ++i) {
        xl[i] = static_cast<ld>(x[i]);
    }

    const auto& ideal = eos.ideal();
    const auto& residual = eos.residual();

    // --- fundamental Helmholtz functions in long double -------------------
    auto a_res = [&](ld cc, ld TT) { return residual.calc_helmholtz(cc, xl.data(), TT); };
    auto a_tot = [&](ld cc, ld TT) {
        return ideal.calc_helmholtz(cc, xl.data(), TT) + residual.calc_helmholtz(cc, xl.data(), TT);
    };

    const ld hc = c0 * 1e-3L;
    const ld hT = T0 * 1e-3L;

    // residual c-derivatives (at fixed T0)
    const ld a_res_c = central_diff([&](ld cc) { return a_res(cc, T0); }, c0, hc);
    const ld a_res_cc =
        central_diff([&](ld cc) { return central_diff([&](ld c2) { return a_res(c2, T0); }, cc, hc); }, c0, hc);
    // residual mixed c-T derivative: d/dT ( da_res/dc )
    const ld a_res_cT =
        central_diff([&](ld TT) { return central_diff([&](ld cc) { return a_res(cc, TT); }, c0, hc); }, T0, hT);

    // total T-derivatives (at fixed c0)
    const ld a_tot_0 = a_tot(c0, T0);
    const ld a_tot_T = central_diff([&](ld TT) { return a_tot(c0, TT); }, T0, hT);
    const ld a_tot_TT =
        central_diff([&](ld TT) { return central_diff([&](ld T2) { return a_tot(c0, T2); }, TT, hT); }, T0, hT);

    // --- reference thermodynamic quantities -------------------------------
    const ld p_ref = (c0 * R * T0) + (c0 * c0 * a_res_c);
    const ld dp_dc_ref = (R * T0) + (2 * c0 * a_res_c) + (c0 * c0 * a_res_cc);
    const ld dp_dT_ref = (c0 * R) + (c0 * c0 * a_res_cT);
    const ld u_ref = a_tot_0 - (T0 * a_tot_T);
    const ld s_ref = -a_tot_T;
    const ld cv_ref = -T0 * a_tot_TT;
    const ld h_ref = u_ref + (p_ref / c0);
    const ld g_ref = a_tot_0 + (p_ref / c0);
    const ld cp_ref = cv_ref + (T0 * dp_dT_ref * dp_dT_ref / (c0 * c0 * dp_dc_ref));
    const ld w2_ref = cp_ref * dp_dc_ref / (static_cast<ld>(effective_molar_mass) * cv_ref);

    // --- framework (double, Enzyme) values --------------------------------
    std::array<double, N> xarr = x;
    std::span<const double, N> xs{xarr};

    check_rel("helmholtz", ge::calc_helmholtz(eos, c, xs, T), static_cast<double>(a_tot_0), 1e-12);
    check_rel("pressure", ge::calc_pressure(eos, c, xs, T), static_cast<double>(p_ref), 1e-8);
    check_rel("internal_energy", ge::calc_internal_energy(eos, c, xs, T), static_cast<double>(u_ref), 1e-8);
    check_rel("enthalpy", ge::calc_enthalpy(eos, c, xs, T), static_cast<double>(h_ref), 1e-8);
    check_rel("entropy", ge::calc_entropy(eos, c, xs, T), static_cast<double>(s_ref), 1e-8);
    check_rel("gibbs", ge::calc_gibbs(eos, c, xs, T), static_cast<double>(g_ref), 1e-8);
    check_rel("dp_dc", ge::calc_dp_dc(eos, c, xs, T), static_cast<double>(dp_dc_ref), 1e-6);
    check_rel("dp_dT", ge::calc_dp_dT(eos, c, xs, T), static_cast<double>(dp_dT_ref), 1e-6);
    check_rel("cv", ge::calc_cv(eos, c, xs, T), static_cast<double>(cv_ref), 1e-6);
    check_rel("cp", ge::calc_cp(eos, c, xs, T), static_cast<double>(cp_ref), 1e-6);
    check_rel("sound_speed_squared", ge::calc_sound_speed_squared(eos, c, xs, T, effective_molar_mass),
              static_cast<double>(w2_ref), 1e-6);

    // --- partial-molar quantities (reverse-mode autodiff) -----------------
    // State in partial concentrations: rho_i = x_i * c.
    std::array<ld, N> rhol{};
    std::array<double, N> rho{};
    for (std::size_t i = 0; i < N; ++i) {
        rhol[i] = xl[i] * c0;
        rho[i] = x[i] * c;
    }

    auto Psi_tot = [&](std::array<ld, N> r) {
        std::array<ld, N> oi{};
        std::array<ld, N> orr{};
        ideal.calc_partial_helmholtz(r.data(), T0, oi.data());
        residual.calc_partial_helmholtz(r.data(), T0, orr.data());
        ld s{0};
        for (std::size_t i = 0; i < N; ++i) {
            s += oi[i] + orr[i];
        }
        return s;
    };
    auto Psi_res = [&](std::array<ld, N> r) {
        std::array<ld, N> orr{};
        residual.calc_partial_helmholtz(r.data(), T0, orr.data());
        ld s{0};
        for (std::size_t i = 0; i < N; ++i) {
            s += orr[i];
        }
        return s;
    };

    // mu_i = d Psi_tot / d rho_i ; mu_i^res = d Psi_res / d rho_i
    std::array<ld, N> mu_ref{};
    std::array<ld, N> mu_res_ref{};
    for (std::size_t i = 0; i < N; ++i) {
        const ld hr = rhol[i] * 1e-3L;
        mu_ref[i] = central_diff(
            [&](ld v) {
                auto r = rhol;
                r[i] = v;
                return Psi_tot(r);
            },
            rhol[i], hr);
        mu_res_ref[i] = central_diff(
            [&](ld v) {
                auto r = rhol;
                r[i] = v;
                return Psi_res(r);
            },
            rhol[i], hr);
    }

    std::span<const double, N> rhos{rho};
    std::array<double, N> chem{};
    ge::calc_chemical_potential(eos, rhos, T, std::span<double, N>{chem});
    for (std::size_t i = 0; i < N; ++i) {
        check_rel(std::format("chemical_potential[{}]", i), chem[i], static_cast<double>(mu_ref[i]), 1e-7);
    }

    // Z = p / (c R T)
    const ld Z_ref = p_ref / (c0 * R * T0);
    std::array<double, N> logphi{};
    ge::calc_log_fugacity_coeff(eos, c, xs, T, rhos, std::span<double, N>{logphi});
    for (std::size_t i = 0; i < N; ++i) {
        const ld ref = (mu_res_ref[i] / (R * T0)) - std::log(Z_ref);
        check_rel(std::format("log_fugacity_coeff[{}]", i), logphi[i], static_cast<double>(ref), 1e-7);
    }

    std::array<double, N> fug{};
    ge::calc_fugacity(eos, rhos, T, std::span<double, N>{fug});
    for (std::size_t i = 0; i < N; ++i) {
        const ld ref = rhol[i] * R * T0 * std::exp(mu_res_ref[i] / (R * T0));
        check_rel(std::format("fugacity[{}]", i), fug[i], static_cast<double>(ref), 1e-7);
    }
}

} // namespace synthesize_test
