//
// Unit tests for the van der Waals residual model (fugacity::VanDerWaals).
//
// Universal identities, derivatives, wrappers, preconditions, dilute-limit
// behavior, deterministic sampling, and fixed/dynamic equivalence are
// registered through support/eos_test_suite.hpp. What remains here is
// vdW-specific:
//   - the residual Helmholtz energy against an independent long-double
//     reference built directly from the spec formulas
//       a0_ii = 27 (R T_c)^2 / (64 P_c),   b_ii = R T_c / (8 P_c),
//       a_ij  = (1 - k_ij) sqrt(a0_ii a0_jj)   (vdW: T-independent, m_ii = 0),
//       a_m   = sum_ij x_i x_j a_ij,           b_m = sum_i x_i b_ii,
//       a_r   = -R T ln(1 - b_m c) - a_m c     (psi_2 = c for vdW)
//   - the critical-point identities p(T_c, c_c) = P_c and dp/dc = 0 at
//     c_c = 1 / (3 b), which pin down the a0/b parameter construction.
//
#include "support/eos_test_suite.hpp"
#include "support/numeric_checks.hpp"
#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/residual_models/van_der_waals.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

using namespace boost::ut;
using namespace fugacity_test;

namespace {

namespace fug = fugacity;

template<std::size_t N> using Input = typename fug::VanDerWaals<N>::SpeciesInput;

// N2 and CO2 critical data.
constexpr Input<2> n2{.T_c = 126.192, .P_c = 3.3958e6};
constexpr Input<2> co2{.T_c = 304.1282, .P_c = 7.3773e6};

constexpr std::array<Input<2>, 2> binary_inputs{{n2, co2}};
constexpr std::array<Input<1>, 1> unary_inputs{{{.T_c = n2.T_c, .P_c = n2.P_c}}};

// Asymmetric on purpose: only the symmetrized average can affect a_m, and the
// reference below applies each k_ij literally in the full double sum.
constexpr std::array<double, 4> binary_kij{0.0, 0.10, 0.04, 0.0};

// A ConstantCp ideal part so the residual can be paired into a full EoS.
template<std::size_t N> auto make_ideal()
{
    std::array<typename fug::ConstantCp<N>::SpeciesInput, N> in{};
    for (std::size_t i = 0; i < N; ++i) {
        in[i] = {.T_ref = 298.15,
                 .p_ref = 1.0e5,
                 .c_p = 29.1 + (2.0 * static_cast<double>(i)),
                 .h_ref = 1000.0 * static_cast<double>(i),
                 .s_ref = 150.0 + (10.0 * static_cast<double>(i))};
    }
    return fug::ConstantCp<N>(in);
}

template<std::size_t N> auto make_dynamic_ideal()
{
    std::vector<fug::ConstantCp<>::SpeciesInput> inputs(N);
    for (std::size_t i = 0; i < N; ++i) {
        inputs[i] = {.T_ref = 298.15,
                     .p_ref = 1.0e5,
                     .c_p = 29.1 + (2.0 * static_cast<double>(i)),
                     .h_ref = 1000.0 * static_cast<double>(i),
                     .s_ref = 150.0 + (10.0 * static_cast<double>(i))};
    }
    return fug::ConstantCp<>{std::span<const fug::ConstantCp<>::SpeciesInput>{inputs}};
}

auto make_fixed_unary_eos() { return fug::EoS{make_ideal<1>(), fug::VanDerWaals<1>{unary_inputs}}; }

auto make_dynamic_unary_eos()
{
    using Residual = fug::VanDerWaals<std::dynamic_extent>;
    const std::vector<Residual::SpeciesInput> inputs{{.T_c = n2.T_c, .P_c = n2.P_c}};
    return fug::EoS{make_dynamic_ideal<1>(), Residual{std::span<const Residual::SpeciesInput>{inputs}}};
}

auto make_fixed_binary_eos() { return fug::EoS{make_ideal<2>(), fug::VanDerWaals<2>{binary_inputs, binary_kij}}; }

auto make_dynamic_binary_eos()
{
    using Residual = fug::VanDerWaals<std::dynamic_extent>;
    const std::vector<Residual::SpeciesInput> inputs{{.T_c = n2.T_c, .P_c = n2.P_c}, {.T_c = co2.T_c, .P_c = co2.P_c}};
    const std::vector<double> kij(binary_kij.begin(), binary_kij.end());
    return fug::EoS{make_dynamic_ideal<2>(),
                   Residual{std::span<const Residual::SpeciesInput>{inputs}, std::span<const double>{kij}}};
}

std::vector<eos_test_state> unary_contract_states()
{
    return {{.c = 100.0, .x = {1.0}, .T = 300.0, .effective_molar_mass = 0.028, .label = "gas"},
            {.c = 8000.0, .x = {1.0}, .T = 320.0, .effective_molar_mass = 0.028, .label = "dense"}};
}

std::vector<eos_test_state> binary_contract_states()
{
    return {{.c = 150.0, .x = {0.3, 0.7}, .T = 310.0, .effective_molar_mass = 0.036, .label = "gas"},
            {.c = 5000.0, .x = {0.6, 0.4}, .T = 350.0, .effective_molar_mass = 0.033, .label = "dense"}};
}

constexpr eos_valid_domain unary_valid_domain{.c_min = 1e-8,
                                              .c_max = 9000.0,
                                              .T_min = 220.0,
                                              .T_max = 500.0,
                                              .minimum_mole_fraction = 0.01,
                                              .seed = 0xC0FFEE,
                                              .random_samples = 50};

constexpr eos_valid_domain binary_valid_domain{.c_min = 1e-8,
                                               .c_max = 9000.0,
                                               .T_min = 240.0,
                                               .T_max = 500.0,
                                               .minimum_mole_fraction = 0.01,
                                               .seed = 0xC0FFEE,
                                               .random_samples = 50};

// --- Independent long-double reference, straight from the spec -------------
constexpr long double Rld = fug::ideal_gas_constant<long double>;

long double vdw_a0(const long double T_c, const long double P_c)
{
    return 27.0L * (Rld * T_c) * (Rld * T_c) / (64.0L * P_c);
}

long double vdw_b(const long double T_c, const long double P_c) { return Rld * T_c / (8.0L * P_c); }

// Molar residual Helmholtz energy a_r [J/mol] for n species.
template<std::size_t N>
long double ref_helmholtz(const std::array<Input<N>, N>& in, const std::array<double, N * N>& kij, long double c,
                          const std::array<double, N>& x, long double T)
{
    long double am = 0.0L;
    long double bm = 0.0L;
    for (std::size_t i = 0; i < N; ++i) {
        bm += static_cast<long double>(x[i]) * vdw_b(in[i].T_c, in[i].P_c);
        for (std::size_t j = 0; j < N; ++j) {
            const long double aij = (1.0L - static_cast<long double>(kij[(i * N) + j])) *
                                    std::sqrt(vdw_a0(in[i].T_c, in[i].P_c) * vdw_a0(in[j].T_c, in[j].P_c));
            am += static_cast<long double>(x[i]) * static_cast<long double>(x[j]) * aij;
        }
    }
    return (-Rld * T * std::log(1.0L - (bm * c))) - (am * c);
}

} // namespace

int main()
{
    suite<"van_der_waals_unary_contracts"> unary_contracts = [] {
        auto dynamic_eos = make_dynamic_unary_eos();
        const auto fixture = eos_test_fixture{.contribution = dynamic_eos.residual(),
                                              .eos = dynamic_eos,
                                              .states = unary_contract_states(),
                                              .domain = unary_valid_domain};
        register_eos_contract_tests(fixture);
        register_residual_contract_tests(
            fixture, {.dilute_concentration = 1e-8, .dilute_tolerance = {.abs = 1e-6, .rel = 1e-6}});
        register_static_dynamic_equivalence_tests(make_fixed_unary_eos(), make_dynamic_unary_eos(),
                                                  unary_contract_states());
    };

    suite<"van_der_waals"> vdw_suite = [] {
        auto dynamic_eos = make_dynamic_binary_eos();
        const auto fixture = eos_test_fixture{.contribution = dynamic_eos.residual(),
                                              .eos = dynamic_eos,
                                              .states = binary_contract_states(),
                                              .domain = binary_valid_domain};
        register_eos_contract_tests(fixture);
        register_residual_contract_tests(
            fixture, {.dilute_concentration = 1e-8, .dilute_tolerance = {.abs = 1e-6, .rel = 1e-6}});
        register_static_dynamic_equivalence_tests(make_fixed_binary_eos(), make_dynamic_binary_eos(),
                                                  binary_contract_states());

        // ===================================================================
        // Pure species: a_r against the closed-form reference over a sweep of
        // states. Pins down a0_ii, b_ii, and the psi_1/psi_2 assembly.
        // ===================================================================
        "pure species residual helmholtz matches closed form"_test = [] {
            const fug::VanDerWaals<1> model(unary_inputs);
            for (const double c : {1.0, 100.0, 5000.0, 15000.0}) {
                for (const double T : {120.0, 300.0, 500.0}) {
                    const long double ref = ref_helmholtz<1>(unary_inputs, {0.0}, c, {1.0}, T);
                    const std::array<double, 1> x{1.0};
                    check_rel("a_r (pure)", model.calc_helmholtz(c, x.data(), T), static_cast<double>(ref), 1e-12);
                }
            }
        };

        // ===================================================================
        // Binary mixture with asymmetric k_ij: the reference applies k_12 and
        // k_21 literally in the full double sum, so this verifies both the
        // mixing rules and that the implementation's internal symmetrization
        // is exact.
        // ===================================================================
        "binary mixture with asymmetric kij matches closed form"_test = [] {
            const fug::VanDerWaals<2> model(binary_inputs, binary_kij);
            for (const double c : {50.0, 2000.0, 9000.0}) {
                for (const double T : {250.0, 320.0, 450.0}) {
                    for (const std::array<double, 2> x : {std::array{0.3, 0.7}, std::array{0.85, 0.15}}) {
                        const long double ref = ref_helmholtz<2>(binary_inputs, binary_kij, c, x, T);
                        check_rel("a_r (binary)", model.calc_helmholtz(c, x.data(), T), static_cast<double>(ref),
                                  1e-12);
                    }
                }
            }
        };

        // ===================================================================
        // Cross-validation against NIST's teqp library (v0.23.1): golden
        // values of alphar = a_r / (R T) computed with teqp.vdWEOS(Tc, pc)
        // .get_Ar00(T, rho, z) for the same critical data as above. Both
        // libraries use the full CODATA gas constant.
        // ===================================================================
        "alphar matches teqp reference values"_test = [] {
            struct Ref {
                double T;      // [K]
                double c;      // [mol/m^3]
                double alphar; // a_r / (R T) [-]
            };
            const double R = fug::ideal_gas_constant<double>;

            const fug::VanDerWaals<1> pure(unary_inputs); // N2
            constexpr std::array<Ref, 4> pure_refs{{
                {.T = 500.0, .c = 1.0, .alphar = 5.7246696369883375e-06},
                {.T = 300.0, .c = 100.0, .alphar = -0.0016133301849820926},
                {.T = 300.0, .c = 5000.0, .alphar = -0.059582547773369049},
                {.T = 120.0, .c = 15000.0, .alphar = -1.1902211525514714},
            }};
            const std::array<double, 1> x1{1.0};
            for (const Ref& r : pure_refs) {
                check_rel("alphar vs teqp (pure N2)", pure.calc_helmholtz(r.c, x1.data(), r.T) / (R * r.T), r.alphar,
                          1e-13);
            }

            // teqp's vdWEOS mixes with sqrt(a_i a_j) and no k_ij.
            const fug::VanDerWaals<2> binary(binary_inputs); // N2 + CO2, kij = 0
            constexpr std::array<Ref, 3> binary_refs{{
                {.T = 300.0, .c = 100.0, .alphar = -0.0063338824515883491},
                {.T = 250.0, .c = 2000.0, .alphar = -0.16509357279216791},
                {.T = 450.0, .c = 9000.0, .alphar = -0.1648018426592458},
            }};
            const std::array<double, 2> x2{0.4, 0.6};
            for (const Ref& r : binary_refs) {
                check_rel("alphar vs teqp (binary)", binary.calc_helmholtz(r.c, x2.data(), r.T) / (R * r.T), r.alphar,
                          1e-13);
            }
        };

        // ===================================================================
        // Omitting kij must equal passing an all-zero matrix.
        // ===================================================================
        "kij defaults to zero"_test = [] {
            const fug::VanDerWaals<2> defaulted(binary_inputs);
            const fug::VanDerWaals<2> zeros(binary_inputs, std::array<double, 4>{});
            const std::array<double, 2> x{0.4, 0.6};
            for (const double c : {100.0, 4000.0}) {
                check_rel("a_r (kij default)", defaulted.calc_helmholtz(c, x.data(), 300.0),
                          zeros.calc_helmholtz(c, x.data(), 300.0), 1e-15);
            }
        };

        // ===================================================================
        // Dynamic constructor with an *omitted* kij span takes the
        // `kij.empty()` branch in init() (the static default `{}` is an all-zero
        // matrix, which is non-empty, so only the dynamic path reaches it). The
        // result must equal passing an explicit all-zero matrix.
        // ===================================================================
        "dynamic empty kij equals zero matrix"_test = [] {
            using DynInput = fug::VanDerWaals<>::SpeciesInput;
            std::vector<DynInput> in_dyn;
            in_dyn.reserve(binary_inputs.size());
            for (const auto& in : binary_inputs) {
                in_dyn.push_back({.T_c = in.T_c, .P_c = in.P_c});
            }
            const fug::VanDerWaals<> empty_kij{std::span<const DynInput>{in_dyn}}; // kij defaulted to empty span
            const std::vector<double> zeros_dyn(4, 0.0);
            const fug::VanDerWaals<> zero_kij{std::span<const DynInput>{in_dyn}, std::span<const double>{zeros_dyn}};
            const std::array<double, 2> x{0.4, 0.6};
            for (const double c : {100.0, 4000.0}) {
                for (const double T : {250.0, 400.0}) {
                    check_rel("a_r (empty kij == zero matrix)", empty_kij.calc_helmholtz(c, x.data(), T),
                              zero_kij.calc_helmholtz(c, x.data(), T), 1e-15);
                }
            }
        };

#ifndef NDEBUG
        // ===================================================================
        // A non-empty kij whose size is not n*n violates the constructor
        // precondition (FUGACITY_ASSERT), which throws std::logic_error in a
        // debug build.
        // ===================================================================
        "wrong-sized kij throws"_test = [] {
            using DynInput = fug::VanDerWaals<>::SpeciesInput;
            std::vector<DynInput> in_dyn;
            in_dyn.reserve(binary_inputs.size());
            for (const auto& in : binary_inputs) {
                in_dyn.push_back({.T_c = in.T_c, .P_c = in.P_c});
            }
            const std::vector<double> kij_bad(3, 0.0); // size 3, not 2*2 = 4
            expect(throws<std::logic_error>([&] {
                const fug::VanDerWaals<> bad{std::span<const DynInput>{in_dyn}, std::span<const double>{kij_bad}};
                (void)bad;
            }));
        };
#endif

        // ===================================================================
        // Pressure-explicit form of vdW for a pure species:
        //     p = c R T / (1 - b c) - a c^2
        // The framework obtains p from the Helmholtz residual via autodiff, so
        // agreement here verifies the a_r assembly against the textbook EoS.
        // ===================================================================
        "pure species pressure matches pressure-explicit form"_test = [] {
            const fug::EoS eos{make_ideal<1>(), fug::VanDerWaals<1>(unary_inputs)};
            const double R = fug::ideal_gas_constant<double>;
            const double a = static_cast<double>(vdw_a0(n2.T_c, n2.P_c));
            const double b = static_cast<double>(vdw_b(n2.T_c, n2.P_c));
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};
            for (const double c : {1.0, 100.0, 5000.0, 15000.0}) {
                for (const double T : {120.0, 300.0, 500.0}) {
                    const double p_ref = (c * R * T / (1.0 - (b * c))) - (a * c * c);
                    check_rel("p == cRT/(1-bc) - a c^2", fug::calc_pressure(eos, c, xs, T), p_ref, 1e-9);
                }
            }
        };

        // ===================================================================
        // vdW critical-point identities for a pure species: at T = T_c and
        // c_c = 1/(3 b) the model must reproduce p = P_c with dp/dc = 0.
        // This is the strongest check that a0 and b were formed correctly.
        // ===================================================================
        "pure species reproduces its critical point"_test = [] {
            const fug::EoS eos{make_ideal<1>(), fug::VanDerWaals<1>(unary_inputs)};
            const double b = static_cast<double>(vdw_b(n2.T_c, n2.P_c));
            const double c_c = 1.0 / (3.0 * b);
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};
            check_rel("p(T_c, c_c) == P_c", fug::calc_pressure(eos, c_c, xs, n2.T_c), n2.P_c, 1e-9);
            const double dpdc = fug::calc_dp_dc(eos, c_c, xs, n2.T_c);
            // dp/dc vanishes at the critical point; scale by R*T_c since the
            // individual terms it cancels from are O(R*T_c).
            expect(std::abs(dpdc) <= 1e-6 * fug::ideal_gas_constant<double> * n2.T_c)
                << "dp/dc at critical point: " << dpdc;
        };
    };

    return ::boost::ut::cfg<>.run();
}
