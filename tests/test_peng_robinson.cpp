//
// Unit tests for the Peng-Robinson residual model (fugacity::PengRobinson).
//
// Universal identities, derivatives, wrappers, preconditions, dilute-limit
// behavior, deterministic sampling, and fixed/dynamic equivalence are
// registered through support/eos_test_suite.hpp. What remains here is
// PR-specific, checked against an independent long-double reference built
// directly from the spec formulas:
//
//   eta_c   = 1 / (1 + (4 - sqrt8)^(1/3) + (4 + sqrt8)^(1/3))
//   Omega_a = (8 + 40 eta_c) / (49 - 37 eta_c),  Omega_b = eta_c / (3 + eta_c)
//   a0_ii   = Omega_a (R T_c)^2 / P_c,           b_ii = Omega_b R T_c / P_c
//   m_ii    = omega-branched correlation (breakpoint at omega = 0.491)
//   a_ii(T) = a0_ii [1 + m_ii (1 - sqrt(T/T_c))]^2
//   a_ij    = (1 - k_ij) sqrt(a_ii a_jj)   <- always +sqrt: |alpha_i||alpha_j|
//   a_r     = -R T ln(1 - b_m c)
//             - a_m ln((D1 b_m c + 1)/(D2 b_m c + 1)) / (b_m (D1 - D2)),
//   D1,D2   = 1 +- sqrt2.
//
// Notable PR-specific cases:
//   - both branches of the m(omega) correlation,
//   - a mixture state where one species has alpha < 0 (T >> T_c of that
//     species), which distinguishes |alpha_i||alpha_j| from alpha_i alpha_j,
//   - the critical-point identities p(T_c, c_c) = P_c, dp/dc = 0 at
//     c_c = eta_c / b (equivalently Z_c = Omega_b / eta_c ~ 0.3074),
//   - the pressure-explicit PR form p = cRT/(1-bc) - a(T) c^2/(1+2bc-(bc)^2).
//
#include "support/eos_test_suite.hpp"
#include "support/numeric_checks.hpp"
#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/residual_models/peng_robinson.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

using namespace boost::ut;
using namespace fugacity_test;

namespace {

namespace fug = fugacity;

template<std::size_t N> using Input = typename fug::PengRobinson<N>::SpeciesInput;

// CH4 and CO2 critical data; CH4's small T_c makes alpha < 0 reachable.
constexpr Input<2> ch4{.T_c = 190.564, .P_c = 4.5992e6, .omega = 0.011};
constexpr Input<2> co2{.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394};
// Fictitious heavy species exercising the omega > 0.491 branch of m(omega).
constexpr Input<1> heavy{.T_c = 650.0, .P_c = 1.5e6, .omega = 0.60};

constexpr std::array<Input<2>, 2> binary_inputs{{ch4, co2}};
constexpr std::array<Input<1>, 1> unary_inputs{{{.T_c = co2.T_c, .P_c = co2.P_c, .omega = co2.omega}}};
constexpr std::array<Input<1>, 1> heavy_inputs{{heavy}};

// Asymmetric on purpose; the reference applies each k_ij literally.
constexpr std::array<double, 4> binary_kij{0.0, 0.12, 0.06, 0.0};

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

auto make_fixed_unary_eos() { return fug::EoS{make_ideal<1>(), fug::PengRobinson<1>{unary_inputs}}; }

auto make_dynamic_unary_eos()
{
    using Residual = fug::PengRobinson<std::dynamic_extent>;
    const std::vector<Residual::SpeciesInput> inputs{{.T_c = co2.T_c, .P_c = co2.P_c, .omega = co2.omega}};
    return fug::EoS{make_dynamic_ideal<1>(), Residual{std::span<const Residual::SpeciesInput>{inputs}}};
}

auto make_fixed_binary_eos() { return fug::EoS{make_ideal<2>(), fug::PengRobinson<2>{binary_inputs, binary_kij}}; }

auto make_dynamic_binary_eos()
{
    using Residual = fug::PengRobinson<std::dynamic_extent>;
    const std::vector<Residual::SpeciesInput> inputs{
        {.T_c = ch4.T_c, .P_c = ch4.P_c, .omega = ch4.omega},
        {.T_c = co2.T_c, .P_c = co2.P_c, .omega = co2.omega},
    };
    const std::vector<double> kij(binary_kij.begin(), binary_kij.end());
    return fug::EoS{make_dynamic_ideal<2>(),
                   Residual{std::span<const Residual::SpeciesInput>{inputs}, std::span<const double>{kij}}};
}

std::vector<eos_test_state> unary_contract_states()
{
    return {{.c = 100.0, .x = {1.0}, .T = 320.0, .effective_molar_mass = 0.044, .label = "gas"},
            {.c = 8000.0, .x = {1.0}, .T = 340.0, .effective_molar_mass = 0.044, .label = "dense"}};
}

std::vector<eos_test_state> binary_contract_states()
{
    return {{.c = 150.0, .x = {0.3, 0.7}, .T = 310.0, .effective_molar_mass = 0.030, .label = "gas"},
            {.c = 5000.0, .x = {0.6, 0.4}, .T = 350.0, .effective_molar_mass = 0.030, .label = "dense"}};
}

constexpr eos_valid_domain unary_valid_domain{.c_min = 1e-8,
                                              .c_max = 9000.0,
                                              .T_min = 240.0,
                                              .T_max = 600.0,
                                              .minimum_mole_fraction = 0.01,
                                              .seed = 0xC0FFEE,
                                              .random_samples = 50};

constexpr eos_valid_domain binary_valid_domain{.c_min = 1e-8,
                                               .c_max = 9000.0,
                                               .T_min = 240.0,
                                               .T_max = 600.0,
                                               .minimum_mole_fraction = 0.01,
                                               .seed = 0xC0FFEE,
                                               .random_samples = 50};

// --- Independent long-double reference, straight from the spec -------------
constexpr long double Rld = fug::ideal_gas_constant<long double>;

long double pr_eta_c()
{
    const long double s8 = std::sqrt(8.0L);
    return 1.0L / (1.0L + std::cbrt(4.0L - s8) + std::cbrt(4.0L + s8));
}

long double pr_omega_a() { return (8.0L + (40.0L * pr_eta_c())) / (49.0L - (37.0L * pr_eta_c())); }
long double pr_omega_b() { return pr_eta_c() / (3.0L + pr_eta_c()); }

long double pr_a0(const long double T_c, const long double P_c)
{
    return pr_omega_a() * (Rld * T_c) * (Rld * T_c) / P_c;
}

long double pr_b(const long double T_c, const long double P_c) { return pr_omega_b() * Rld * T_c / P_c; }

long double pr_m(const long double omega)
{
    if (omega <= 0.491L) {
        return 0.37464L + (1.54226L * omega) - (0.26992L * omega * omega);
    }
    return 0.379642L + (1.48503L * omega) - (0.164423L * omega * omega) + (0.016666L * omega * omega * omega);
}

// Temperature-dependent pure attractive parameter a_ii(T) >= 0.
template<class In> long double pr_aii(const In& in, const long double T)
{
    const long double alpha = 1.0L + (pr_m(in.omega) * (1.0L - std::sqrt(T / in.T_c)));
    return pr_a0(in.T_c, in.P_c) * alpha * alpha;
}

// Molar residual Helmholtz energy a_r [J/mol]; the cross term takes the
// literal +sqrt(a_ii a_jj), i.e. |alpha_i| |alpha_j|.
template<std::size_t N>
long double ref_helmholtz(const std::array<Input<N>, N>& in, const std::array<double, N * N>& kij, long double c,
                          const std::array<double, N>& x, long double T)
{
    constexpr long double d1 = 1.0L + std::numbers::sqrt2_v<long double>;
    constexpr long double d2 = 1.0L - std::numbers::sqrt2_v<long double>;

    long double am = 0.0L;
    long double bm = 0.0L;
    for (std::size_t i = 0; i < N; ++i) {
        bm += static_cast<long double>(x[i]) * pr_b(in[i].T_c, in[i].P_c);
        for (std::size_t j = 0; j < N; ++j) {
            const long double aij =
                (1.0L - static_cast<long double>(kij[(i * N) + j])) * std::sqrt(pr_aii(in[i], T) * pr_aii(in[j], T));
            am += static_cast<long double>(x[i]) * static_cast<long double>(x[j]) * aij;
        }
    }
    const long double psi1 = -std::log(1.0L - (bm * c));
    const long double psi2 = std::log(((d1 * bm * c) + 1.0L) / ((d2 * bm * c) + 1.0L)) / (bm * (d1 - d2));
    return (Rld * T * psi1) - (am * psi2);
}

} // namespace

int main()
{
    suite<"peng_robinson_unary_contracts"> unary_contracts = [] {
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

    suite<"peng_robinson"> pr_suite = [] {
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
        // Delta_1/Delta_2 are compile-time constants with the PR values.
        // ===================================================================
        "delta constants"_test = [] {
            static_assert(fug::PengRobinson<1>::delta1 == 1.0 + std::numbers::sqrt2);
            static_assert(fug::PengRobinson<1>::delta2 == 1.0 - std::numbers::sqrt2);
            expect(fug::PengRobinson<2>::delta1 > 2.41);
        };

        // ===================================================================
        // Pure species (omega <= 0.491 branch): a_r against the closed-form
        // reference. Pins down eta_c, Omega_a/b, a0, b, m and the psi assembly.
        // ===================================================================
        "pure species residual helmholtz matches closed form"_test = [] {
            const fug::PengRobinson<1> model(unary_inputs);
            for (const double c : {1.0, 100.0, 5000.0, 15000.0}) {
                for (const double T : {220.0, 304.1282, 500.0}) {
                    const long double ref = ref_helmholtz<1>(unary_inputs, {0.0}, c, {1.0}, T);
                    const std::array<double, 1> x{1.0};
                    check_rel("a_r (pure)", model.calc_helmholtz(c, x.data(), T), static_cast<double>(ref), 1e-12);
                }
            }
        };

        // ===================================================================
        // omega > 0.491 must select the quartic m(omega) correlation; the two
        // branches differ by ~1% in m at omega = 0.6, well above the tolerance.
        // ===================================================================
        "high-omega species uses the omega > 0.491 correlation"_test = [] {
            const fug::PengRobinson<1> model(heavy_inputs);
            const std::array<double, 1> x{1.0};
            for (const double T : {400.0, 650.0}) {
                const long double ref = ref_helmholtz<1>(heavy_inputs, {0.0}, 800.0, {1.0}, T);
                check_rel("a_r (heavy)", model.calc_helmholtz(800.0, x.data(), T), static_cast<double>(ref), 1e-12);
            }
        };

        // ===================================================================
        // Binary mixture with asymmetric k_ij against the literal double sum.
        // ===================================================================
        "binary mixture with asymmetric kij matches closed form"_test = [] {
            const fug::PengRobinson<2> model(binary_inputs, binary_kij);
            for (const double c : {50.0, 2000.0, 9000.0}) {
                for (const double T : {230.0, 320.0, 450.0}) {
                    for (const std::array<double, 2> x : {std::array{0.3, 0.7}, std::array{0.85, 0.15}}) {
                        const long double ref = ref_helmholtz<2>(binary_inputs, binary_kij, c, x, T);
                        check_rel("a_r (binary)", model.calc_helmholtz(c, x.data(), T), static_cast<double>(ref),
                                  1e-12);
                    }
                }
            }
        };

        // ===================================================================
        // Negative-alpha regime: at T = 2600 K, CH4 (m ~ 0.392) has
        // alpha = 1 + m(1 - sqrt(T/T_c)) < 0 while CO2's alpha stays positive.
        // The cross term must follow sqrt(a_11 a_22) = |alpha_1||alpha_2| > 0;
        // a signed alpha_1 alpha_2 implementation gets the wrong sign here.
        // ===================================================================
        "mixture cross term uses |alpha| when one alpha is negative"_test = [] {
            const double T = 2600.0;
            // Sanity-check the state really is in the negative-alpha regime.
            const long double alpha_ch4 = 1.0L + (pr_m(ch4.omega) * (1.0L - std::sqrt(T / ch4.T_c)));
            expect(alpha_ch4 < 0.0L) << "test state must have alpha(CH4) < 0";

            const fug::PengRobinson<2> model(binary_inputs, binary_kij);
            for (const double c : {20.0, 500.0, 3000.0}) {
                for (const std::array<double, 2> x : {std::array{0.5, 0.5}, std::array{0.9, 0.1}}) {
                    const long double ref = ref_helmholtz<2>(binary_inputs, binary_kij, c, x, T);
                    check_rel("a_r (alpha < 0)", model.calc_helmholtz(c, x.data(), T), static_cast<double>(ref), 1e-12);
                }
            }
        };

        // ===================================================================
        // Cross-validation against NIST's teqp library (v0.23.1): golden
        // values of alphar = a_r / (R T) computed with
        // teqp.canonical_PR(Tc, pc, acentric, kmat).get_Ar00(T, rho, z) for
        // the same critical data as above. The binary case uses
        // kmat = [[0, 0.09], [0.09, 0]], the symmetrized counterpart of this
        // file's asymmetric binary_kij (0.12, 0.06), so it also cross-checks
        // the internal symmetrization. Both libraries use the full CODATA gas
        // constant.
        // ===================================================================
        "alphar matches teqp reference values"_test = [] {
            struct Ref {
                double T;      // [K]
                double c;      // [mol/m^3]
                double alphar; // a_r / (R T) [-]
            };
            const double R = fug::ideal_gas_constant<double>;

            const fug::PengRobinson<1> pure(unary_inputs); // CO2
            constexpr std::array<Ref, 4> pure_refs{{
                {.T = 500.0, .c = 1.0, .alphar = -3.4438456980806519e-05},
                {.T = 300.0, .c = 100.0, .alphar = -0.013328528270659428},
                {.T = 250.0, .c = 5000.0, .alphar = -0.82160353956496968},
                {.T = 320.0, .c = 15000.0, .alphar = -1.1205590999551662},
            }};
            const std::array<double, 1> x1{1.0};
            for (const Ref& r : pure_refs) {
                check_rel("alphar vs teqp (pure CO2)", pure.calc_helmholtz(r.c, x1.data(), r.T) / (R * r.T), r.alphar,
                          1e-13);
            }

            const fug::PengRobinson<2> binary(binary_inputs, binary_kij); // CH4 + CO2
            constexpr std::array<Ref, 3> binary_refs{{
                {.T = 300.0, .c = 100.0, .alphar = -0.010240637172220948},
                {.T = 250.0, .c = 2000.0, .alphar = -0.27582580656291911},
                {.T = 450.0, .c = 9000.0, .alphar = -0.19037527397951859},
            }};
            const std::array<double, 2> x2{0.3, 0.7};
            for (const Ref& r : binary_refs) {
                check_rel("alphar vs teqp (binary)", binary.calc_helmholtz(r.c, x2.data(), r.T) / (R * r.T), r.alphar,
                          1e-13);
            }
        };

        // ===================================================================
        // Omitting kij must equal passing an all-zero matrix.
        // ===================================================================
        "kij defaults to zero"_test = [] {
            const fug::PengRobinson<2> defaulted(binary_inputs);
            const fug::PengRobinson<2> zeros(binary_inputs, std::array<double, 4>{});
            const std::array<double, 2> x{0.4, 0.6};
            for (const double c : {100.0, 4000.0}) {
                check_rel("a_r (kij default)", defaulted.calc_helmholtz(c, x.data(), 300.0),
                          zeros.calc_helmholtz(c, x.data(), 300.0), 1e-15);
            }
        };

        // ===================================================================
        // Dynamic constructor with an *omitted* kij span takes the
        // `kij.empty()` branch in BaseCubic (the static default `{}` is an
        // all-zero matrix, which is non-empty, so only the dynamic path reaches
        // it). The result must equal passing an explicit all-zero matrix.
        // ===================================================================
        "dynamic empty kij equals zero matrix"_test = [] {
            using DynInput = fug::PengRobinson<>::SpeciesInput;
            std::vector<DynInput> in_dyn;
            in_dyn.reserve(binary_inputs.size());
            for (const auto& in : binary_inputs) {
                in_dyn.push_back({.T_c = in.T_c, .P_c = in.P_c, .omega = in.omega});
            }
            const fug::PengRobinson<> empty_kij{std::span<const DynInput>{in_dyn}}; // kij defaulted to empty span
            const std::vector<double> zeros_dyn(4, 0.0);
            const fug::PengRobinson<> zero_kij{std::span<const DynInput>{in_dyn}, std::span<const double>{zeros_dyn}};
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
        // A non-empty kij whose size is not n*n violates the BaseCubic
        // precondition (FUGACITY_ASSERT), which throws std::logic_error in a
        // debug build.
        // ===================================================================
        "wrong-sized kij throws"_test = [] {
            using DynInput = fug::PengRobinson<>::SpeciesInput;
            std::vector<DynInput> in_dyn;
            in_dyn.reserve(binary_inputs.size());
            for (const auto& in : binary_inputs) {
                in_dyn.push_back({.T_c = in.T_c, .P_c = in.P_c, .omega = in.omega});
            }
            const std::vector<double> kij_bad(3, 0.0); // size 3, not 2*2 = 4
            expect(throws<std::logic_error>([&] {
                const fug::PengRobinson<> bad{std::span<const DynInput>{in_dyn}, std::span<const double>{kij_bad}};
                (void)bad;
            }));
        };
#endif

        // ===================================================================
        // Pressure-explicit form of PR for a pure species:
        //     p = c R T / (1 - b c) - a(T) c^2 / (1 + 2 b c - (b c)^2)
        // The framework obtains p from the Helmholtz residual via autodiff, so
        // agreement here verifies the a_r assembly against the textbook EoS.
        // ===================================================================
        "pure species pressure matches pressure-explicit form"_test = [] {
            const fug::EoS eos{make_ideal<1>(), fug::PengRobinson<1>(unary_inputs)};
            const double R = fug::ideal_gas_constant<double>;
            const double b = static_cast<double>(pr_b(co2.T_c, co2.P_c));
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};
            for (const double c : {1.0, 100.0, 5000.0, 15000.0}) {
                for (const double T : {220.0, 320.0, 500.0}) {
                    const double a_T = static_cast<double>(pr_aii(unary_inputs[0], T));
                    const double bc = b * c;
                    const double p_ref = (c * R * T / (1.0 - bc)) - (a_T * c * c / (1.0 + (2.0 * bc) - (bc * bc)));
                    check_rel("p (pressure-explicit PR)", fug::calc_pressure(eos, c, xs, T), p_ref, 1e-9);
                }
            }
        };

        // ===================================================================
        // PR critical-point identities for a pure species: at T = T_c and
        // c_c = eta_c / b the model must reproduce p = P_c with dp/dc = 0
        // (equivalently Z_c = Omega_b / eta_c ~ 0.3074).
        // ===================================================================
        "pure species reproduces its critical point"_test = [] {
            const fug::EoS eos{make_ideal<1>(), fug::PengRobinson<1>(unary_inputs)};
            const double b = static_cast<double>(pr_b(co2.T_c, co2.P_c));
            const double c_c = static_cast<double>(pr_eta_c()) / b;
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};
            check_rel("p(T_c, c_c) == P_c", fug::calc_pressure(eos, c_c, xs, co2.T_c), co2.P_c, 1e-9);
            const double dpdc = fug::calc_dp_dc(eos, c_c, xs, co2.T_c);
            expect(std::abs(dpdc) <= 1e-6 * fug::ideal_gas_constant<double> * co2.T_c)
                << "dp/dc at critical point: " << dpdc;
        };
    };

    return ::boost::ut::cfg<>.run();
}
