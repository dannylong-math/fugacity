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
#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/residual_models/peng_robinson.hpp"
#include "support/eos_test_state.hpp"
#include "support/eos_test_suite.hpp"
#include "support/numeric_checks.hpp"

#include <array>
#include <bit>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

using namespace boost::ut;
using namespace fugacity_test;

namespace {

namespace fug = fugacity;

template<std::size_t N> using Input = fug::PengRobinson<N>::SpeciesInput;

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

// Deterministic, deliberately nonuniform mixtures used to characterize the
// runtime-size path around the row-blocking threshold and its tails. The
// asymmetry of kij and unique values at every species index make first/last
// and block-boundary indexing errors observable.
constexpr std::uint64_t blocking_characterization_seed = 0xB10C5EED;

template<std::size_t N> struct blocking_case {
    struct species_data {
        double T_c;
        double P_c;
        double omega;
        double c_p;
    };

    std::array<species_data, N> species{};
    std::array<double, N * N> kij{};
    std::array<double, N> x{};
    std::array<double, N> rho{};
    double c{};
    double T{};
};

template<std::size_t N> blocking_case<N> make_blocking_case()
{
    // A reproducible characterization sequence is intentional, not security-sensitive randomness.
    // NOLINTNEXTLINE(bugprone-random-generator-seed)
    std::mt19937_64 generator{blocking_characterization_seed};
    std::uniform_real_distribution<double> unit{0.0, 1.0};
    blocking_case<N> test_case;
    double x_sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        test_case.species[i] = {.T_c = 180.0 + (420.0 * unit(generator)) + static_cast<double>(i),
                                .P_c = 3.5e6 + (4.0e6 * unit(generator)) + (1000.0 * static_cast<double>(i)),
                                .omega = 0.02 + (0.65 * unit(generator)),
                                .c_p = 24.0 + (18.0 * unit(generator)) + (0.125 * static_cast<double>(i))};
        test_case.x[i] = 0.1 + unit(generator) + (0.01 * static_cast<double>(i));
        x_sum += test_case.x[i];
    }
    for (double& xi : test_case.x) {
        xi /= x_sum;
    }
    test_case.c = 1250.0;
    test_case.T = 415.0;
    for (std::size_t i = 0; i < N; ++i) {
        test_case.rho[i] = test_case.c * test_case.x[i];
        for (std::size_t j = 0; j < N; ++j) {
            test_case.kij[(i * N) + j] =
                (i == j) ? 0.0 : -0.08 + (0.20 * unit(generator)) + (1.0e-5 * static_cast<double>((i * N) + j));
        }
    }
    return test_case;
}

template<std::size_t N> auto make_fixed_pr(const blocking_case<N>& test_case)
{
    using Model = fug::PengRobinson<N>;
    std::array<typename Model::SpeciesInput, N> inputs{};
    for (std::size_t i = 0; i < N; ++i) {
        inputs[i] = {
            .T_c = test_case.species[i].T_c, .P_c = test_case.species[i].P_c, .omega = test_case.species[i].omega};
    }
    return Model{inputs, test_case.kij};
}

template<std::size_t N> auto make_dynamic_pr(const blocking_case<N>& test_case)
{
    using Model = fug::PengRobinson<>;
    std::vector<Model::SpeciesInput> inputs(N);
    for (std::size_t i = 0; i < N; ++i) {
        inputs[i] = {
            .T_c = test_case.species[i].T_c, .P_c = test_case.species[i].P_c, .omega = test_case.species[i].omega};
    }
    return Model{std::span<const typename Model::SpeciesInput>{inputs}, std::span<const double>{test_case.kij}};
}

} // namespace

// Test-only copy of the pre-optimization runtime loop. Keeping runtime-sized
// storage and loop bounds in the oracle avoids compiler unrolling differences
// between fixed and dynamic model types when Enzyme differentiates the code.
// External linkage is required because Enzyme instantiates differentiation intrinsics with this type.
// NOLINTNEXTLINE(misc-use-internal-linkage)
class legacy_runtime_pr {
public:
    template<std::size_t N>
    explicit legacy_runtime_pr(const blocking_case<N>& test_case) : b_(N), p_(N), q_(N), a_(N * N)
    {
        constexpr double R = fug::ideal_gas_constant<double>;
        const double s8 = std::sqrt(8.0);
        const double eta_c = 1.0 / (1.0 + std::cbrt(4.0 - s8) + std::cbrt(4.0 + s8));
        const double omega_a = (8.0 + (40.0 * eta_c)) / (49.0 - (37.0 * eta_c));
        const double omega_b = eta_c / (3.0 + eta_c);
        std::vector<double> a0(N);
        for (std::size_t i = 0; i < N; ++i) {
            const double w = test_case.species[i].omega;
            const double m = w <= 0.491 ? 0.37464 + (1.54226 * w) - (0.26992 * w * w)
                                        : 0.379642 + (1.48503 * w) - (0.164423 * w * w) + (0.016666 * w * w * w);
            const double RTc = R * test_case.species[i].T_c;
            a0[i] = omega_a * RTc * RTc / test_case.species[i].P_c;
            b_[i] = omega_b * RTc / test_case.species[i].P_c;
            p_[i] = 1.0 + m;
            q_[i] = m / std::sqrt(test_case.species[i].T_c);
        }
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                const double k_sym = 0.5 * (test_case.kij[(i * N) + j] + test_case.kij[(j * N) + i]);
                a_[(i * N) + j] = (1.0 - k_sym) * std::sqrt(a0[i] * a0[j]);
            }
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return b_.size(); }

    template<std::floating_point Number>
    [[nodiscard]] Number calc_helmholtz(const Number c, const Number* x, const Number T) const
    {
        constexpr double d1 = 1.0 + std::numbers::sqrt2;
        constexpr double d2 = 1.0 - std::numbers::sqrt2;
        const Number R = fug::ideal_gas_constant<Number>;
        const std::size_t n = size();
        const Number sT = std::sqrt(T);
        Number am{0};
        Number bm{0};
        for (std::size_t i = 0; i < n; ++i) {
            bm += x[i] * b_[i];
            const Number ti = x[i] * std::abs(p_[i] - (q_[i] * sT));
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += x[j] * std::abs(p_[j] - (q_[j] * sT)) * a_[(i * n) + j];
            }
            am += ti * row;
        }
        const Number bc = bm * c;
        const Number psi1 = -std::log(Number{1} - bc);
        const Number psi2 = std::log(((d1 * bc) + Number{1}) / ((d2 * bc) + Number{1})) / (bm * (d1 - d2));
        return (R * T * psi1) - (am * psi2);
    }

    template<std::floating_point Number>
    [[nodiscard]] Number calc_helmholtz_density(const Number* rho_i, const Number T) const
    {
        constexpr double d1 = 1.0 + std::numbers::sqrt2;
        constexpr double d2 = 1.0 - std::numbers::sqrt2;
        const Number R = fug::ideal_gas_constant<Number>;
        const std::size_t n = size();
        const Number sT = std::sqrt(T);
        Number c{0};
        Number bc{0};
        Number ac{0};
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
            bc += rho_i[i] * b_[i];
            const Number ti = rho_i[i] * std::abs(p_[i] - (q_[i] * sT));
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += rho_i[j] * std::abs(p_[j] - (q_[j] * sT)) * a_[(i * n) + j];
            }
            ac += ti * row;
        }
        const Number psi1 = -std::log(Number{1} - bc);
        return (R * T * c * psi1) -
               (ac * std::log(((d1 * bc) + Number{1}) / ((d2 * bc) + Number{1})) / (bc * (d1 - d2)));
    }

    template<std::floating_point Number>
    void calc_partial_helmholtz(const Number* rho_i, const Number T, Number* out) const
    {
        Number c{0};
        for (std::size_t i = 0; i < size(); ++i) {
            c += rho_i[i];
        }
        const Number scale = calc_helmholtz_density(rho_i, T) / c;
        for (std::size_t i = 0; i < size(); ++i) {
            out[i] = rho_i[i] * scale;
        }
    }

private:
    std::vector<double> b_;
    std::vector<double> p_;
    std::vector<double> q_;
    std::vector<double> a_;
};

namespace {

template<std::size_t N> auto make_dynamic_blocking_eos(const blocking_case<N>& test_case)
{
    using Ideal = fug::ConstantCp<>;
    std::vector<Ideal::SpeciesInput> inputs(N);
    for (std::size_t i = 0; i < N; ++i) {
        inputs[i] = {.T_ref = 298.15,
                     .p_ref = 1.0e5 + (100.0 * static_cast<double>(i)),
                     .c_p = test_case.species[i].c_p,
                     .h_ref = 125.0 * static_cast<double>(i),
                     .s_ref = 130.0 + (0.5 * static_cast<double>(i))};
    }
    return fug::EoS{Ideal{std::span<const typename Ideal::SpeciesInput>{inputs}}, make_dynamic_pr(test_case)};
}

template<std::size_t N> auto make_legacy_blocking_eos(const blocking_case<N>& test_case)
{
    using Ideal = fug::ConstantCp<>;
    std::vector<Ideal::SpeciesInput> inputs(N);
    for (std::size_t i = 0; i < N; ++i) {
        inputs[i] = {.T_ref = 298.15,
                     .p_ref = 1.0e5 + (100.0 * static_cast<double>(i)),
                     .c_p = test_case.species[i].c_p,
                     .h_ref = 125.0 * static_cast<double>(i),
                     .s_ref = 130.0 + (0.5 * static_cast<double>(i))};
    }
    return fug::EoS{Ideal{std::span<const typename Ideal::SpeciesInput>{inputs}}, legacy_runtime_pr{test_case}};
}

#ifdef __NO_MATH_ERRNO__
constexpr bool relaxed_arithmetic_build = true;
#else
constexpr bool relaxed_arithmetic_build = false;
#endif

template<std::floating_point Number>
void expect_relaxed_arithmetic_close(const std::string_view quantity, const Number actual, const Number expected)
{
    const Number magnitude = std::abs(expected) > Number{1} ? std::abs(expected) : Number{1};
    const Number tolerance = Number{4096} * std::numeric_limits<Number>::epsilon() * magnitude;
    expect(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) <= tolerance)
        << quantity << ": actual=" << actual << ", expected=" << expected << ", tolerance=" << tolerance;
}

void expect_same_double(const std::string_view quantity, const double actual, const double expected)
{
    if constexpr (relaxed_arithmetic_build) {
        // release-max explicitly permits reassociation and changed rounding; it is
        // an additional compatibility build, not the numerical baseline.
        expect_relaxed_arithmetic_close(quantity, actual, expected);
        return;
    }
    const auto actual_bits = std::bit_cast<std::uint64_t>(actual);
    const auto expected_bits = std::bit_cast<std::uint64_t>(expected);
    const auto ulp_distance = actual_bits > expected_bits ? actual_bits - expected_bits : expected_bits - actual_bits;
    expect(actual_bits == expected_bits) << quantity << ": actual=" << actual << ", expected=" << expected
                                         << ", ulp distance=" << ulp_distance;
}

std::uint64_t ordered_double_bits(const double value)
{
    constexpr std::uint64_t sign_bit = std::uint64_t{1} << 63U;
    const auto bits = std::bit_cast<std::uint64_t>(value);
    return (bits & sign_bit) != 0U ? ~bits : bits | sign_bit;
}

void expect_within_ulps(const std::string_view quantity, const double actual, const double expected,
                        const std::uint64_t max_ulps)
{
    if constexpr (relaxed_arithmetic_build) {
        expect_relaxed_arithmetic_close(quantity, actual, expected);
        return;
    }
    const auto actual_bits = ordered_double_bits(actual);
    const auto expected_bits = ordered_double_bits(expected);
    const auto ulp_distance = actual_bits > expected_bits ? actual_bits - expected_bits : expected_bits - actual_bits;
    expect(std::isfinite(actual) && std::isfinite(expected) && ulp_distance <= max_ulps)
        << quantity << ": actual=" << actual << ", expected=" << expected << ", ulp distance=" << ulp_distance
        << ", budget=" << max_ulps;
}

template<std::floating_point Number>
void expect_same_number(const std::string_view quantity, const Number actual, const Number expected)
{
    if constexpr (std::same_as<Number, double>) {
        expect_same_double(quantity, actual, expected);
    }
    else if constexpr (relaxed_arithmetic_build) {
        expect_relaxed_arithmetic_close(quantity, actual, expected);
    }
    else {
        expect(actual == expected) << quantity << ": actual=" << actual << ", expected=" << expected;
    }
}

template<std::size_t N, std::floating_point Number> void check_molar_blocking_characterization()
{
    const auto test_case = make_blocking_case<N>();
    const auto fixed = make_fixed_pr(test_case);
    const auto dynamic = make_dynamic_pr(test_case);
    std::array<Number, N> x{};
    std::array<Number, N> rho{};
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = static_cast<Number>(test_case.x[i]);
        rho[i] = static_cast<Number>(test_case.rho[i]);
    }
    const auto c = static_cast<Number>(test_case.c);
    const auto T = static_cast<Number>(test_case.T);
    expect_same_number("molar runtime/static", dynamic.calc_helmholtz(c, x.data(), T),
                       fixed.calc_helmholtz(c, x.data(), T));
    expect_same_number("density runtime/static", dynamic.calc_helmholtz_density(rho.data(), T),
                       fixed.calc_helmholtz_density(rho.data(), T));

    std::array<Number, N> fixed_partial{};
    std::array<Number, N> dynamic_partial{};
    fixed.calc_partial_helmholtz(rho.data(), T, fixed_partial.data());
    dynamic.calc_partial_helmholtz(rho.data(), T, dynamic_partial.data());
    for (std::size_t i = 0; i < N; ++i) {
        expect_same_number("partial runtime/static", dynamic_partial[i], fixed_partial[i]);
    }
    if constexpr (std::same_as<Number, double>) {
        const legacy_runtime_pr legacy{test_case};
        expect_same_double("molar runtime/legacy", dynamic.calc_helmholtz(c, x.data(), T),
                           legacy.calc_helmholtz(c, x.data(), T));
        expect_same_double("density runtime/legacy", dynamic.calc_helmholtz_density(rho.data(), T),
                           legacy.calc_helmholtz_density(rho.data(), T));
        std::array<double, N> legacy_partial{};
        legacy.calc_partial_helmholtz(rho.data(), T, legacy_partial.data());
        for (std::size_t i = 0; i < N; ++i) {
            expect_same_double("partial runtime/legacy", dynamic_partial[i], legacy_partial[i]);
        }
    }
}

template<std::size_t N> void check_blocking_public_properties()
{
    const auto test_case = make_blocking_case<N>();
    const auto legacy = make_legacy_blocking_eos(test_case);
    const auto dynamic = make_dynamic_blocking_eos(test_case);
    auto x = test_case.x;
    // These ULP budgets are local regression limits for the blocked molar path.
    // The independent scientific-oracle tests below retain their established tolerances.
    expect_same_double("pressure runtime/legacy", fug::calc_pressure(dynamic, test_case.c, x, test_case.T),
                       fug::calc_pressure(legacy, test_case.c, x, test_case.T));
    expect_within_ulps("cp runtime/legacy", fug::calc_cp(dynamic, test_case.c, x, test_case.T),
                       fug::calc_cp(legacy, test_case.c, x, test_case.T), 16);
    expect_within_ulps("sound speed runtime/legacy",
                       fug::calc_sound_speed_squared(dynamic, test_case.c, x, test_case.T, 0.035),
                       fug::calc_sound_speed_squared(legacy, test_case.c, x, test_case.T, 0.035), 16);

    const double invT = 1.0 / test_case.T;
    const auto check_lambda = [&](const std::string_view name, const double actual, const double expected) {
        expect_within_ulps(name, actual, expected, 32);
    };
    check_lambda("lambda(0,0) runtime/legacy",
                 fug::detail::calc_lambda<0, 0>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<0, 0>(legacy.residual(), test_case.c, x.data(), invT));
    check_lambda("lambda(1,0) runtime/legacy",
                 fug::detail::calc_lambda<1, 0>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<1, 0>(legacy.residual(), test_case.c, x.data(), invT));
    check_lambda("lambda(0,1) runtime/legacy",
                 fug::detail::calc_lambda<0, 1>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<0, 1>(legacy.residual(), test_case.c, x.data(), invT));
    check_lambda("lambda(0,2) runtime/legacy",
                 fug::detail::calc_lambda<0, 2>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<0, 2>(legacy.residual(), test_case.c, x.data(), invT));
    check_lambda("lambda(1,1) runtime/legacy",
                 fug::detail::calc_lambda<1, 1>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<1, 1>(legacy.residual(), test_case.c, x.data(), invT));
    check_lambda("lambda(2,0) runtime/legacy",
                 fug::detail::calc_lambda<2, 0>(dynamic.residual(), test_case.c, x.data(), invT),
                 fug::detail::calc_lambda<2, 0>(legacy.residual(), test_case.c, x.data(), invT));

    std::array<double, N> legacy_fugacity{};
    std::array<double, N> dynamic_fugacity{};
    fug::calc_fugacity(legacy, std::span<const double, N>{test_case.rho}, test_case.T,
                       std::span<double, N>{legacy_fugacity});
    fug::calc_fugacity(dynamic, std::span<const double, N>{test_case.rho}, test_case.T,
                       std::span<double, N>{dynamic_fugacity});
    for (std::size_t i = 0; i < N; ++i) {
        expect_same_double("fugacity runtime/legacy", dynamic_fugacity[i], legacy_fugacity[i]);
    }
}

#if defined(NDEBUG) && !defined(__FAST_MATH__) && !defined(__NO_MATH_ERRNO__)

blocking_case<16> make_cp_dx_finite_difference_case()
{
    constexpr std::size_t n = 16;
    constexpr std::uint64_t seed = 0xF1D1FFB10C5EED11;
    // A reproducible regression sequence is intentional, not security-sensitive randomness.
    // NOLINTNEXTLINE(bugprone-random-generator-seed)
    std::mt19937_64 generator{seed};
    std::uniform_real_distribution<double> unit{0.0, 1.0};
    blocking_case<n> test_case;
    double x_sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        test_case.species[i] = {.T_c = 180.0 + (420.0 * unit(generator)) + static_cast<double>(i),
                                .P_c = 3.5e6 + (4.0e6 * unit(generator)) + (1000.0 * static_cast<double>(i)),
                                .omega = 0.02 + (0.65 * unit(generator)),
                                .c_p = 24.0 + (18.0 * unit(generator)) + (0.125 * static_cast<double>(i))};
        test_case.x[i] = 0.1 + unit(generator) + (0.01 * static_cast<double>(i));
        x_sum += test_case.x[i];
    }
    for (double& xi : test_case.x) {
        xi /= x_sum;
    }
    test_case.c = 1250.0;
    test_case.T = 415.0;
    for (std::size_t i = 0; i < n; ++i) {
        test_case.rho[i] = test_case.c * test_case.x[i];
        for (std::size_t j = 0; j < n; ++j) {
            test_case.kij[(i * n) + j] =
                (i == j) ? 0.0 : -0.08 + (0.20 * unit(generator)) + (1.0e-5 * static_cast<double>((i * n) + j));
        }
    }
    return test_case;
}

void check_cp_dx_against_five_point_finite_difference()
{
    constexpr std::size_t n = 16;
    constexpr double step = 5.0e-4;
    constexpr double max_scaled_error = 2.0e-7;
    const auto test_case = make_cp_dx_finite_difference_case();
    const auto eos = make_dynamic_blocking_eos(test_case);
    auto x = test_case.x;
    std::array<double, n> gradient{};
    fug::calc_cp_dx(eos, test_case.c, x, test_case.T, gradient);

    for (std::size_t direction = 0; direction < n; ++direction) {
        const auto cp_at_offset = [&](const double offset) {
            auto perturbed_x = test_case.x;
            perturbed_x[direction] += offset;
            return fug::calc_cp(eos, test_case.c, perturbed_x, test_case.T);
        };
        const double finite_difference = (-cp_at_offset(2.0 * step) + (8.0 * cp_at_offset(step)) -
                                          (8.0 * cp_at_offset(-step)) + cp_at_offset(-2.0 * step)) /
                                         (12.0 * step);
        const double magnitude = std::abs(gradient[direction]) > std::abs(finite_difference)
                                     ? std::abs(gradient[direction])
                                     : std::abs(finite_difference);
        const double scale = magnitude > 1.0 ? magnitude : 1.0;
        const double scaled_error = std::abs(gradient[direction] - finite_difference) / scale;
        expect(std::isfinite(gradient[direction]) && std::isfinite(finite_difference) &&
               scaled_error <= max_scaled_error)
            << "direction=" << direction << ", Enzyme=" << gradient[direction]
            << ", finite difference=" << finite_difference << ", scaled error=" << scaled_error
            << ", budget=" << max_scaled_error;
    }
}

#endif

} // namespace

// Test entry points intentionally let assertion failures escape to the runner.
// NOLINTNEXTLINE(bugprone-exception-escape)
int main()
{
    const suite<"peng_robinson_unary_contracts"> unary_contracts = [] {
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

    const suite<"peng_robinson"> pr_suite = [] {
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

        "runtime kernels preserve behavior around molar blocking boundaries"_test = [] {
            check_molar_blocking_characterization<7, double>();
            check_molar_blocking_characterization<8, double>();
            check_molar_blocking_characterization<9, double>();
            check_molar_blocking_characterization<10, double>();
            check_molar_blocking_characterization<15, double>();
            check_molar_blocking_characterization<16, double>();
            check_molar_blocking_characterization<17, double>();
            check_molar_blocking_characterization<23, double>();
            check_molar_blocking_characterization<24, double>();
            check_molar_blocking_characterization<25, double>();
            check_molar_blocking_characterization<31, double>();
            check_molar_blocking_characterization<32, double>();
            check_molar_blocking_characterization<33, double>();
        };

        "non-double and static kernels remain scalar-path canaries"_test = [] {
            check_molar_blocking_characterization<9, float>();
            check_molar_blocking_characterization<10, float>();
            check_molar_blocking_characterization<17, float>();
            check_molar_blocking_characterization<9, long double>();
            check_molar_blocking_characterization<10, long double>();
            check_molar_blocking_characterization<17, long double>();
        };

        "runtime Enzyme properties preserve legacy order at dispatch and tail boundaries"_test = [] {
            check_blocking_public_properties<9>();
            check_blocking_public_properties<10>();
            check_blocking_public_properties<17>();
        };

#if defined(NDEBUG) && !defined(__FAST_MATH__) && !defined(__NO_MATH_ERRNO__)
        // Enzyme's reverse composition derivative is optimization-sensitive at
        // Debug/O1. Keep this independent oracle on the strict Release/O3
        // numerical baseline; release-max explicitly permits reassociation.
        "runtime calc_cp_dx matches five-point finite differences"_test = [] {
            check_cp_dx_against_five_point_finite_difference();
        };
#endif

        "runtime mixture kernels are covariant under species reversal"_test = [] {
            constexpr std::size_t n = 17;
            const auto original_case = make_blocking_case<n>();
            auto reversed_case = original_case;
            for (std::size_t i = 0; i < n; ++i) {
                reversed_case.species[i] = original_case.species[n - 1 - i];
                reversed_case.x[i] = original_case.x[n - 1 - i];
                reversed_case.rho[i] = original_case.rho[n - 1 - i];
                for (std::size_t j = 0; j < n; ++j) {
                    reversed_case.kij[(i * n) + j] = original_case.kij[((n - 1 - i) * n) + (n - 1 - j)];
                }
            }
            const auto original = make_dynamic_pr(original_case);
            const auto reversed = make_dynamic_pr(reversed_case);
            check_rel("permuted molar",
                      reversed.calc_helmholtz(reversed_case.c, reversed_case.x.data(), reversed_case.T),
                      original.calc_helmholtz(original_case.c, original_case.x.data(), original_case.T), 2e-15);
            check_rel("permuted density", reversed.calc_helmholtz_density(reversed_case.rho.data(), reversed_case.T),
                      original.calc_helmholtz_density(original_case.rho.data(), original_case.T), 2e-15);
            std::array<double, n> original_partial{};
            std::array<double, n> reversed_partial{};
            original.calc_partial_helmholtz(original_case.rho.data(), original_case.T, original_partial.data());
            reversed.calc_partial_helmholtz(reversed_case.rho.data(), reversed_case.T, reversed_partial.data());
            for (std::size_t i = 0; i < n; ++i) {
                check_rel("permuted partial", reversed_partial[i], original_partial[n - 1 - i], 2e-15);
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
            const auto b = static_cast<double>(pr_b(co2.T_c, co2.P_c));
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};
            for (const double c : {1.0, 100.0, 5000.0, 15000.0}) {
                for (const double T : {220.0, 320.0, 500.0}) {
                    const auto a_T = static_cast<double>(pr_aii(unary_inputs[0], T));
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
            const auto b = static_cast<double>(pr_b(co2.T_c, co2.P_c));
            const auto c_c = static_cast<double>(pr_eta_c()) / b;
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
