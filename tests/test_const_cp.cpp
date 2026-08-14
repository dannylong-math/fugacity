//
// Unit tests for the constant-cp ideal-gas model (synthesize::ConstantCp).
//
// Universal, ideal-gas, derivative, wrapper, precondition, deterministic-sweep,
// and static/dynamic checks are registered through support/eos_test_suite.hpp.
// What remains here is genuinely ConstantCp-specific: the per-species kernel
// algebra and the caloric properties checked against the *physical* reference
// data (h_ref, s_ref, c_p), which the generic helpers cannot know.
//
// Physical references (pure ideal gas, constant isobaric molar heat capacity
// c_p, referenced to a state (T_ref, p_ref)):
//     h(T)    = h_ref + c_p (T - T_ref)
//     s(T, c) = s_ref + (c_p - R) ln(T / T_ref) - R ln(c / c_ref)
//     c_p     = c_p  (constant),   with c_ref = p_ref / (R T_ref).
// At the reference state (T = T_ref, c = c_ref): h = h_ref and s = s_ref.
//
#include "derivative_test_harness.hpp"
#include "support/eos_test_suite.hpp"
#include "synthesize/core/core_calculations.hpp"
#include "synthesize/core/eos_pair.hpp"
#include "synthesize/core/numbers.hpp"
#include "synthesize/ideal_models/const_cp.hpp"
#include "synthesize/residual_models/no_residual.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

using namespace boost::ut;
using namespace synthesize_test;

namespace {

namespace ge = synthesize;

template<std::size_t N> using Input = typename ge::ConstantCp<N>::SpeciesInput;

// Build a complete EoS: a constant-cp ideal contribution + a vanishing residual.
template<std::size_t N> auto make_const_cp_eos(const std::array<Input<N>, N>& in)
{
    return ge::EoS<ge::ConstantCp<N>, ge::NoResidual<N>>{ge::ConstantCp<N>(in), ge::NoResidual<N>{}};
}

// A representative two-component parameter set. The two species deliberately use
// *different* reference temperatures and pressures to exercise per-species
// reference handling.
constexpr std::array<Input<2>, 2> binary_inputs{{
    {/*T_ref*/ .T_ref=300.0, /*p_ref*/ .p_ref=1.0e5, /*c_p*/ .c_p=29.1, /*h_ref*/ .h_ref=1500.0, /*s_ref*/ .s_ref=191.0},
    {/*T_ref*/ .T_ref=320.0, /*p_ref*/ .p_ref=9.0e4, /*c_p*/ .c_p=33.6, /*h_ref*/ .h_ref=-2200.0, /*s_ref*/ .s_ref=189.0},
}};

std::vector<ge::ConstantCp<>::SpeciesInput> dynamic_binary_inputs()
{
    std::vector<ge::ConstantCp<>::SpeciesInput> result;
    result.reserve(binary_inputs.size());
    for (const auto& input : binary_inputs) {
        result.push_back({.T_ref=input.T_ref, .p_ref=input.p_ref, .c_p=input.c_p, .h_ref=input.h_ref, .s_ref=input.s_ref});
    }
    return result;
}

auto make_dynamic_const_cp_eos()
{
    const auto inputs = dynamic_binary_inputs();
    return ge::EoS{ge::ConstantCp<>{std::span<const ge::ConstantCp<>::SpeciesInput>{inputs}},
                   ge::NoResidual<std::dynamic_extent>{inputs.size()}};
}

std::vector<eos_test_state> const_cp_contract_states()
{
    return {{.c = 40.0, .x = {0.5, 0.5}, .T = 280.0, .effective_molar_mass = 0.031, .label = "low-density"},
            {.c = 100.0, .x = {0.4, 0.6}, .T = 330.0, .effective_molar_mass = 0.031, .label = "representative"},
            {.c = 250.0, .x = {0.7, 0.3}, .T = 425.0, .effective_molar_mass = 0.031, .label = "warm"}};
}

constexpr eos_valid_domain const_cp_valid_domain{.c_min = 20.0,
                                                 .c_max = 260.0,
                                                 .T_min = 250.0,
                                                 .T_max = 600.0,
                                                 .minimum_mole_fraction = 0.01,
                                                 .seed = 0xC0FFEE,
                                                 .random_samples = 50};

// Single-species reference data for the caloric checks.
constexpr double T_ref = 300.0;
constexpr double p_ref = 1.0e5;
constexpr double cp_in = 29.1;
constexpr double h_ref = -1234.0;
constexpr double s_ref = 205.0;

constexpr std::array<Input<1>, 1> unary_inputs{{{.T_ref=T_ref, .p_ref=p_ref, .c_p=cp_in, .h_ref=h_ref, .s_ref=s_ref}}};

} // namespace

int main()
{
    suite<"const_cp"> const_cp = [] {
        const double R = ge::ideal_gas_constant<double>;

        auto dynamic_eos = make_dynamic_const_cp_eos();
        const auto fixture = eos_test_fixture{.contribution = dynamic_eos.ideal(),
                                              .eos = dynamic_eos,
                                              .states = const_cp_contract_states(),
                                              .domain = const_cp_valid_domain};
        register_eos_contract_tests(fixture);
        register_ideal_gas_contract_tests(fixture);
        register_static_dynamic_equivalence_tests(make_const_cp_eos<2>(binary_inputs), make_dynamic_const_cp_eos(),
                                                  const_cp_contract_states());

        // ===================================================================
        // ConstantCp-specific caloric properties vs the PHYSICAL reference data
        // (the generic helpers don't know h_ref/s_ref/c_p). This is the guard
        // that the Helmholtz formula itself is correct.
        // ===================================================================
        "single-species reference state"_test = [&] {
            auto eos = make_const_cp_eos<1>(unary_inputs);
            const double c_ref = p_ref / (R * T_ref);
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};

            check_rel("calc_cp       == c_p", ge::calc_cp(eos, c_ref, xs, T_ref), cp_in, 1e-9);
            check_rel("calc_enthalpy == h_ref", ge::calc_enthalpy(eos, c_ref, xs, T_ref), h_ref, 1e-9);
            check_rel("calc_entropy  == s_ref", ge::calc_entropy(eos, c_ref, xs, T_ref), s_ref, 1e-9);
        };

        "single-species off-reference caloric"_test = [&] {
            auto eos = make_const_cp_eos<1>(unary_inputs);
            const double c_ref = p_ref / (R * T_ref);
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};

            for (const double T : {275.0, 350.0, 425.0}) {
                for (const double c : {c_ref, 0.5 * c_ref, 3.0 * c_ref}) {
                    const double h_expected = h_ref + (cp_in * (T - T_ref));
                    const double s_expected = s_ref + ((cp_in - R) * std::log(T / T_ref)) - (R * std::log(c / c_ref));
                    check_rel("calc_enthalpy(T)", ge::calc_enthalpy(eos, c, xs, T), h_expected, 1e-9);
                    check_rel("calc_cp(T,c)", ge::calc_cp(eos, c, xs, T), cp_in, 1e-9);
                    check_rel("calc_entropy(T,c)", ge::calc_entropy(eos, c, xs, T), s_expected, 1e-9);
                }
            }
        };
    };

    return ::boost::ut::cfg<>.run();
}
