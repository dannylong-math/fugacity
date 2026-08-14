//
// Unit tests for the NASA-7 polynomial ideal-gas model (synthesize::Nasa7).
//
// Universal, ideal-gas, derivative, wrapper, precondition, deterministic-sweep,
// and static/dynamic checks are registered through support/eos_test_suite.hpp.
// What remains here is genuinely NASA7-specific: the per-species kernel algebra
// and the caloric properties checked against the *physical* NASA-7 relations,
// which the generic helpers cannot know.
//
// NASA-7 standard-state relations for a single species (T in [T_low, T_high]),
// with R the molar gas constant:
//
//   c_p(T) / R       = a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4
//   h(T)  / (R T)    = a0 + a1/2 T + a2/3 T^2 + a3/4 T^3 + a4/5 T^4 + a5/T
//   s(T)  / R        = a0 ln(T) + a1 T + a2/2 T^2 + a3/3 T^3 + a4/4 T^4 + a6
//
// The standard-state entropy s(T) above is the value at the reference
// concentration c_ref = p_ref / (R T_ref); at a general concentration c the
// model entropy carries the ideal-gas mixing term  -R ln(c / c_ref).
// Enthalpy and c_p are pressure/concentration independent for an ideal gas.
//
#include "derivative_test_harness.hpp"
#include "support/eos_test_suite.hpp"
#include "synthesize/core/core_calculations.hpp"
#include "synthesize/core/eos_pair.hpp"
#include "synthesize/core/numbers.hpp"
#include "synthesize/ideal_models/nasa7.hpp"
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

template<std::size_t N> using Input = typename ge::Nasa7<N>::SpeciesInput;

// Build a complete EoS: a NASA-7 ideal contribution + a vanishing residual.
template<std::size_t N> auto make_nasa7_eos(const std::array<Input<N>, N>& in)
{
    return ge::EoS<ge::Nasa7<N>, ge::NoResidual<N>>{ge::Nasa7<N>(in), ge::NoResidual<N>{}};
}

// Representative low-temperature (200-1000 K) NASA-7 coefficient sets for N2 and
// O2. The two species deliberately use *different* reference temperatures and
// pressures to exercise per-species reference handling.
constexpr std::array<Input<2>, 2> binary_inputs{{
    {/*a0*/ .a0 = 3.53100528,
     /*a1*/ .a1 = -1.23660988e-4,
     /*a2*/ .a2 = -5.02999433e-7,
     /*a3*/ .a3 = 2.43530612e-9,
     /*a4*/ .a4 = -1.40881235e-12,
     /*a5*/ .a5 = -1046.97628,
     /*a6*/ .a6 = 2.96747038,
     /*T_ref*/ .T_ref = 300.0,
     /*p_ref*/ .p_ref = 1.0e5},
    {/*a0*/ .a0 = 3.78245636,
     /*a1*/ .a1 = -2.99673416e-3,
     /*a2*/ .a2 = 9.84730201e-6,
     /*a3*/ .a3 = -9.68129509e-9,
     /*a4*/ .a4 = 3.24372837e-12,
     /*a5*/ .a5 = -1063.94356,
     /*a6*/ .a6 = 3.65767573,
     /*T_ref*/ .T_ref = 320.0,
     /*p_ref*/ .p_ref = 9.0e4},
}};

std::vector<ge::Nasa7<>::SpeciesInput> dynamic_binary_inputs()
{
    std::vector<ge::Nasa7<>::SpeciesInput> result;
    result.reserve(binary_inputs.size());
    for (const auto& input : binary_inputs) {
        result.push_back(
            {.a0=input.a0, .a1=input.a1, .a2=input.a2, .a3=input.a3, .a4=input.a4, .a5=input.a5, .a6=input.a6, .T_ref=input.T_ref, .p_ref=input.p_ref});
    }
    return result;
}

auto make_dynamic_nasa7_eos()
{
    const auto inputs = dynamic_binary_inputs();
    return ge::EoS{ge::Nasa7<>{std::span<const ge::Nasa7<>::SpeciesInput>{inputs}},
                   ge::NoResidual<std::dynamic_extent>{inputs.size()}};
}

std::vector<eos_test_state> nasa7_contract_states()
{
    return {{.c = 40.0, .x = {0.5, 0.5}, .T = 280.0, .effective_molar_mass = 0.029, .label = "low-density"},
            {.c = 100.0, .x = {0.4, 0.6}, .T = 300.0, .effective_molar_mass = 0.029, .label = "representative"},
            {.c = 250.0, .x = {0.7, 0.3}, .T = 410.0, .effective_molar_mass = 0.029, .label = "warm"}};
}

constexpr eos_valid_domain nasa7_valid_domain{.c_min = 20.0,
                                              .c_max = 260.0,
                                              .T_min = 250.0,
                                              .T_max = 600.0,
                                              .minimum_mole_fraction = 0.01,
                                              .seed = 0xC0FFEE,
                                              .random_samples = 50};

// Single-species reference data (N2) for the caloric checks.
constexpr double T_ref = 300.0;
constexpr double p_ref = 1.0e5;

constexpr std::array<Input<1>, 1> unary_inputs{{{/*a0*/ .a0 = 3.53100528,
                                                 /*a1*/ .a1 = -1.23660988e-4,
                                                 /*a2*/ .a2 = -5.02999433e-7,
                                                 /*a3*/ .a3 = 2.43530612e-9,
                                                 /*a4*/ .a4 = -1.40881235e-12,
                                                 /*a5*/ .a5 = -1046.97628,
                                                 /*a6*/ .a6 = 2.96747038,
                                                 .T_ref = T_ref,
                                                 .p_ref = p_ref}}};

// --- Closed-form NASA-7 reference relations (in double) --------------------
double nasa7_cp(const Input<1>& in, double T)
{
    const double R = ge::ideal_gas_constant<double>;
    return R * (in.a0 + (in.a1 * T) + (in.a2 * T * T) + (in.a3 * T * T * T) + (in.a4 * T * T * T * T));
}

double nasa7_enthalpy(const Input<1>& in, double T)
{
    const double R = ge::ideal_gas_constant<double>;
    return R * T *
           ((in.a4 / 5 * T * T * T * T) + (in.a3 / 4 * T * T * T) + (in.a2 / 3 * T * T) + (in.a1 / 2 * T) + in.a0 +
            (in.a5 / T));
}

// Standard-state entropy: the model entropy at the reference concentration c_ref.
double nasa7_entropy_std(const Input<1>& in, double T)
{
    const double R = ge::ideal_gas_constant<double>;
    return R * ((in.a4 / 4 * T * T * T * T) + (in.a3 / 3 * T * T * T) + (in.a2 / 2 * T * T) + (in.a1 * T) +
                (in.a0 * std::log(T)) + in.a6);
}

} // namespace

int main()
{
    suite<"nasa7"> nasa7 = [] {
        const double R = ge::ideal_gas_constant<double>;

        auto dynamic_eos = make_dynamic_nasa7_eos();
        const auto fixture = eos_test_fixture{.contribution = dynamic_eos.ideal(),
                                              .eos = dynamic_eos,
                                              .states = nasa7_contract_states(),
                                              .domain = nasa7_valid_domain};
        register_eos_contract_tests(fixture);
        register_ideal_gas_contract_tests(fixture);
        register_static_dynamic_equivalence_tests(make_nasa7_eos<2>(binary_inputs), make_dynamic_nasa7_eos(),
                                                  nasa7_contract_states());

        // ===================================================================
        // NASA7-specific caloric properties vs the closed-form NASA-7 relations
        // at the reference state (T = T_ref, c = c_ref). This is the guard that
        // the Helmholtz formula itself reproduces the intended polynomial.
        // ===================================================================
        "single-species reference state"_test = [&] {
            auto eos = make_nasa7_eos<1>(unary_inputs);
            const Input<1>& in = unary_inputs[0];
            const double c_ref = p_ref / (R * T_ref);
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};

            check_rel("calc_cp       == NASA7 c_p", ge::calc_cp(eos, c_ref, xs, T_ref), nasa7_cp(in, T_ref), 1e-9);
            check_rel("calc_enthalpy == NASA7 h", ge::calc_enthalpy(eos, c_ref, xs, T_ref), nasa7_enthalpy(in, T_ref),
                      1e-9);
            check_rel("calc_entropy  == NASA7 s", ge::calc_entropy(eos, c_ref, xs, T_ref), nasa7_entropy_std(in, T_ref),
                      1e-9);
        };

        // ===================================================================
        // NASA7-specific caloric properties off the reference state: the
        // polynomial relations hold at any T. The standard-state entropy is
        // referenced to the standard pressure p_ref, so at a general state the
        // entropy carries the ideal-gas pressure term -R ln(p / p_ref), with
        // p = c R T  =>  p / p_ref = (c T) / (c_ref T_ref). (Enthalpy and c_p
        // are pressure/concentration independent for an ideal gas.)
        // ===================================================================
        "single-species off-reference caloric"_test = [&] {
            auto eos = make_nasa7_eos<1>(unary_inputs);
            const Input<1>& in = unary_inputs[0];
            const double c_ref = p_ref / (R * T_ref);
            const std::array<double, 1> x{1.0};
            const std::span<const double, 1> xs{x};

            for (const double T : {250.0, 350.0, 600.0}) {
                for (const double c : {c_ref, 0.5 * c_ref, 3.0 * c_ref}) {
                    const double h_expected = nasa7_enthalpy(in, T);
                    const double cp_expected = nasa7_cp(in, T);
                    const double s_expected = nasa7_entropy_std(in, T) - (R * std::log((c * T) / (c_ref * T_ref)));
                    check_rel("calc_enthalpy(T)", ge::calc_enthalpy(eos, c, xs, T), h_expected, 1e-9);
                    check_rel("calc_cp(T,c)", ge::calc_cp(eos, c, xs, T), cp_expected, 1e-9);
                    check_rel("calc_entropy(T,c)", ge::calc_entropy(eos, c, xs, T), s_expected, 1e-9);
                }
            }
        };
    };

    return ::boost::ut::cfg<>.run();
}
