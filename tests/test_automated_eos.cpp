#include "support/eos_test_suite.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/residual_models/no_residual.hpp"
#include "fugacity/residual_models/peng_robinson.hpp"
#include "fugacity/residual_models/van_der_waals.hpp"

#include <array>
#include <boost/ut.hpp>
#include <span>
#include <type_traits>
#include <vector>

using namespace boost::ut;
using namespace fugacity_test;

namespace {
namespace fug = fugacity;

using FixedIdeal = fug::ConstantCp<2>;
using DynamicIdeal = fug::ConstantCp<std::dynamic_extent>;
using FixedUnaryIdeal = fug::ConstantCp<1>;

constexpr std::array<FixedIdeal::SpeciesInput, 2> fixed_inputs{{
    {.T_ref = 300.0, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 1500.0, .s_ref = 191.0},
    {.T_ref = 320.0, .p_ref = 9.0e4, .c_p = 33.6, .h_ref = -2200.0, .s_ref = 189.0},
}};

std::vector<DynamicIdeal::SpeciesInput> make_dynamic_inputs()
{
    return {{.T_ref = 300.0, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 1500.0, .s_ref = 191.0},
            {.T_ref = 320.0, .p_ref = 9.0e4, .c_p = 33.6, .h_ref = -2200.0, .s_ref = 189.0}};
}

auto make_fixed_eos() { return fug::EoS{FixedIdeal{fixed_inputs}, fug::NoResidual<2>{}}; }

auto make_dynamic_eos()
{
    const auto inputs = make_dynamic_inputs();
    return fug::EoS{DynamicIdeal{std::span<const DynamicIdeal::SpeciesInput>{inputs}},
                   fug::NoResidual<std::dynamic_extent>{inputs.size()}};
}

auto make_fixed_unary_eos()
{
    constexpr std::array<FixedUnaryIdeal::SpeciesInput, 1> inputs{{
        {.T_ref = 300.0, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 1500.0, .s_ref = 191.0},
    }};
    return fug::EoS{FixedUnaryIdeal{inputs}, fug::NoResidual<1>{}};
}

auto make_dynamic_unary_eos()
{
    const std::vector<DynamicIdeal::SpeciesInput> inputs{
        {.T_ref = 300.0, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 1500.0, .s_ref = 191.0}};
    return fug::EoS{DynamicIdeal{std::span<const DynamicIdeal::SpeciesInput>{inputs}},
                   fug::NoResidual<std::dynamic_extent>{inputs.size()}};
}

auto make_fixed_vdw_eos()
{
    using Residual = fug::VanDerWaals<2>;
    const std::array<Residual::SpeciesInput, 2> inputs{
        {{.T_c = 126.192, .P_c = 3.3958e6}, {.T_c = 304.1282, .P_c = 7.3773e6}}};
    return fug::EoS{FixedIdeal{fixed_inputs}, Residual{inputs}};
}

auto make_dynamic_vdw_eos()
{
    using Residual = fug::VanDerWaals<std::dynamic_extent>;
    const auto ideal_inputs = make_dynamic_inputs();
    const std::vector<Residual::SpeciesInput> residual_inputs{{.T_c = 126.192, .P_c = 3.3958e6},
                                                              {.T_c = 304.1282, .P_c = 7.3773e6}};
    return fug::EoS{DynamicIdeal{std::span<const DynamicIdeal::SpeciesInput>{ideal_inputs}},
                   Residual{std::span<const Residual::SpeciesInput>{residual_inputs}}};
}

auto make_fixed_pr_eos()
{
    using Residual = fug::PengRobinson<2>;
    const std::array<Residual::SpeciesInput, 2> inputs{{
        {.T_c = 190.564, .P_c = 4.5992e6, .omega = 0.011},
        {.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394},
    }};
    constexpr std::array<double, 4> kij{0.0, 0.09, 0.09, 0.0};
    return fug::EoS{FixedIdeal{fixed_inputs}, Residual{inputs, kij}};
}

auto make_dynamic_pr_eos()
{
    using Residual = fug::PengRobinson<std::dynamic_extent>;
    const auto ideal_inputs = make_dynamic_inputs();
    const std::vector<Residual::SpeciesInput> residual_inputs{
        {.T_c = 190.564, .P_c = 4.5992e6, .omega = 0.011},
        {.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394},
    };
    constexpr std::array<double, 4> kij{0.0, 0.09, 0.09, 0.0};
    return fug::EoS{DynamicIdeal{std::span<const DynamicIdeal::SpeciesInput>{ideal_inputs}},
                   Residual{std::span<const Residual::SpeciesInput>{residual_inputs}, std::span{kij}}};
}

std::vector<eos_test_state> make_states()
{
    return {{.c = 15.0, .x = {0.25, 0.75}, .T = 260.0, .effective_molar_mass = 0.031, .label = "low-density"},
            {.c = 100.0, .x = {0.4, 0.6}, .T = 330.0, .effective_molar_mass = 0.031, .label = "representative"},
            {.c = 240.0, .x = {0.8, 0.2}, .T = 580.0, .effective_molar_mass = 0.031, .label = "warm"}};
}

constexpr eos_valid_domain valid_domain{.c_min = 10.0,
                                        .c_max = 250.0,
                                        .T_min = 240.0,
                                        .T_max = 600.0,
                                        .minimum_mole_fraction = 0.01,
                                        .seed = 0xC0FFEE,
                                        .random_samples = 50};
} // namespace

int main()
{
    suite<"automated_unary_eos_contracts"> unary = [] {
        auto dynamic_eos = make_dynamic_unary_eos();
        const eos_valid_domain domain{.c_min = 10.0,
                                      .c_max = 250.0,
                                      .T_min = 240.0,
                                      .T_max = 600.0,
                                      .minimum_mole_fraction = 0.01,
                                      .seed = 0xC0FFEE,
                                      .random_samples = 50};
        const std::vector<eos_test_state> states{
            {.c = 100.0, .x = {1.0}, .T = 330.0, .effective_molar_mass = 0.028, .label = "unary"}};
        const auto fixture = eos_test_fixture{
            .contribution = dynamic_eos.ideal(), .eos = dynamic_eos, .states = states, .domain = domain};
        register_eos_contract_tests(fixture);
        register_ideal_gas_contract_tests(fixture);
        register_static_dynamic_equivalence_tests(make_fixed_unary_eos(), make_dynamic_unary_eos(), states);
    };

    suite<"automated_eos_contracts"> automated = [] {
        auto dynamic_eos = make_dynamic_eos();
        auto fixture = eos_test_fixture{
            .contribution = dynamic_eos.ideal(), .eos = dynamic_eos, .states = make_states(), .domain = valid_domain};
        register_eos_contract_tests(fixture);
        register_ideal_gas_contract_tests(fixture);
        register_static_dynamic_equivalence_tests(make_fixed_eos(), make_dynamic_eos(), make_states());

        "multiprecision derivative oracle and ADL math"_test = [] {
            const double point = 120.0;
            const double temperature = 350.0;
            const auto estimate = multiprecision_first_derivative(
                [temperature](const auto& concentration) {
                    using Number = std::remove_cvref_t<decltype(concentration)>;
                    const Number R{"8.31446261815324"};
                    return R * Number{temperature} * test_math::log(concentration);
                },
                point, 1.0, 10.0, 250.0);
            const auto actual = static_cast<double>(estimate.value);
            const double expected = fug::ideal_gas_constant<double> * temperature / point;
            const auto state = make_eos_test_state(point, std::array{0.4, 0.6}, temperature, "multiprecision-oracle");
            check_close("d(RT log(c))/dc", actual, expected, {.abs = 1e-13, .rel = 1e-13}, state,
                        static_cast<double>(estimate.step), static_cast<double>(estimate.error));

            const auto second = multiprecision_second_derivative(
                [temperature](const auto& concentration) {
                    using Number = std::remove_cvref_t<decltype(concentration)>;
                    return Number{"8.31446261815324"} * Number{temperature} * test_math::log(concentration);
                },
                point, 1.0, 10.0, 250.0);
            check_close("d2(RT log(c))/dc2", static_cast<double>(second.value), -expected / point,
                        {.abs = 1e-13, .rel = 1e-13}, state, static_cast<double>(second.step),
                        static_cast<double>(second.error));

            const auto mixed = multiprecision_mixed_derivative(
                [](const auto& concentration, const auto& T) {
                    return concentration * T + concentration * concentration * T * T;
                },
                2.0, 3.0, 1.0, 1.0, 1.0, 4.0, 1.0, 5.0);
            check_close("d2(cT+c2T2)/dcdT", static_cast<double>(mixed.value), 25.0, {.abs = 1e-20, .rel = 1e-20}, state,
                        static_cast<double>(mixed.step), static_cast<double>(mixed.error));
        };
    };

    suite<"automated_residual_contracts"> residual = [] {
        auto dynamic_eos = make_dynamic_vdw_eos();
        const eos_valid_domain domain{.c_min = 1e-8,
                                      .c_max = 5000.0,
                                      .T_min = 250.0,
                                      .T_max = 500.0,
                                      .minimum_mole_fraction = 0.01,
                                      .seed = 0xC0FFEE,
                                      .random_samples = 50};
        const std::vector<eos_test_state> states{
            {.c = 25.0, .x = {0.2, 0.8}, .T = 270.0, .effective_molar_mass = 0.036, .label = "dilute"},
            {.c = 500.0, .x = {0.4, 0.6}, .T = 330.0, .effective_molar_mass = 0.036, .label = "representative"},
            {.c = 4000.0, .x = {0.75, 0.25}, .T = 450.0, .effective_molar_mass = 0.036, .label = "dense-valid"}};
        const auto fixture = eos_test_fixture{
            .contribution = dynamic_eos.residual(), .eos = dynamic_eos, .states = states, .domain = domain};
        register_eos_contract_tests(fixture);
        register_residual_contract_tests(
            fixture, {.dilute_concentration = 1e-8, .dilute_tolerance = {.abs = 1e-6, .rel = 1e-6}});
        register_static_dynamic_equivalence_tests(make_fixed_vdw_eos(), make_dynamic_vdw_eos(), states);
    };

    suite<"automated_peng_robinson_contracts"> peng_robinson = [] {
        auto dynamic_eos = make_dynamic_pr_eos();
        const eos_valid_domain domain{.c_min = 1e-8,
                                      .c_max = 4000.0,
                                      .T_min = 240.0,
                                      .T_max = 600.0,
                                      .minimum_mole_fraction = 0.01,
                                      .seed = 0xC0FFEE,
                                      .random_samples = 50};
        const std::vector<eos_test_state> states{
            {.c = 25.0, .x = {0.2, 0.8}, .T = 270.0, .effective_molar_mass = 0.036, .label = "dilute"},
            {.c = 500.0, .x = {0.4, 0.6}, .T = 330.0, .effective_molar_mass = 0.036, .label = "representative"},
            {.c = 3000.0, .x = {0.75, 0.25}, .T = 500.0, .effective_molar_mass = 0.036, .label = "dense-valid"}};
        const auto fixture = eos_test_fixture{
            .contribution = dynamic_eos.residual(), .eos = dynamic_eos, .states = states, .domain = domain};
        register_eos_contract_tests(fixture);
        register_residual_contract_tests(
            fixture, {.dilute_concentration = 1e-8, .dilute_tolerance = {.abs = 1e-6, .rel = 1e-6}});
        register_static_dynamic_equivalence_tests(make_fixed_pr_eos(), make_dynamic_pr_eos(), states);
    };

    return ::boost::ut::cfg<>.run();
}
