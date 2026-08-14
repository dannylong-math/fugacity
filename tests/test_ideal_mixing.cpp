#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/ideal_models/nasa7.hpp"
#include "fugacity/ideal_models/nasa9.hpp"
#include "fugacity/residual_models/no_residual.hpp"
#include "support/numeric_checks.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <span>

using namespace boost::ut;
using namespace fugacity_test;

namespace {
namespace fug = fugacity;

template<class UnaryIdeal, class BinaryIdeal>
void check_ideal_mixing(const UnaryIdeal& unary_ideal, const BinaryIdeal& binary_ideal)
{
    const fug::EoS unary{unary_ideal, fug::NoResidual<1>{}};
    const fug::EoS binary{binary_ideal, fug::NoResidual<2>{}};

    constexpr double c = 120.0;
    constexpr double T = 375.0;
    constexpr double R = fug::ideal_gas_constant<double>;
    constexpr std::array<double, 1> unary_x{1.0};
    constexpr std::array<double, 2> binary_x{0.25, 0.75};
    constexpr std::array<double, 1> unary_rho{c};
    constexpr std::array<double, 2> binary_rho{c * binary_x[0], c * binary_x[1]};

    const double sum_xlnx = (binary_x[0] * std::log(binary_x[0])) + (binary_x[1] * std::log(binary_x[1]));
    const double expected_molar_mixing = R * T * sum_xlnx;
    const double expected_density_mixing = c * expected_molar_mixing;

    const double unary_a = unary_ideal.calc_helmholtz(c, unary_x.data(), T);
    const double binary_a = binary_ideal.calc_helmholtz(c, binary_x.data(), T);
    check_rel("molar Helmholtz mixing", binary_a - unary_a, expected_molar_mixing, 1e-11);

    const double unary_psi = unary_ideal.calc_helmholtz_density(unary_rho.data(), T);
    const double binary_psi = binary_ideal.calc_helmholtz_density(binary_rho.data(), T);
    check_rel("Helmholtz-density mixing", binary_psi - unary_psi, expected_density_mixing, 1e-11);

    std::array<double, 1> unary_partial{};
    std::array<double, 2> binary_partial{};
    unary_ideal.calc_partial_helmholtz(unary_rho.data(), T, unary_partial.data());
    binary_ideal.calc_partial_helmholtz(binary_rho.data(), T, binary_partial.data());
    check_rel("partial Helmholtz-density mixing", binary_partial[0] + binary_partial[1] - unary_partial[0],
              expected_density_mixing, 1e-11);

    const std::span<const double, 1> unary_x_span{unary_x};
    const std::span<const double, 2> binary_x_span{binary_x};
    check_rel("entropy of mixing",
              fug::calc_entropy(binary, c, binary_x_span, T) - fug::calc_entropy(unary, c, unary_x_span, T),
              -R * sum_xlnx, 1e-10);
    check_rel("Gibbs energy of mixing",
              fug::calc_gibbs(binary, c, binary_x_span, T) - fug::calc_gibbs(unary, c, unary_x_span, T),
              expected_molar_mixing, 1e-11);

    std::array<double, 2> chemical_potential{};
    fug::calc_chemical_potential(binary, std::span<const double, 2>{binary_rho}, T,
                                 std::span<double, 2>{chemical_potential});
    check_rel("ideal chemical-potential difference", chemical_potential[0] - chemical_potential[1],
              R * T * std::log(binary_x[0] / binary_x[1]), 1e-10);

    constexpr std::array<double, 2> pure_binary_x{1.0, 0.0};
    const double pure_binary_a = binary_ideal.calc_helmholtz(c, pure_binary_x.data(), T);
    expect(std::isfinite(pure_binary_a)) << "zero-composition Helmholtz contribution must be finite";
    check_rel("binary pure limit", pure_binary_a, unary_a, 1e-12);
}

template<class UnaryIdeal, class BinaryIdeal, class UnaryInput, class BinaryInput>
void run_model_test(const UnaryInput& unary_input, const BinaryInput& binary_input)
{
    check_ideal_mixing(UnaryIdeal{unary_input}, BinaryIdeal{binary_input});
}
} // namespace

int main()
{
    suite<"constant-cp ideal mixing"> constant_cp = [] {
        "physical mixing terms"_test = [] {
            using Unary = fug::ConstantCp<1>;
            using Binary = fug::ConstantCp<2>;
            constexpr Unary::SpeciesInput species{
                .T_ref = 300.0, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 1500.0, .s_ref = 191.0};
            constexpr std::array<Unary::SpeciesInput, 1> unary_input{species};
            constexpr std::array<Binary::SpeciesInput, 2> binary_input{{
                {.T_ref = species.T_ref,
                 .p_ref = species.p_ref,
                 .c_p = species.c_p,
                 .h_ref = species.h_ref,
                 .s_ref = species.s_ref},
                {.T_ref = species.T_ref,
                 .p_ref = species.p_ref,
                 .c_p = species.c_p,
                 .h_ref = species.h_ref,
                 .s_ref = species.s_ref},
            }};
            run_model_test<Unary, Binary>(unary_input, binary_input);
        };
    };

    suite<"NASA7 ideal mixing"> nasa7 = [] {
        "physical mixing terms"_test = [] {
            using Unary = fug::Nasa7<1>;
            using Binary = fug::Nasa7<2>;
            constexpr Unary::SpeciesInput species{.a0 = 3.53100528,
                                                  .a1 = -1.23660988e-4,
                                                  .a2 = -5.02999433e-7,
                                                  .a3 = 2.43530612e-9,
                                                  .a4 = -1.40881235e-12,
                                                  .a5 = -1046.97628,
                                                  .a6 = 2.96747038,
                                                  .T_ref = 300.0,
                                                  .p_ref = 1.0e5};
            constexpr std::array<Unary::SpeciesInput, 1> unary_input{species};
            constexpr Binary::SpeciesInput binary_species{.a0 = species.a0,
                                                          .a1 = species.a1,
                                                          .a2 = species.a2,
                                                          .a3 = species.a3,
                                                          .a4 = species.a4,
                                                          .a5 = species.a5,
                                                          .a6 = species.a6,
                                                          .T_ref = species.T_ref,
                                                          .p_ref = species.p_ref};
            constexpr std::array<Binary::SpeciesInput, 2> binary_input{binary_species, binary_species};
            run_model_test<Unary, Binary>(unary_input, binary_input);
        };
    };

    suite<"NASA9 ideal mixing"> nasa9 = [] {
        "physical mixing terms"_test = [] {
            using Unary = fug::Nasa9<1>;
            using Binary = fug::Nasa9<2>;
            constexpr Unary::SpeciesInput species{.a0 = 2.210371497e4,
                                                  .a1 = -3.818461820e2,
                                                  .a2 = 6.082738360,
                                                  .a3 = -8.530914410e-3,
                                                  .a4 = 1.384646189e-5,
                                                  .a5 = -9.625793620e-9,
                                                  .a6 = 2.519705809e-12,
                                                  .a7 = 7.108460860e2,
                                                  .a8 = -1.076003744e1,
                                                  .T_ref = 300.0,
                                                  .p_ref = 1.0e5};
            constexpr std::array<Unary::SpeciesInput, 1> unary_input{species};
            constexpr Binary::SpeciesInput binary_species{.a0 = species.a0,
                                                          .a1 = species.a1,
                                                          .a2 = species.a2,
                                                          .a3 = species.a3,
                                                          .a4 = species.a4,
                                                          .a5 = species.a5,
                                                          .a6 = species.a6,
                                                          .a7 = species.a7,
                                                          .a8 = species.a8,
                                                          .T_ref = species.T_ref,
                                                          .p_ref = species.p_ref};
            constexpr std::array<Binary::SpeciesInput, 2> binary_input{binary_species, binary_species};
            run_model_test<Unary, Binary>(unary_input, binary_input);
        };
    };

    return ::boost::ut::cfg<>.run();
}
