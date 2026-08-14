//
// Unit tests for fugacity::NoResidual, an identically zero residual model.
//
// Universal identities, derivatives, wrappers, preconditions, ideal-gas laws,
// deterministic sampling, and fixed/dynamic equivalence are registered through
// support/eos_test_suite.hpp. The explicit tests below cover the model's own
// zero-output and extent behavior.
//
#include "support/eos_test_suite.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/residual_models/no_residual.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cstddef>
#include <span>
#include <tuple>
#include <vector>

using namespace boost::ut;
using namespace fugacity_test;
using fugacity::NoResidual;

namespace {

namespace fug = fugacity;

using FixedIdeal = fug::ConstantCp<2>;
using DynamicIdeal = fug::ConstantCp<std::dynamic_extent>;

constexpr std::array<FixedIdeal::SpeciesInput, 2> ideal_inputs{{
    {.T_ref = 298.15, .p_ref = 1.0e5, .c_p = 29.1, .h_ref = 0.0, .s_ref = 191.0},
    {.T_ref = 310.0, .p_ref = 9.0e4, .c_p = 33.6, .h_ref = 1200.0, .s_ref = 205.0},
}};

std::vector<DynamicIdeal::SpeciesInput> dynamic_ideal_inputs()
{
    std::vector<DynamicIdeal::SpeciesInput> result;
    result.reserve(ideal_inputs.size());
    for (const auto& input : ideal_inputs) {
        result.push_back(
            {.T_ref = input.T_ref, .p_ref = input.p_ref, .c_p = input.c_p, .h_ref = input.h_ref, .s_ref = input.s_ref});
    }
    return result;
}

auto make_fixed_eos() { return fug::EoS{FixedIdeal{ideal_inputs}, NoResidual<2>{}}; }

auto make_dynamic_eos()
{
    const auto inputs = dynamic_ideal_inputs();
    return fug::EoS{DynamicIdeal{std::span<const DynamicIdeal::SpeciesInput>{inputs}},
                   NoResidual<std::dynamic_extent>{inputs.size()}};
}

std::vector<eos_test_state> contract_states()
{
    return {{.c = 1.0, .x = {0.2, 0.8}, .T = 220.0, .effective_molar_mass = 0.031, .label = "dilute"},
            {.c = 100.0, .x = {0.4, 0.6}, .T = 300.0, .effective_molar_mass = 0.031, .label = "representative"},
            {.c = 5000.0, .x = {0.7, 0.3}, .T = 500.0, .effective_molar_mass = 0.031, .label = "dense"}};
}

constexpr eos_valid_domain valid_domain{.c_min = 0.5,
                                        .c_max = 5500.0,
                                        .T_min = 200.0,
                                        .T_max = 550.0,
                                        .minimum_mole_fraction = 0.01,
                                        .seed = 0xC0FFEE,
                                        .random_samples = 50};

} // namespace

int main()
{
    suite<"no_residual"> s = [] {
        auto dynamic_eos = make_dynamic_eos();
        const auto fixture = eos_test_fixture{.contribution = dynamic_eos.residual(),
                                              .eos = dynamic_eos,
                                              .states = contract_states(),
                                              .domain = valid_domain};
        register_eos_contract_tests(fixture);
        register_ideal_gas_contract_tests(fixture);
        register_static_dynamic_equivalence_tests(make_fixed_eos(), make_dynamic_eos(), contract_states());

        // -------------------------------------------------------------------
        // Compile-time known size: default-constructible, fixed size().
        // -------------------------------------------------------------------
        "compile-time size: construction and zero output"_test = []<typename Number> {
            constexpr std::size_t N = 3;
            const NoResidual<N> res{};            // only a default ctor is needed
            expect(eq(NoResidual<N>::size(), N)); // size() is static for a compile-time size

            const std::array<Number, N> rho{Number{1.0}, Number{2.0}, Number{3.0}};
            const Number c{6.0};
            const Number T{300.0};

            expect(eq(res.calc_helmholtz(c, rho.data(), T), Number{0}));
            expect(eq(res.calc_helmholtz_density(rho.data(), T), Number{0}));

            // Pre-fill with garbage to confirm every entry is overwritten.
            std::array<Number, N> out{Number{7}, Number{-8}, Number{9}};
            res.calc_partial_helmholtz(rho.data(), T, out.data());
            for (std::size_t i = 0; i < N; ++i) {
                expect(eq(out[i], Number{0}));
            }
        } | std::tuple<float, double, long double>{};

        // -------------------------------------------------------------------
        // Runtime known size: size passed to the constructor, then honoured by
        // size() and by calc_partial_helmholtz's loop bound.
        // -------------------------------------------------------------------
        "runtime size: construction and zero output"_test = []<typename Number> {
            constexpr std::size_t N = 4;
            const NoResidual<std::dynamic_extent> res{N}; // runtime size to ctor
            expect(eq(res.size(), N));

            const std::array<Number, N> rho{Number{1.0}, Number{2.0}, Number{3.0}, Number{4.0}};
            const Number c{10.0};
            const Number T{275.0};

            expect(eq(res.calc_helmholtz(c, rho.data(), T), Number{0}));
            expect(eq(res.calc_helmholtz_density(rho.data(), T), Number{0}));

            std::array<Number, N> out{Number{7}, Number{-8}, Number{9}, Number{-1}};
            res.calc_partial_helmholtz(rho.data(), T, out.data());
            for (std::size_t i = 0; i < N; ++i) {
                expect(eq(out[i], Number{0}));
            }
        } | std::tuple<float, double, long double>{};
    };

    return ::boost::ut::cfg<>.run();
}
