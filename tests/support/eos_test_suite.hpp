#pragma once

#include "eos_contracts.hpp"

#include <boost/ut.hpp>
#include <utility>
#include <vector>

namespace synthesize_test {

template<class Contribution, class EoSPair> struct eos_test_fixture {
    Contribution contribution;
    EoSPair eos;
    std::vector<eos_test_state> states;
    eos_valid_domain domain;
};

template<class Fixture> void validate_fixture(const Fixture& fixture)
{
    if (fixture.contribution.size() != fixture.eos.size()) {
        throw std::invalid_argument("fixture contribution and complete EoS have different component counts");
    }
    if (fixture.states.empty()) {
        throw std::invalid_argument("an EoS test fixture needs at least one curated state");
    }
    for (const auto& state : fixture.states) {
        validate_fixture_state(state, fixture.domain, fixture.eos.size());
    }
}

template<class Fixture> void register_eos_contract_tests(Fixture fixture)
{
    using namespace boost::ut;
    validate_fixture(fixture);

    "contribution contracts (curated states)"_test = [fixture] {
        for (const auto& state : fixture.states) {
            check_contribution_contracts(fixture.contribution, state);
        }
    };
    "complete EoS identities (curated states)"_test = [fixture] {
        for (const auto& state : fixture.states) {
            check_complete_eos_identities(fixture.eos, state);
            check_property_wrapper_contracts(fixture.eos, state);
        }
    };
    "public API preconditions"_test = [fixture] { check_public_preconditions(fixture.eos, fixture.states.front()); };
    "all public first derivatives (adaptive finite differences)"_test = [fixture] {
        for (const auto& state : fixture.states) {
            check_public_first_derivatives(fixture.eos, state, fixture.domain);
        }
    };
    "deterministic valid-domain sweep"_test = [fixture] {
        const double molar_mass = fixture.states.front().effective_molar_mass;
        for (const auto& state : sample_valid_states(fixture.domain, fixture.eos.size(), molar_mass)) {
            validate_fixture_state(state, fixture.domain, fixture.eos.size());
            check_contribution_contracts(fixture.contribution, state);
            check_complete_eos_identities(fixture.eos, state);
            check_property_wrapper_contracts(fixture.eos, state);
        }
    };
}

template<class Fixture> void register_ideal_gas_contract_tests(Fixture fixture)
{
    using namespace boost::ut;
    validate_fixture(fixture);
    "ideal-gas contracts"_test = [fixture] {
        for (const auto& state : fixture.states) {
            check_ideal_gas_contracts(fixture.eos, state, fixture.domain);
        }
    };
}

template<class Fixture> void register_residual_contract_tests(Fixture fixture, residual_contract_options options)
{
    using namespace boost::ut;
    validate_fixture(fixture);
    if (!(options.dilute_concentration >= fixture.domain.c_min &&
          options.dilute_concentration <= fixture.domain.c_max)) {
        throw std::invalid_argument("dilute-limit concentration is outside the fixture's declared valid domain");
    }
    "residual dilute-limit contracts"_test = [fixture, options] {
        for (const auto& state : fixture.states) {
            check_residual_dilute_limit(fixture.eos, state, options);
        }
    };
}

template<class StaticEoS, class DynamicEoS>
void register_static_dynamic_equivalence_tests(StaticEoS fixed, DynamicEoS dynamic, std::vector<eos_test_state> states)
{
    using namespace boost::ut;
    "static/dynamic extent equivalence"_test = [fixed = std::move(fixed), dynamic = std::move(dynamic),
                                                states = std::move(states)] {
        for (const auto& state : states) {
            check_static_dynamic_equivalence(fixed, dynamic, state);
        }
    };
}

} // namespace synthesize_test
