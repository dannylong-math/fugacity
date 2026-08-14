#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iterator>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace synthesize_test {

struct eos_test_state {
    double c{};
    std::vector<double> x;
    double T{};
    double effective_molar_mass{0.02};
    std::string label;
};

template<class X>
eos_test_state make_eos_test_state(double c, const X& x, double T, std::string_view label,
                                   double effective_molar_mass = 0.02)
{
    return {c, {std::begin(x), std::end(x)}, T, effective_molar_mass, std::string{label}};
}

struct eos_valid_domain {
    double c_min{};
    double c_max{};
    double T_min{};
    double T_max{};
    double minimum_mole_fraction{1e-4};
    std::uint32_t seed{0xC0FFEE};
    std::size_t random_samples{50};

    [[nodiscard]] bool contains(const eos_test_state& state) const
    {
        if (!std::isfinite(state.c) || !std::isfinite(state.T) || state.c < c_min || state.c > c_max ||
              state.T < T_min || state.T > T_max || !std::isfinite(state.effective_molar_mass) ||
              state.effective_molar_mass <= 0.0 ||
            state.x.empty()) {
            return false;
        }
        const double sum = std::accumulate(state.x.begin(), state.x.end(), 0.0);
        return std::ranges::all_of(
                   state.x, [&](double xi) { return std::isfinite(xi) && xi >= minimum_mole_fraction && xi <= 1.0; }) &&
               std::abs(sum - 1.0) <= 1e-12;
    }
};

inline std::string describe(const eos_test_state& state)
{
    std::string composition;
    for (std::size_t i = 0; i < state.x.size(); ++i) {
        composition += std::format("{}{}", i == 0 ? "" : ",", state.x[i]);
    }
    return std::format("state='{}' c={} T={} x=[{}] M={}", state.label, state.c, state.T, composition,
                       state.effective_molar_mass);
}

inline void validate_fixture_state(const eos_test_state& state, const eos_valid_domain& domain, std::size_t eos_size)
{
    if (state.x.size() != eos_size) {
        throw std::invalid_argument(
            std::format("{} has {} components, but the EoS has {}", describe(state), state.x.size(), eos_size));
    }
    if (!domain.contains(state)) {
        throw std::invalid_argument(
            std::format("fixture state lies outside its declared valid domain: {}", describe(state)));
    }
}

inline std::vector<eos_test_state> sample_valid_states(const eos_valid_domain& domain, std::size_t component_count,
                                                       double effective_molar_mass = 0.02)
{
    if (component_count == 0 || domain.c_min <= 0.0 || domain.c_max < domain.c_min || domain.T_min <= 0.0 ||
        domain.T_max < domain.T_min || domain.minimum_mole_fraction * static_cast<double>(component_count) >= 1.0) {
        throw std::invalid_argument("invalid EoS test domain");
    }

    std::mt19937 generator(domain.seed);
    std::uniform_real_distribution<double> c_distribution(domain.c_min, domain.c_max);
    std::uniform_real_distribution<double> T_distribution(domain.T_min, domain.T_max);
    std::exponential_distribution<double> composition_distribution(1.0);
    std::vector<eos_test_state> result;
    result.reserve(domain.random_samples);

    for (std::size_t sample = 0; sample < domain.random_samples; ++sample) {
        std::vector<double> weights(component_count);
        for (double& value : weights) {
            value = composition_distribution(generator);
        }
        const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
        const double remainder = 1.0 - domain.minimum_mole_fraction * static_cast<double>(component_count);
        for (double& value : weights) {
            value = domain.minimum_mole_fraction + remainder * value / total;
        }
        result.push_back({.c=c_distribution(generator), .x=std::move(weights), .T=T_distribution(generator),
                          .effective_molar_mass=effective_molar_mass, .label=std::format("random(seed={},sample={})", domain.seed, sample)});
    }
    return result;
}

} // namespace synthesize_test
