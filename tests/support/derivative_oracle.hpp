#pragma once

#include "eos_test_state.hpp"

#include <algorithm>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/ut.hpp>
#include <cmath>
#include <format>
#include <limits>
#include <string_view>
#include <utility>

namespace synthesize_test {

using multiprecision_float = boost::multiprecision::cpp_dec_float_50;

// Unqualified calls let ordinary overloads come from std and multiprecision
// overloads come from argument-dependent lookup.
namespace test_math {
template<class Number> auto log(const Number& value)
{
    using std::log;
    return log(value);
}
template<class Number> auto sqrt(const Number& value)
{
    using std::sqrt;
    return sqrt(value);
}
template<class Number> auto exp(const Number& value)
{
    using std::exp;
    return exp(value);
}
} // namespace test_math

struct tolerance {
    double abs{1e-8};
    double rel{1e-7};
};

inline bool close(double actual, double expected, tolerance limits)
{
    return std::abs(actual - expected) <= limits.abs + limits.rel * std::abs(expected);
}

inline void check_close(std::string_view property, double actual, double expected, tolerance limits,
                        const eos_test_state& state, double step = 0.0, double oracle_error = 0.0)
{
    using namespace boost::ut;
    const double absolute_error = std::abs(actual - expected);
    const double relative_error = absolute_error / std::max(std::abs(expected), limits.abs);
    expect(close(actual, expected, limits))
        << std::format("{}: {} actual={:.16g} expected={:.16g} abs_error={:.3e} rel_error={:.3e} "
                       "abs_tol={:.3e} rel_tol={:.3e} step={:.3e} oracle_error={:.3e}",
                       property, describe(state), actual, expected, absolute_error, relative_error, limits.abs,
                       limits.rel, step, oracle_error);
}

template<class Number> struct derivative_estimate {
    Number value{};
    Number step{};
    Number error{};
    bool stable{};
};

namespace detail {
template<class Number, class Function> Number first_stencil(Function&& function, Number point, Number step)
{
    const auto evaluate = [&](Number value) -> Number { return Number{function(value)}; };
    return (-evaluate(point + Number{2} * step) + Number{8} * evaluate(point + step) -
            Number{8} * evaluate(point - step) + evaluate(point - Number{2} * step)) /
           (Number{12} * step);
}

template<class Number, class Function> Number second_stencil(Function&& function, Number point, Number step)
{
    const auto evaluate = [&](Number value) -> Number { return Number{function(value)}; };
    return (-evaluate(point + Number{2} * step) + Number{16} * evaluate(point + step) - Number{30} * evaluate(point) +
            Number{16} * evaluate(point - step) - evaluate(point - Number{2} * step)) /
           (Number{12} * step * step);
}

template<class Number, class Function>
Number mixed_stencil(Function&& function, Number first, Number second, Number first_step, Number second_step)
{
    const auto evaluate = [&](Number first_value, Number second_value) -> Number {
        return Number{function(first_value, second_value)};
    };
    return (evaluate(first + first_step, second + second_step) - evaluate(first + first_step, second - second_step) -
            evaluate(first - first_step, second + second_step) + evaluate(first - first_step, second - second_step)) /
           (Number{4} * first_step * second_step);
}
} // namespace detail

template<class Number, class Function>
derivative_estimate<Number> adaptive_first_derivative(Function&& function, Number point, Number scale, Number lower,
                                                      Number upper)
{
    using std::abs;
    Number step_factor{};
    if constexpr (std::numeric_limits<Number>::digits10 > 30) {
        step_factor = Number{"1e-8"};
    }
    else {
        step_factor = Number{1e-3};
    }
    const Number magnitude{abs(point)};
    Number step = (magnitude > scale ? magnitude : scale) * step_factor;
    const Number lower_distance{(point - lower) / Number{2.1}};
    const Number upper_distance{(upper - point) / Number{2.1}};
    const Number boundary_step = lower_distance < upper_distance ? lower_distance : upper_distance;
    step = step < boundary_step ? step : boundary_step;
    derivative_estimate<Number> best{};
    best.error = (std::numeric_limits<Number>::max)();
    if (!(step > Number{0})) {
        return best;
    }

    for (int attempt = 0; attempt < 7; ++attempt) {
        const Number coarse = detail::first_stencil<Number>(function, point, step);
        const Number fine = detail::first_stencil<Number>(function, point, step / Number{2});
        const Number error = abs(fine - coarse) / Number{15};
        if (error < best.error) {
            best = {fine, step / Number{2}, error, true};
        }
        step /= Number{2};
    }
    return best;
}

template<class Number, class Function>
derivative_estimate<Number> adaptive_second_derivative(Function&& function, Number point, Number scale, Number lower,
                                                       Number upper)
{
    using std::abs;
    Number step_factor{};
    if constexpr (std::numeric_limits<Number>::digits10 > 30) {
        step_factor = Number{"1e-6"};
    }
    else {
        step_factor = Number{3e-3};
    }
    const Number magnitude{abs(point)};
    Number step = (magnitude > scale ? magnitude : scale) * step_factor;
    const Number lower_distance{(point - lower) / Number{2.1}};
    const Number upper_distance{(upper - point) / Number{2.1}};
    const Number boundary_step = lower_distance < upper_distance ? lower_distance : upper_distance;
    step = step < boundary_step ? step : boundary_step;
    derivative_estimate<Number> best{};
    best.error = (std::numeric_limits<Number>::max)();
    if (!(step > Number{0})) {
        return best;
    }

    for (int attempt = 0; attempt < 7; ++attempt) {
        const Number coarse = detail::second_stencil<Number>(function, point, step);
        const Number fine = detail::second_stencil<Number>(function, point, step / Number{2});
        const Number error = abs(fine - coarse) / Number{15};
        if (error < best.error) {
            best = {fine, step / Number{2}, error, true};
        }
        step /= Number{2};
    }
    return best;
}

template<class Number, class Function>
derivative_estimate<Number> adaptive_mixed_derivative(Function&& function, Number first, Number second,
                                                      Number first_scale, Number second_scale, Number first_lower,
                                                      Number first_upper, Number second_lower, Number second_upper)
{
    using std::abs;
    Number step_factor{};
    if constexpr (std::numeric_limits<Number>::digits10 > 30) {
        step_factor = Number{"1e-10"};
    }
    else {
        step_factor = Number{1e-3};
    }
    const Number first_magnitude{abs(first)};
    const Number second_magnitude{abs(second)};
    Number first_step = (first_magnitude > first_scale ? first_magnitude : first_scale) * step_factor;
    Number second_step = (second_magnitude > second_scale ? second_magnitude : second_scale) * step_factor;
    const Number first_boundary =
        ((first - first_lower) < (first_upper - first) ? (first - first_lower) : (first_upper - first)) / Number{1.1};
    const Number second_boundary =
        ((second - second_lower) < (second_upper - second) ? (second - second_lower) : (second_upper - second)) /
        Number{1.1};
    first_step = first_step < first_boundary ? first_step : first_boundary;
    second_step = second_step < second_boundary ? second_step : second_boundary;
    derivative_estimate<Number> best{};
    best.error = (std::numeric_limits<Number>::max)();
    if (!(first_step > Number{0} && second_step > Number{0})) {
        return best;
    }

    for (int attempt = 0; attempt < 7; ++attempt) {
        const Number coarse = detail::mixed_stencil<Number>(function, first, second, first_step, second_step);
        const Number fine =
            detail::mixed_stencil<Number>(function, first, second, first_step / Number{2}, second_step / Number{2});
        const Number error = abs(fine - coarse) / Number{3};
        if (error < best.error) {
            best = {fine, (first_step < second_step ? first_step : second_step) / Number{2}, error, true};
        }
        first_step /= Number{2};
        second_step /= Number{2};
    }
    return best;
}

template<class Function>
auto multiprecision_first_derivative(Function&& function, double point, double scale, double lower, double upper)
{
    return adaptive_first_derivative<multiprecision_float>(std::forward<Function>(function),
                                                           multiprecision_float{point}, multiprecision_float{scale},
                                                           multiprecision_float{lower}, multiprecision_float{upper});
}

template<class Function>
auto multiprecision_second_derivative(Function&& function, double point, double scale, double lower, double upper)
{
    return adaptive_second_derivative<multiprecision_float>(std::forward<Function>(function),
                                                            multiprecision_float{point}, multiprecision_float{scale},
                                                            multiprecision_float{lower}, multiprecision_float{upper});
}

template<class Function>
auto multiprecision_mixed_derivative(Function&& function, double first, double second, double first_scale,
                                     double second_scale, double first_lower, double first_upper, double second_lower,
                                     double second_upper)
{
    return adaptive_mixed_derivative<multiprecision_float>(
        std::forward<Function>(function), multiprecision_float{first}, multiprecision_float{second},
        multiprecision_float{first_scale}, multiprecision_float{second_scale}, multiprecision_float{first_lower},
        multiprecision_float{first_upper}, multiprecision_float{second_lower}, multiprecision_float{second_upper});
}

} // namespace synthesize_test
