#pragma once

#include <algorithm>
#include <boost/ut.hpp>
#include <cmath>
#include <format>
#include <string_view>

namespace synthesize_test {

// Relative-error expectation for model-specific reference-value checks. The
// scale floor makes the check useful for values close to zero as well.
inline void check_rel(std::string_view name, double actual, double expected, double reltol)
{
    using namespace boost::ut;
    const double scale = std::max(1.0, std::abs(expected));
    const double relative_error = std::abs(actual - expected) / scale;
    expect(relative_error <= reltol) << std::format("{}: actual={:.12g} expected={:.12g} rel_err={:.3e} (tol={:.1e})",
                                                    name, actual, expected, relative_error, reltol);
}

} // namespace synthesize_test
