#pragma once
///
/// File ``horner.hpp``.
/// Compile-time-sized polynomial evaluation.
///

#include "fugacity/core/attributes.hpp"

#include <array>
#include <concepts>

namespace fugacity {
///
/// Evaluate a degree-``N`` polynomial at ``x``.
///
/// Evaluates :math:`\sum_{k=0}^{N} \text{coeffs}[k]\,x^k` using Horner's method.
///
///
/// :tparam N: Degree of the polynomial (``coeffs`` holds ``N`` + 1 values).
/// :tparam Number: A floating-point type.
/// :param coeffs: Coefficients in ascending power order (``coeffs[k]`` multiplies :math:`x^k`).
/// :param x: The evaluation point.
/// :returns: The polynomial value.
///
/// .. code-block:: cpp
///
///    std::array<double, 3> coeffs{1.0, 2.0, 3.0}; // 1 + 2x + 3x^2
///    double y = fugacity::eval_polynomial<2>(coeffs, 2.0); // 17.0
///
///
/// \ingroup core
template<int N, std::floating_point Number>
[[nodiscard]] [[clang::always_inline]] constexpr Number eval_polynomial(std::array<Number, N + 1>& coeffs, Number x)
{
    Number p = coeffs[N];
    for (int idx = N - 1; idx >= 0; --idx) {
        p *= x;
        p += coeffs[idx];
    }
    return p;
}
} // namespace fugacity
