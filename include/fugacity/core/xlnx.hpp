#pragma once

#include "fugacity/core/core_calculations.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
namespace fugacity {
namespace detail {
//
// Internal implementation detail.
// Compile-time binomial coefficient :math:`\binom{n}{k}`.
//
// Evaluated with the multiplicative formula using the symmetry
//
// :math:`\binom{n}{k} = \binom{n}{n-k}` to keep the running product small.
//
//
// :param n: The size of the set.
// :param k: The number of elements chosen.
// :returns: :math:`\binom{n}{k}`, or :math:`0` when ``k`` > ``n``.
//
// \ingroup core
consteval unsigned long long choose(unsigned int n, unsigned int k)
{
    if (k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }

    k = std::min(k, n - k);

    unsigned long long result = 1;
    for (unsigned int i = 1; i <= k; ++i) {
        result *= (n - k + i);
        result /= i;
    }
    return result;
}

//
// Internal implementation detail.
// Coefficients of the inner polynomial of the order-``N`` smoothstep.
//
// The order-``N`` smoothstep is
//
// .. math::
//
//      S_N(x) = x^{N+1} \sum_{n=0}^{N} c_n\, x^{n},
//      \qquad c_n = (-1)^n \binom{N+n}{n}\binom{2N+1}{N-n},
//
// the unique degree-:math:`(2N+1)` polynomial with :math:`S_N(0)=0`, :math:`S_N(1)=1`
// and vanishing first ``N`` derivatives at both endpoints. This returns the
// coefficients :math:`c_0,\dots,c_N` of the inner sum (the leading :math:`x^{N+1}`
// factor is applied separately in detail::smooth_step). The coefficients
// alternate in sign, so a signed integer type is required.
//
//
// :tparam N: The smoothness order (the result is ``N`` + 1 coefficients).
// :returns: The coefficients :math:`c_0,\dots,c_N`, in ascending power order.
//
// \ingroup core
template<int N> consteval std::array<long long, N + 1> get_smooth_step_coeffs()
{
    std::array<long long, N + 1> coeffs{};
    for (int n = 0; n <= N; ++n) {
        long long coeff = n % 2 == 0 ? 1 : -1;
        coeff *= static_cast<long long>(choose(N + n, n));
        coeff *= static_cast<long long>(choose((2 * N) + 1, N - n));
        coeffs[n] = coeff;
    }
    return coeffs;
}

//
// Internal implementation detail.
// Generalized smoothstep :math:`S_N` of order ``N``.
//
// Evaluates :math:`S_N(x) = x^{N+1} \sum_{n=0}^{N} c_n x^n` (see
// detail::get_smooth_step_coeffs) using Horner's method for the inner sum and
// detail::fast_pow for the :math:`x^{N+1}` factor. On :math:`[0, 1]` it ramps
// monotonically from :math:`S_N(0) = 0` to :math:`S_N(1) = 1` with its first ``N``
// derivatives vanishing at both endpoints, so it is the :math:`C^N` "S-curve" used
// to smooth the otherwise non-differentiable kink of xlnx() at the origin.
//
//
// :tparam N: The smoothness order (number of vanishing endpoint
//                derivatives). Must be in :math:`[0, 6]`.
// :tparam Number: A floating-point type.
// :param x: The argument; the smoothstep behaviour is intended for
//                :math:`x \in [0, 1]`.
// :returns: :math:`S_N(x)`.
//
// \ingroup core
template<int N, std::floating_point Number> constexpr Number smooth_step(const Number x)
{
    static_assert(N >= 0, "The smoothness order N must be non-negative.");
    static_assert(N < 7, "Using a large value for N is almost certainly unnecessary.");
    constexpr auto coeffs = get_smooth_step_coeffs<N>();
    // Horner's method evaluation of the inner polynomial sum_{n} c_n x^n.
    auto result = static_cast<Number>(coeffs[N]);
    for (int idx = N - 1; idx >= 0; --idx) {
        result *= x;
        result += static_cast<Number>(coeffs[idx]);
    }
    // Apply the overall x^{N+1} factor carried by the smoothstep.
    result *= fast_pow<Number, N + 1>(x);
    return result;
}

} // namespace detail
///
/// Evaluate :math:`x\ln x` with value zero for nonpositive arguments.
///
/// This definition uses :math:`\lim_{x\to0^+}x\ln x=0` at the origin and does
/// not call ``std::log`` for :math:`x\le0`.
///
///
/// :tparam Number: A floating-point type.
/// :param x: The argument.
/// :returns: :math:`x \ln(x)` for :math:`x > 0`, otherwise :math:`0`.
///
/// \id exact
/// \ingroup core
template<std::floating_point Number> Number xlnx(const Number x)
{
    return x <= Number{0} ? Number{0} : x * std::log(x);
}

///
/// Evaluate a smooth extension of :math:`x\ln x` at the origin.
///
/// ``Continuity == 0`` is equivalent to the unsmoothed ``xlnx`` overload. For
/// ``Continuity >= 1``, a smoothstep modifies the function over
/// :math:`0<x<\varepsilon`, where :math:`\varepsilon` is machine epsilon. The
/// result has ``Continuity`` continuous derivatives and equals :math:`x\ln x`
/// for :math:`x\ge\varepsilon`.
///
///
/// :tparam Continuity: The order of continuity to enforce at the origin
///                    (number of continuous derivatives). Must be non-negative;
///                    see detail::smooth_step for the supported upper bound.
/// :tparam Number: A floating-point type.
/// :param x: The argument.
/// :returns: The :math:`C^{\text{Continuity}}` extension of :math:`x \ln(x)`.
///
/// \id smooth
/// \ingroup core
template<int Continuity, std::floating_point Number> Number xlnx(const Number x)
{
    static_assert(Continuity >= 0, "The continuity parameter must be non-negative!");
    if constexpr (Continuity == 0) {
        return xlnx(x);
    }
    else {
        constexpr auto eps = std::numeric_limits<Number>::epsilon();
        if (x <= Number{0}) {
            return Number{0};
        }
        if (x >= eps) {
            return x * std::log(x);
        }
        return detail::smooth_step<Continuity, Number>(x / eps) * x * std::log(x);
    }
}
} // namespace fugacity
