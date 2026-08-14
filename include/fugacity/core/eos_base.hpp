#pragma once

///
/// Component-count storage and the ideal-model tag.
///

#include "fugacity/core/assertions.hpp"
#include "fugacity/core/attributes.hpp"

#include <cstddef>
#include <span>
#include <utility>
namespace fugacity {

///
/// Store the number of components in an equation-of-state model.
///
/// Specify ``N`` when the component count is known at compile time. Use
/// ``std::dynamic_extent`` when the count is determined at run time.
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \id fixed-size
/// \ingroup core
template<std::size_t N = std::dynamic_extent> class BaseEoS {
public:
    constexpr BaseEoS() noexcept = default;

    ///
    /// Check an explicit component count against ``N``.
    ///
    /// :param n: Component count. Must equal ``N``. A mismatch throws
    ///           ``std::logic_error`` in debug builds; the check is omitted in
    ///           release builds.
    ///
    constexpr explicit BaseEoS([[maybe_unused]] const std::size_t n) { FUGACITY_ASSERT(n == N); }

    ///
    /// Number of chemical components.
    ///
    /// :returns: ``N``
    ///
    [[nodiscard]] static constexpr std::size_t size() noexcept { return N; }

    ///
    /// Call ``f(i)`` for each component index in ascending order.
    ///
    /// :tparam F: Type callable as ``f(std::size_t)``.
    /// :param f: Callable applied to the indices ``0`` through ``size() - 1``.
    ///
    /// .. code-block:: cpp
    ///
    ///    // Zero every entry in a component-sized buffer.
    ///    model.for_each_component([&](std::size_t i) { out[i] = 0.0; });
    ///
    ///
    template<class F> [[clang::always_inline]] constexpr void for_each_component(F f) const
    {
        [&]<std::size_t... Is>(std::index_sequence<Is...>) { (f(Is), ...); }(std::make_index_sequence<N>{});
    }
};

///
/// Store a component count that is determined at run time.
///
/// \id runtime-size
/// \ingroup core
template<> class BaseEoS<std::dynamic_extent> {
public:
    ///
    /// Construct with a component count.
    ///
    /// :param n: Number of components.
    ///
    constexpr explicit BaseEoS(const std::size_t n) noexcept : n_(n) {}

    ///
    /// Number of chemical components.
    ///
    /// :returns: The stored component count.
    ///
    [[nodiscard]] constexpr std::size_t size() const noexcept { return n_; }

    ///
    /// Call ``f(i)`` for each component index in ascending order.
    ///
    /// :tparam F: Type callable as ``f(std::size_t)``.
    /// :param f: Callable applied to the indices ``0`` through ``size() - 1``.
    ///
    template<class F> [[clang::always_inline]] constexpr void for_each_component(F f) const
    {
        for (std::size_t idx = 0; idx < n_; ++idx) {
            f(idx);
        }
    }

private:
    std::size_t n_; // Number of components.
};

///
/// Mark a model as an ideal-gas contribution.
///
/// Publicly derive an ideal model from this class so that it satisfies
/// :cpp:concept:`fugacity::IdealEoS`. Residual models must not derive from it.
///
/// \ingroup core
class BaseIdealEoS {};
} // namespace fugacity
