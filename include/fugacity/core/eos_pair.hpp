#pragma once

///
/// Pair an ideal model with a residual model.
///

#include "fugacity/core/assertions.hpp"
#include "fugacity/core/concepts.hpp"
namespace fugacity {

///
/// Combine an ideal contribution and a residual contribution into a complete
/// equation of state.
///
/// For molar Helmholtz energy,
///
/// .. math::
///
///    a(c,\boldsymbol{x},T)
///    = a^\mathrm{ideal}(c,\boldsymbol{x},T)
///    + a^\mathrm{res}(c,\boldsymbol{x},T).
///
/// Pass the resulting object to the ``calc_*`` property functions.
///
/// :tparam Ideal: Ideal contribution satisfying :cpp:concept:`fugacity::IdealEoS`.
/// :tparam Residual: Residual contribution satisfying :cpp:concept:`fugacity::ResidualEoS`.
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual> class EoS {
public:
    using ideal_type = Ideal;       ///< Type of the ideal contribution.
    using residual_type = Residual; ///< Type of the residual contribution.

    ///
    /// Construct a complete equation of state.
    ///
    /// :param ideal: Ideal contribution.
    /// :param residual: Residual contribution.
    /// :precondition: ``ideal.size() == residual.size()``. A mismatch throws
    ///                ``std::logic_error`` in debug builds; the check is omitted
    ///                in release builds.
    ///
    EoS(Ideal ideal, Residual residual)
#ifdef NDEBUG
        // In release the precondition check is elided, so the constructor can
        // only throw if a move constructor does. In debug FUGACITY_ASSERT may
        // throw std::logic_error, so the constructor is left potentially-throwing.
        noexcept(std::is_nothrow_move_constructible_v<Ideal> && std::is_nothrow_move_constructible_v<Residual>)
#endif
        :
        ideal_{std::move(ideal)}, residual_{std::move(residual)}
    {
        FUGACITY_ASSERT(ideal_.size() == residual_.size());
    }

    /// Return the ideal contribution.
    const Ideal& ideal() const noexcept { return ideal_; }
    /// Return the residual contribution.
    const Residual& residual() const noexcept { return residual_; }
    ///
    /// Number of chemical components.
    ///
    /// :returns: Component count (equal for both contributions).
    ///
    [[nodiscard]] constexpr std::size_t size() const noexcept { return ideal_.size(); }

private:
    Ideal ideal_;       // The ideal contribution.
    Residual residual_; // The residual contribution.
};
} // namespace fugacity
