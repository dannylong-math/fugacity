#pragma once

/**
 * @file eos_pair.hpp
 * @brief Couples an ideal and a residual model into a single equation of state.
 */

#include "fugacity/core/assertions.hpp"
#include "fugacity/core/concepts.hpp"
namespace fugacity {

/**
 * @brief A complete equation of state, formed by pairing an ideal contribution
 *        with a residual (departure) contribution.
 *
 * Thermodynamic properties are obtained by passing an instance of this class to
 * the free functions in core_calculations.hpp. The total reduced Helmholtz
 * energy is the sum of the ideal and residual parts, and the property routines
 * combine the two contributions as the relevant thermodynamics dictates.
 *
 * @tparam Ideal    A model satisfying fugacity::IdealEoS.
 * @tparam Residual A model satisfying fugacity::ResidualEoS.
 */
template<IdealEoS Ideal, ResidualEoS Residual> class EoS {
public:
    using ideal_type = Ideal;       ///< Type of the ideal contribution.
    using residual_type = Residual; ///< Type of the residual contribution.

    /**
     * @brief Construct from an ideal and a residual model.
     * @param ideal    The ideal contribution.
     * @param residual The residual contribution.
     * @pre  `ideal.size() == residual.size()`: both models must describe the same
     *       number of components. In a debug build a mismatch throws
     *       @c std::logic_error (via FUGACITY_ASSERT); in a release build the
     *       check is elided, so the @c noexcept specification is preserved.
     */
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

    /// @brief Access the ideal contribution.
    const Ideal& ideal() const noexcept { return ideal_; }
    /// @brief Access the residual contribution.
    const Residual& residual() const noexcept { return residual_; }
    /**
     * @brief Number of chemical components.
     * @return Component count (equal for both contributions).
     */
    [[nodiscard]] constexpr std::size_t size() const noexcept { return ideal_.size(); }

private:
    Ideal ideal_;       ///< The ideal contribution.
    Residual residual_; ///< The residual contribution.
};
} // namespace fugacity