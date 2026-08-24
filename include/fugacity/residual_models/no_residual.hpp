#pragma once

#include "fugacity/core/concepts.hpp"

#include <span>

namespace fugacity {

///
/// Zero residual model.
///
/// The residual molar Helmholtz energy and Helmholtz energy density are
///
/// \f[
///
/// a^\mathrm{res}(c,\boldsymbol{x},T)=0,\qquad
/// \Psi^\mathrm{res}(\boldsymbol{\rho},T)=0.
///
/// \f]
/// \code{.cpp}
///
/// fugacity::NoResidual<2> residual;
/// const std::array<double, 2> x{0.4, 0.6};
/// const double a_res = residual.calc_helmholtz(40.0, x.data(), 300.0);
///
/// \endcode
/// \tparam N Component count, or `std::dynamic_extent` for a runtime count.
///
/// \ingroup residual-models
template<std::size_t N> class NoResidual : public BaseEoS<N> {
public:
    ///
    /// Construct a model whose component count is known at compile time.
    ///
    constexpr NoResidual() noexcept
        requires(N != std::dynamic_extent)
    = default;

    ///
    /// Construct a model with a runtime component count.
    ///
    /// \param n Component count.
    ///
    constexpr explicit NoResidual(const std::size_t n) noexcept
        requires(N == std::dynamic_extent)
        : BaseEoS<N>(n)
    {
    }

    ///
    /// Return the molar residual Helmholtz energy [J/mol].
    ///
    /// \return `Number{0}`.
    ///
    template<std::floating_point Number>
    Number calc_helmholtz(const Number /*c*/, const Number* /*x*/, const Number /*T*/) const
    {
        return Number{0};
    }

    ///
    /// Return the residual Helmholtz energy density [J/m^3].
    ///
    /// \return `Number{0}`.
    ///
    template<std::floating_point Number>
    Number calc_helmholtz_density(const Number* /*rho_i*/, const Number /*T*/) const
    {
        return Number{0};
    }

    ///
    /// Set every per-component residual Helmholtz energy density to zero.
    ///
    /// \param out Output array of length `size()` [J/m^3].
    ///
    template<std::floating_point Number>
    void calc_partial_helmholtz(const Number* /*rho_i*/, const Number /*T*/, Number* out) const
    {
        this->for_each_component([&](std::size_t idx) { out[idx] = Number{0}; });
    }
};

static_assert(ResidualEoS<NoResidual<2>>, "NoResidual must satisfy the ResidualEoS concept.");
static_assert(ResidualEoS<NoResidual<std::dynamic_extent>>, "NoResidual must satisfy the ResidualEoS concept.");

} // namespace fugacity
