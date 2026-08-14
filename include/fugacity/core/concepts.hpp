#pragma once
///
/// Equation-of-state model concepts.
///
#include "fugacity/core/eos_base.hpp"

#include <concepts>
#include <cstddef>
namespace fugacity {

///
/// Portable spelling of a restricted pointer qualifier.
///
/// The macro expands to the compiler-specific ``restrict`` keyword when one is
/// available and otherwise expands to nothing.
///
#if defined(__GNUC__) || defined(__clang__)
#define FUGACITY_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define FUGACITY_RESTRICT __restrict
#else
#define FUGACITY_RESTRICT
#endif

///
/// Require the common interface for an ideal or residual model.
///
/// Implement the three calculation functions as templates over a floating-point
/// number type. The concept checks the following interface with ``double``:
///
/// Required members:
///
/// - ``calc_helmholtz(c, x, T)`` returns molar Helmholtz energy :math:`a`
///   [J/mol]. Here ``c`` is molar concentration [mol/m^3], ``x`` points to the
///   mole fractions [-], and ``T`` is temperature [K].
///
/// - ``calc_helmholtz_density(rho_i, T)`` returns Helmholtz energy density
///   :math:`\Psi` [J/m^3]. ``rho_i`` points to partial molar concentrations
///   [mol/m^3].
///
/// - ``calc_partial_helmholtz(rho_i, T, out)`` writes a per-component
///   decomposition :math:`\Psi_i` [J/m^3] satisfying
///   :math:`\sum_i\Psi_i=\Psi`.
///
/// - ``size()`` returns the component count.
///
/// Keep the molar and density forms consistent:
///
/// .. math::
///
///    \Psi(\boldsymbol{\rho},T)=c\,a(c,\boldsymbol{x},T),\qquad
///    c=\sum_i\rho_i,\qquad x_i=\rho_i/c.
///
/// \ingroup core
template<class E>
concept EquationOfState =
    requires(const E& eos, const double mole_concentration, const double* mole_fractions,
             const double* partial_mole_concentrations, const double temperature, double* out_array) {
        { eos.calc_partial_helmholtz(partial_mole_concentrations, temperature, out_array) } -> std::same_as<void>;
        { eos.calc_helmholtz_density(partial_mole_concentrations, temperature) } -> std::convertible_to<double>;
        { eos.calc_helmholtz(mole_concentration, mole_fractions, temperature) } -> std::convertible_to<double>;
        { eos.size() } -> std::convertible_to<std::size_t>;
    };

///
/// Match an equation-of-state model marked as an ideal contribution.
///
/// \ingroup core
template<class E>
concept IdealEoS = std::derived_from<E, BaseIdealEoS> && EquationOfState<E>;

///
/// Match an equation-of-state model that is not marked as ideal.
///
/// \ingroup core
template<class E>
concept ResidualEoS = EquationOfState<E> && !IdealEoS<E>;

} // namespace fugacity
