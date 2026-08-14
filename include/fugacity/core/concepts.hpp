#pragma once
///
/// File ``concepts.hpp``.
/// Concepts that define the interface every equation-of-state model must
///        provide, and the ideal/residual classification.
///
#include "fugacity/core/eos_base.hpp"

#include <concepts>
#include <cstddef>
namespace fugacity {

///
/// Macro ``FUGACITY_RESTRICT``.
/// Portable spelling of the C99 ``restrict`` pointer qualifier.
///
/// Expands to the compiler's restrict keyword (``__restrict__`` on GCC/Clang,
/// ``__restrict`` on MSVC) or to nothing on unknown compilers. Used to promise that
/// the buffers passed to the hot autodiff routines do not alias.
///
#if defined(__GNUC__) || defined(__clang__)
#define FUGACITY_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define FUGACITY_RESTRICT __restrict
#else
#define FUGACITY_RESTRICT
#endif

///
/// Interface required of every equation-of-state model (ideal or residual).
///
/// A conforming type ``E`` must expose the following const member functions. The
/// concept is checked with ``double``, but models are expected to template these on
/// the floating-point type so the autodiff layer can also instantiate them with,
/// e.g., a higher-precision type.
///
/// Required members:
///
/// - ``calc_helmholtz(c, x, T) -> (convertible to double)``
///   Molar Helmholtz energy :math:`a`.
///   - ``c`` molar concentration [mol/m^3]
///   - ``x`` pointer to the mole-fraction array [-]
///   - ``T`` temperature [K]
///   - returns the **molar** Helmholtz energy [J/mol]
///
/// - ``calc_partial_helmholtz(rho_i, T, out) -> void``
///   Per-component decomposition of the Helmholtz energy **density** :math:`\Psi`
///   (such that :math:`\sum_i \text{out}[i] = \Psi`).
///   - ``rho_i`` pointer to the partial-molar-concentration array [mol/m^3]
///   - ``T`` temperature [K]
///   - ``out`` output array, per-component Helmholtz energy density [J/m^3]
///
/// - ``calc_helmholtz_density(rho_i, T) -> (convertible to double)``
///   Total Helmholtz energy **density** :math:`\Psi = \sum_i \Psi_i`, accumulated
///   directly into a scalar (no per-component buffer). This is the scalar the
///   reverse-mode autodiff differentiates to obtain chemical potentials.
///   - ``rho_i`` pointer to the partial-molar-concentration array [mol/m^3]
///   - ``T`` temperature [K]
///   - returns the total Helmholtz energy density [J/m^3]
///
/// - ``size() -> (convertible to std::size_t)``
///   Number of components [-].
///
/// **Note:** The two Helmholtz members use different units on purpose:
///       ``calc_helmholtz`` is an intensive **molar** quantity [J/mol], while
///       ``calc_partial_helmholtz`` returns a **volumetric** density [J/m^3]. They
///       are linked by :math:`\Psi(\rho,T) = c\,a(c,x,T)` with :math:`c=\sum_i\rho_i`
///       and :math:`x_i=\rho_i/c`.
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
/// An ideal-gas contribution: an EquationOfState that also derives from
///        BaseIdealEoS.
///
/// The inheritance requirement is how the library tells ideal models apart from
/// residual ones at compile time.
///
/// \ingroup core
template<class E>
concept IdealEoS = std::derived_from<E, BaseIdealEoS> && EquationOfState<E>;

///
/// A residual (departure) contribution: an EquationOfState that is *not*
///        an IdealEoS.
///
/// Equivalently, a model that does not derive from BaseIdealEoS.
///
/// \ingroup core
template<class E>
concept ResidualEoS = EquationOfState<E> && !IdealEoS<E>;

} // namespace fugacity
