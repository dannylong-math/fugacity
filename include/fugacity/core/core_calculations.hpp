#pragma once
///
/// File ``core_calculations.hpp``.
/// Thermodynamic property calculations built on top of an EoS pair.
///
/// Every property is derived from the reduced molar Helmholtz energy
///
/// :math:`\alpha = a/(RT)` and its derivatives. The derivatives are obtained by
/// automatic differentiation with `Enzyme <https://enzyme.mit.edu>`_:
///
/// - forward mode for the :math:`1/T`- and :math:`c`-derivatives (see
///   detail::calc_alpha / detail::calc_lambda), and
///
/// - reverse mode for the partial-molar derivatives w.r.t. each :math:`\rho_i`
///   (see detail::calc_dPsi_drhoi), used for chemical potentials and fugacities.
///
/// Symbol / unit conventions used throughout:
///
/// - ``c``     molar concentration (molar density) [mol/m^3]
/// - ``x``     mole fractions [-]
/// - ``T``     temperature [K]
/// - ``invT``  inverse temperature :math:`1/T` [1/K]
/// - ``rho_i`` partial molar concentrations [mol/m^3]
/// - ``R``     gas constant [J/(mol K)]
///
#include "fugacity/core/assertions.hpp"
#include "fugacity/core/concepts.hpp"
#include "fugacity/core/eos_pair.hpp"
#include "fugacity/core/numbers.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>

// NOLINTBEGIN
// Enzyme autodiff requires a few global definitions
inline int enzyme_dup;
inline int enzyme_dupnoneed;
inline int enzyme_out;
inline int enzyme_const;

template<typename return_type, typename... T> return_type __enzyme_fwddiff(void*, T...);

template<typename return_type, typename... T> return_type __enzyme_autodiff(void*, T...);
// NOLINTEND

namespace fugacity {

namespace detail {
//
// Internal implementation detail.
// Compile-time integer power :math:`\text{base}^N` by exponentiation by
//        squaring.
//
// :tparam Number: Floating-point base type.
// :tparam N: Integer exponent. May be negative (computes the reciprocal
//                power) or zero (returns 1).
// :param base: The base. Units are arbitrary; the result carries ``base``'s unit
//                raised to the ``N``.
// :returns: :math:`\text{base}^N`.
// :precondition: For ``N`` < 0, ``base != 0``.
//
// \ingroup core
template<std::floating_point Number, int N> constexpr Number fast_pow(const Number base)
{
    if constexpr (N == 0) {
        return Number{1};
    }
    else if constexpr (N < 0) {
        FUGACITY_ASSERT(base != Number{0});
        return fast_pow<Number, -N>(Number{1} / base);
    }
    else if constexpr (N == 1) {
        return base;
    }
    else if constexpr (N % 2 == 0) {
        Number half = fast_pow<Number, N / 2>(base);
        return half * half;
    }
    else {
        return base * fast_pow<Number, N - 1>(base);
    }
}

//
// Internal implementation detail.
// Mixed derivative of the reduced molar Helmholtz energy
//        :math:`\alpha = a/(RT)`.
//
// Returns
// :math:`\dfrac{\partial^{\,i+j}\alpha}{\partial (1/T)^{i}\,\partial c^{j}}`,
// evaluated by recursively applying Enzyme forward-mode differentiation: an
// :math:`1/T`-derivative while ``i > j``, otherwise a ``c``-derivative. The base case
// (``i == j == 0``) is :math:`\alpha = a/(RT)` itself.
//
// :tparam i: Number of derivatives w.r.t. inverse temperature :math:`1/T`.
// :tparam j: Number of derivatives w.r.t. concentration ``c``.
// :param eos: A single-contribution model (ideal or residual).
// :param c: Molar concentration [mol/m^3].
// :param x: Mole-fraction array [-].
// :param invT: Inverse temperature :math:`1/T` [1/K].
// :returns: The derivative; :math:`\alpha` is dimensionless, so the result has units
//         :math:`\mathrm{K}^{i}\,(\mathrm{m^3/mol})^{j}`.
//
// \ingroup core
template<int i, int j, EquationOfState EoS, std::floating_point Number>
[[nodiscard]] Number calc_alpha(const EoS& eos, const Number c, const Number* x, const Number invT)
{
    static_assert(i >= 0, "The template parameter `i` must be non-negative. It represent the number of derivatives "
                          "with respect to `invT`.");
    static_assert(
        j >= 0,
        "The template parameter `j` must be non-negative. It represent the number of derivatives with respect to `c`.");
    if constexpr (i == 0 && j == 0) {
        const Number T = Number{1} / invT;
        constexpr Number R = ideal_gas_constant<Number>;
        return eos.calc_helmholtz(c, x, T) / (R * T);
    }
    else if constexpr (i > j) {
        Number dinvT{1.};
        return __enzyme_fwddiff<Number>((void*)calc_alpha<i - 1, j, EoS, Number>, enzyme_const, &eos, enzyme_const, c,
                                        enzyme_const, x, enzyme_dup, invT, dinvT);
    }
    else {
        Number dc{1.};
        return __enzyme_fwddiff<Number>((void*)calc_alpha<i, j - 1, EoS, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                        enzyme_const, x, enzyme_const, invT);
    }
}

//
// Internal implementation detail.
// Dimensionless scaled derivative
// :math:`\lambda_{i,j} = (1/T)^{i}\,c^{\,j}\,
// \dfrac{\partial^{\,i+j}\alpha}{\partial (1/T)^{i}\,\partial c^{j}}`.
//
// Multiplying detail::calc_alpha by :math:`(1/T)^i c^j` cancels the units of the
// derivative, giving the dimensionless reduced derivatives that the property
// formulas below are written in terms of.
//
// :tparam i: Number of :math:`1/T`-derivatives.
// :tparam j: Number of ``c``-derivatives.
// :param eos: A single-contribution model (ideal or residual).
// :param c: Molar concentration [mol/m^3].
// :param x: Mole-fraction array [-].
// :param invT: Inverse temperature :math:`1/T` [1/K].
// :returns: :math:`\lambda_{i,j}` [-] (dimensionless).
//
// \ingroup core
template<int i, int j, EquationOfState EoS, std::floating_point Number>
[[nodiscard]] Number calc_lambda(const EoS& eos, const Number c, const Number* x, const Number invT)
{
    // TODO: better assertions?
    static_assert(i >= 0, "The template parameter `i` must be non-negative. It represent the number of derivatives "
                          "with respect to `invT`.");
    static_assert(
        j >= 0,
        "The template parameter `j` must be non-negative. It represent the number of derivatives with respect to `c`.");
    return fast_pow<Number, i>(invT) * fast_pow<Number, j>(c) * calc_alpha<i, j, EoS, Number>(eos, c, x, invT);
}

//
// Internal implementation detail.
// Total Helmholtz energy density :math:`\Psi = \sum_i \Psi_i` for one
//        contribution.
//
// Thin wrapper over the model's ``calc_helmholtz_density``, which accumulates the
// sum directly into a scalar (no per-component buffer). This is the scalar that
// detail::calc_dPsi_drhoi differentiates (reverse mode) to obtain chemical
// potentials.
//
// :param eos: A single-contribution model (ideal or residual).
// :param rho_i: Partial molar concentrations [mol/m^3].
// :param T: Temperature [K].
// :returns: Total Helmholtz energy density :math:`\Psi` [J/m^3].
//
// \ingroup core
template<EquationOfState EoS, std::floating_point Number>
[[nodiscard]] Number calc_Psi(const EoS& eos, const Number* FUGACITY_RESTRICT rho_i, const Number T)
{
    return eos.calc_helmholtz_density(rho_i, T);
}

//
// Internal implementation detail.
// ``i-th`` temperature derivative of the Helmholtz energy density,
//        :math:`\partial^{i}\Psi/\partial T^{i}`.
//
// Computed by recursively applying Enzyme forward-mode differentiation w.r.t.
// ``T``; the base case (``i == 0``) is detail::calc_Psi.
//
// :tparam i: Number of temperature derivatives.
// :param eos: A single-contribution model (ideal or residual).
// :param rho_i: Partial molar concentrations [mol/m^3].
// :param T: Temperature [K].
// :returns: :math:`\partial^{i}\Psi/\partial T^{i}` [J/(m^3 K^i)].
//
// \ingroup core
template<int i, EquationOfState EoS, std::floating_point Number>
[[nodiscard]] Number calc_dPsi_dT(const EoS& eos, const Number* FUGACITY_RESTRICT rho_i, const Number T)
{
    static_assert(i >= 0, "The template parameter `i` must be non-negative. It represent the number of derivatives "
                          "with respect to `T`.");
    if constexpr (i == 0) {
        return calc_Psi(eos, rho_i, T);
    }
    else {
        Number dT{1};
        return __enzyme_fwddiff<Number>((void*)calc_dPsi_dT<i - 1, EoS, Number>, enzyme_const, &eos, enzyme_const,
                                        rho_i, enzyme_dup, T, dT);
    }
}

//
// Internal implementation detail.
// Gradient of the Helmholtz energy density w.r.t. the partial molar
//        concentrations, :math:`\partial\Psi/\partial\rho_i` (a chemical potential).
//
// Uses Enzyme reverse mode on detail::calc_Psi. Currently only the first
// derivative (``i == 1``) is implemented; higher orders would require tensors.
//
// **Warning:** Enzyme reverse mode **accumulates** into ``dPsi_drho``. The caller must
//          zero it beforehand; calling repeatedly with the same buffer sums the
//          contributions (this is how the ideal and residual parts are combined).
//
// :tparam i: Derivative order; must be 1.
// :param eos: A single-contribution model (ideal or residual).
// :param rho_i: Partial molar concentrations [mol/m^3].
// :param T: Temperature [K].
// :param dPsi_drho: Output gradient (length ``eos.size()``), accumulated [J/mol].
//
// \ingroup core
template<int i, EquationOfState EoS, std::floating_point Number>
void calc_dPsi_drhoi(const EoS& eos, const Number* FUGACITY_RESTRICT rho_i, const Number T,
                     Number* FUGACITY_RESTRICT dPsi_drho)
{
    // Only implemented for 1 derivative because higher order would require tensors.
    static_assert(i > 0, "The template parameter `i` must be positive. It represents the number of derivatives with "
                         "respect to each `rho_i`. Use `calc_Psi()` for the '0th' derivative.");
    // NOTE: Enzyme reverse mode ACCUMULATES into `dPsi_drho`, so the caller must
    // zero it before the first call (and may sum several contributions by
    // calling repeatedly with the same buffer).
    if constexpr (i == 1) {
        // `calc_Psi` accumulates the per-component Helmholtz density into a local
        // scalar, so there is no caller-visible intermediate buffer to shadow:
        // `rho_i` is the only active input and the gradient lands in `dPsi_drho`.
        // FIXME: should <Number> be <void> since nothing is returned?
        __enzyme_autodiff<void>((void*)calc_Psi<EoS, Number>, enzyme_const, &eos, enzyme_dup, rho_i, dPsi_drho,
                                enzyme_const, T);
        return;
    }
    else {
        static_assert(false, "Higher order derivatives are not implemented (yet?) because they require tensors!");
        return;
    }
}
} // namespace detail

///
/// Total molar Helmholtz energy :math:`a = a^{\text{ideal}} + a^{\text{res}}`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar Helmholtz energy [J/mol].
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_helmholtz(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_helmholtz(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial a/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar Helmholtz energy [J/(mol K)].
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_helmholtz_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_helmholtz_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial a/\partial c)_{T,\boldsymbol{x}}` [J m^3/mol^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar Helmholtz energy [J m^3/mol^2].
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_helmholtz_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_helmholtz_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial a/\partial x_i` [J/mol].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/mol]. Overwritten.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_helmholtz_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_helmholtz_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_helmholtz(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Handle errors
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    return ideal.calc_helmholtz(c, x, T) + residual.calc_helmholtz(c, x, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_helmholtz_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Handle errors
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_helmholtz<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_helmholtz_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Handle errors
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_helmholtz<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c,
                                    dc, enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_helmholtz_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                       Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    // TODO: Handle errors
    __enzyme_autodiff<void>((void*)calc_helmholtz<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                            enzyme_dup, x, gradient, enzyme_const, T);
}

///
/// Pressure :math:`p = cRT\,(1 + \lambda^{\text{res}}_{0,1})`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Pressure [Pa].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_pressure(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_pressure(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial p/\partial T)_{c,\boldsymbol{x}}`
/// [Pa/K].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of pressure [Pa/K].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_pressure_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_pressure_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial p/\partial c)_{T,\boldsymbol{x}}` [Pa m^3/mol].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of pressure [Pa m^3/mol].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_pressure_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_pressure_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial p/\partial x_i` [Pa].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [Pa]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_pressure_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_pressure_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_pressure(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return c * R * T * (Number{1} + detail::calc_lambda<0, 1>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_pressure_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_pressure<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_pressure_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_pressure<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c,
                                    dc, enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_pressure_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                      Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    __enzyme_autodiff<void>((void*)calc_pressure<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                            enzyme_dup, x, gradient, enzyme_const, T);
}

///
/// Molar internal energy :math:`u = a + Ts`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar internal energy [J/mol].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_internal_energy(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_internal_energy(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial u/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar internal energy [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_internal_energy_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_internal_energy_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial u/\partial c)_{T,\boldsymbol{x}}` [J m^3/mol^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar internal energy [J m^3/mol^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_internal_energy_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_internal_energy_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial u/\partial x_i` [J/mol].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/mol]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_internal_energy_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_internal_energy_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_internal_energy(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * T * (detail::calc_lambda<1, 0>(ideal, c, x, invT) + detail::calc_lambda<1, 0>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_internal_energy_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_internal_energy<Ideal, Residual, Number>, enzyme_const, &eos,
                                    enzyme_const, c, enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_internal_energy_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_internal_energy<Ideal, Residual, Number>, enzyme_const, &eos,
                                    enzyme_dup, c, dc, enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_internal_energy_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                             Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_internal_energy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                            enzyme_dup, x, gradient, enzyme_const, T);
}

///
/// Molar enthalpy :math:`h = u + p/c`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar enthalpy [J/mol].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_enthalpy(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_enthalpy(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial h/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar enthalpy [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_enthalpy_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_enthalpy_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial h/\partial c)_{T,\boldsymbol{x}}` [J m^3/mol^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar enthalpy [J m^3/mol^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_enthalpy_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_enthalpy_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial h/\partial x_i` [J/mol].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/mol]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_enthalpy_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_enthalpy_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_enthalpy(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * T *
           (Number{1} + detail::calc_lambda<0, 1>(residual, c, x, invT) + detail::calc_lambda<1, 0>(ideal, c, x, invT) +
            detail::calc_lambda<1, 0>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_enthalpy_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_enthalpy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_enthalpy_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_enthalpy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c,
                                    dc, enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_enthalpy_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                      Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_enthalpy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                            enzyme_dup, x, gradient, enzyme_const, T);
}

///
/// Molar entropy :math:`s = - (\partial a / \partial T)_v`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar entropy [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_entropy(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_entropy(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial s/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K^2)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar entropy [J/(mol K^2)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_entropy_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_entropy_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial s/\partial c)_{T,\boldsymbol{x}}` [J m^3/(mol^2 K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar entropy [J m^3/(mol^2 K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_entropy_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_entropy_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial s/\partial x_i` [J/(mol K)].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/(mol K)]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_entropy_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_entropy_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_entropy(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * (detail::calc_lambda<1, 0>(ideal, c, x, invT) + detail::calc_lambda<1, 0>(residual, c, x, invT) -
                detail::calc_lambda<0, 0>(ideal, c, x, invT) - detail::calc_lambda<0, 0>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_entropy_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_entropy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_entropy_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_entropy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_entropy_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_entropy<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                            enzyme_dup, x, gradient, enzyme_const, T);
}

///
/// Molar Gibbs energy :math:`g = h - Ts`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar Gibbs energy [J/mol].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_gibbs(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_gibbs(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial g/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar Gibbs energy [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_gibbs_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_gibbs_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial g/\partial c)_{T,\boldsymbol{x}}` [J m^3/mol^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar Gibbs energy [J m^3/mol^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_gibbs_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_gibbs_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial g/\partial x_i` [J/mol].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/mol]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_gibbs_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_gibbs_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_gibbs(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * T *
           (Number{1} + detail::calc_lambda<0, 1>(residual, c, x, invT) + detail::calc_lambda<0, 0>(ideal, c, x, invT) +
            detail::calc_lambda<0, 0>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_gibbs_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_gibbs<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_gibbs_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_gibbs<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_gibbs_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_gibbs<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c, enzyme_dup,
                            x, gradient, enzyme_const, T);
}

///
/// Partial derivative of pressure w.r.t. concentration,
/// :math:`(\partial p/\partial c)_{T,x}`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: :math:`\partial p/\partial c` [Pa m^3/mol] (= J/mol).
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_dp_dc(eos, c, x.data(), T);
}

///
/// Temperature derivative of :math:`(\partial p/\partial c)_{T,\boldsymbol{x}}`
/// [Pa m^3/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Mixed pressure derivative [Pa m^3/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dc_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_dp_dc_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative of
/// :math:`(\partial p/\partial c)_{T,\boldsymbol{x}}` [Pa m^6/mol^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Second concentration derivative of pressure [Pa m^6/mol^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dc_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_dp_dc_dc(eos, c, x.data(), T);
}

///
/// Composition gradient of :math:`(\partial p/\partial c)_{T,\boldsymbol{x}}`
/// [Pa m^3/mol].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [Pa m^3/mol]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_dp_dc_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_dp_dc_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * T *
           (Number{1} + (Number{2} * detail::calc_lambda<0, 1>(residual, c, x, invT)) +
            detail::calc_lambda<0, 2>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dc_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_dp_dc<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dc_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_dp_dc<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_dp_dc_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_dp_dc<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c, enzyme_dup,
                            x, gradient, enzyme_const, T);
}

///
/// Partial derivative of pressure w.r.t. temperature,
/// :math:`(\partial p/\partial T)_{c,x}`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: :math:`\partial p/\partial T` [Pa/K].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_dp_dT(eos, c, x.data(), T);
}

///
/// Temperature derivative of :math:`(\partial p/\partial T)_{c,\boldsymbol{x}}`
/// [Pa/K^2].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Second temperature derivative of pressure [Pa/K^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dT_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_dp_dT_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative of
/// :math:`(\partial p/\partial T)_{c,\boldsymbol{x}}` [Pa m^3/(mol K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Mixed pressure derivative [Pa m^3/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_dp_dT_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_dp_dT_dc(eos, c, x.data(), T);
}

///
/// Composition gradient of :math:`(\partial p/\partial T)_{c,\boldsymbol{x}}`
/// [Pa/K].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [Pa/K]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_dp_dT_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_dp_dT_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return R * c *
           (Number{1} + detail::calc_lambda<0, 1>(residual, c, x, invT) -
            detail::calc_lambda<1, 1>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dT_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_dp_dT<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_dp_dT_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_dp_dT<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_dp_dT_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_dp_dT<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c, enzyme_dup,
                            x, gradient, enzyme_const, T);
}

///
/// Molar isochoric heat capacity :math:`c_v = (\partial u/\partial T)_{c}`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar :math:`c_v` [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cv(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_cv(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial c_v/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K^2)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar :math:`c_v` [J/(mol K^2)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cv_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_cv_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial c_v/\partial c)_{T,\boldsymbol{x}}` [J m^3/(mol^2 K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar :math:`c_v` [J m^3/(mol^2 K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cv_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_cv_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial c_v/\partial x_i` [J/(mol K)].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/(mol K)]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_cv_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_cv_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cv(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    return -R * (detail::calc_lambda<2, 0>(ideal, c, x, invT) + detail::calc_lambda<2, 0>(residual, c, x, invT));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cv_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_cv<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cv_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_cv<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_cv_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_cv<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c, enzyme_dup, x,
                            gradient, enzyme_const, T);
}

///
/// Molar isobaric heat capacity :math:`c_p`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :returns: Molar :math:`c_p` [J/(mol K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cp(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_cp(eos, c, x.data(), T);
}

///
/// Temperature derivative :math:`(\partial c_p/\partial T)_{c,\boldsymbol{x}}`
/// [J/(mol K^2)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :returns: Temperature derivative of molar :math:`c_p` [J/(mol K^2)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cp_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_cp_dT(eos, c, x.data(), T);
}

///
/// Concentration derivative
/// :math:`(\partial c_p/\partial c)_{T,\boldsymbol{x}}` [J m^3/(mol^2 K)].
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :returns: Concentration derivative of molar :math:`c_p` [J m^3/(mol^2 K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_cp_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_cp_dc(eos, c, x.data(), T);
}

///
/// Composition gradient :math:`\partial c_p/\partial x_i` [J/(mol K)].
///
/// Treat the mole fractions as independent coordinates.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param gradient: Output composition gradient [J/(mol K)]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_cp_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_cp_dx(eos, c, x.data(), T, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cp(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    const Ideal& ideal = eos.ideal();
    const Residual& residual = eos.residual();
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invT = Number{1} / T;
    const Number lr_01 = detail::calc_lambda<0, 1>(residual, c, x, invT);
    return R * (-detail::calc_lambda<2, 0>(ideal, c, x, invT) - detail::calc_lambda<2, 0>(residual, c, x, invT) +
                (detail::fast_pow<Number, 2>(Number{1} + lr_01 - detail::calc_lambda<1, 1>(residual, c, x, invT)) /
                 (Number{1} + (Number{2} * lr_01) + detail::calc_lambda<0, 2>(residual, c, x, invT))));
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cp_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_cp<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c,
                                    enzyme_const, x, enzyme_dup, T, dT);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_cp_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_cp<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_dup, c, dc,
                                    enzyme_const, x, enzyme_const, T);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_cp_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_cp<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const, c, enzyme_dup, x,
                            gradient, enzyme_const, T);
}

// FIXME: the effective_molar_mass parameter should be turned into a function that can compute the molar mass
//        This mainly affects taking derivatives correctly, so it is low priority
///
/// Squared speed of sound :math:`w^2 = c_p\,(\partial p/\partial c)/(M\,c_v)`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :param effective_molar_mass: Mixture molar mass :math:`M` [kg/mol].
/// :returns: Squared speed of sound [m^2/s^2].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_sound_speed_squared(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T,
                                const Number effective_molar_mass)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    return calc_sound_speed_squared(eos, c, x.data(), T, effective_molar_mass);
}

///
/// Temperature derivative
/// :math:`(\partial w^2/\partial T)_{c,\boldsymbol{x},M}` [m^2/(s^2 K)].
///
/// Hold the supplied effective molar mass constant.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K].
/// :param effective_molar_mass: Mixture molar mass :math:`M` [kg/mol]. Held constant in the derivative.
/// :returns: Temperature derivative of squared speed of sound [m^2/(s^2 K)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_sound_speed_squared_dT(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T,
                                   const Number effective_molar_mass)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_sound_speed_squared_dT(eos, c, x.data(), T, effective_molar_mass);
}

///
/// Concentration derivative
/// :math:`(\partial w^2/\partial c)_{T,\boldsymbol{x},M}` [m^5/(mol s^2)].
///
/// Hold the supplied effective molar mass constant.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-]. Held constant in the derivative.
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param effective_molar_mass: Mixture molar mass :math:`M` [kg/mol]. Held constant in the derivative.
/// :returns: Concentration derivative of squared speed of sound [m^5/(mol s^2)].
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V>
Number calc_sound_speed_squared_dc(const EoS<Ideal, Residual>& eos, const Number c, V& x, const Number T,
                                   const Number effective_molar_mass)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    return calc_sound_speed_squared_dc(eos, c, x.data(), T, effective_molar_mass);
}

///
/// Composition gradient :math:`\partial w^2/\partial x_i` [m^2/s^2].
///
/// Treat the mole fractions as independent coordinates and hold the supplied
/// effective molar mass constant.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3]. Held constant in the derivative.
/// :param x: Mole fractions [-].
/// :param T: Temperature [K]. Held constant in the derivative.
/// :param effective_molar_mass: Mixture molar mass :math:`M` [kg/mol]. Held constant in the derivative.
/// :param gradient: Output composition gradient [m^2/s^2]. Overwritten.
/// :precondition: ``T > 0``.
/// :precondition: ``x.size() == eos.size()``
/// :precondition: ``gradient.size() == eos.size()``
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, typename V1, typename V2>
void calc_sound_speed_squared_dx(const EoS<Ideal, Residual>& eos, const Number c, V1& x, const Number T,
                                 const Number effective_molar_mass, V2& gradient)
{
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_ASSERT(x.size() == eos.size());
    calc_sound_speed_squared_dx(eos, c, x.data(), T, effective_molar_mass, gradient.data());
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_sound_speed_squared(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                                const Number effective_molar_mass)
{
    // TODO: If I need more derivatives of sound speed, then I might need to compute the effective molar mass inside
    // this function
    // TODO: Consider a custom assertion with a better error message
    FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T);
    Number cp = calc_cp(eos, c, x, T);
    Number cv = calc_cv(eos, c, x, T);
    Number dp_dc = calc_dp_dc(eos, c, x, T);
    return cp * dp_dc / (effective_molar_mass * cv);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_sound_speed_squared_dT(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                                   const Number effective_molar_mass)
{
    Number dT{1};
    return __enzyme_fwddiff<Number>((void*)calc_sound_speed_squared<Ideal, Residual, Number>, enzyme_const, &eos,
                                    enzyme_const, c, enzyme_const, x, enzyme_dup, T, dT, enzyme_const,
                                    effective_molar_mass);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
Number calc_sound_speed_squared_dc(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                                   const Number effective_molar_mass)
{
    Number dc{1};
    return __enzyme_fwddiff<Number>((void*)calc_sound_speed_squared<Ideal, Residual, Number>, enzyme_const, &eos,
                                    enzyme_dup, c, dc, enzyme_const, x, enzyme_const, T, enzyme_const,
                                    effective_molar_mass);
}

template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number>
void calc_sound_speed_squared_dx(const EoS<Ideal, Residual>& eos, const Number c, const Number* x, const Number T,
                                 const Number effective_molar_mass, Number* gradient)
{
    std::fill_n(gradient, eos.size(), Number{0});
    __enzyme_autodiff<void>((void*)calc_sound_speed_squared<Ideal, Residual, Number>, enzyme_const, &eos, enzyme_const,
                            c, enzyme_dup, x, gradient, enzyme_const, T, enzyme_const, effective_molar_mass);
}

// TODO: add derivatives of vector-valued functions
///
/// Chemical potentials of all components, including ideal and residual parts.
///
/// .. math::
///
///    \mu_i=\left(\frac{\partial\Psi}{\partial\rho_i}\right)_{T,\rho_{j\ne i}}.
///
/// :param eos: The equation of state.
/// :param rho_i: Partial molar concentrations [mol/m^3].
/// :param T: Temperature [K].
/// :param chemical_potential: Output chemical potentials [J/mol]. Overwritten (zeroed then filled).
/// :precondition: ``rho_i.size() == eos.size()``
/// :precondition: ``chemical_potential.size() == eos.size()``
/// :precondition: ``T > 0``.
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, std::size_t N>
void calc_chemical_potential(const EoS<Ideal, Residual>& eos, std::span<const Number, N> rho_i, const Number T,
                             std::span<Number, N> chemical_potential)
{
    FUGACITY_ASSERT(rho_i.size() == eos.size());
    FUGACITY_ASSERT(chemical_potential.size() == eos.size());
    // Enzyme reverse mode accumulates into the output, so it must start zeroed;
    // the ideal and residual contributions are then summed in place.
    std::fill(chemical_potential.begin(), chemical_potential.end(), Number{0});
    detail::calc_dPsi_drhoi<1>(eos.ideal(), rho_i.data(), T, chemical_potential.data());
    detail::calc_dPsi_drhoi<1>(eos.residual(), rho_i.data(), T, chemical_potential.data());
}

///
/// Natural logarithms of the fugacity coefficients,
/// :math:`\ln\varphi_i = \mu_i^{\text{res}}/(RT) - \ln Z`.
///
/// :param eos: The equation of state.
/// :param c: Molar concentration [mol/m^3].
/// :param x: Mole fractions [-].
/// :param T: Temperature [K].
/// :param rho_i: Partial molar concentrations [mol/m^3] (should equal ``x*c``).
/// :param log_fug_coeff: Output :math:`\ln\varphi_i` [-] (length ``eos.size()``).
///             Overwritten (zeroed then filled).
/// :precondition: ``x.size() == rho_i.size() == log_fug_coeff.size() == eos.size()``.
/// :precondition: ``rho_i[i] == c*x[i]``.
/// :precondition: ``T > 0`` and :math:`Z>0`.
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, std::size_t N>
void calc_log_fugacity_coeff(const EoS<Ideal, Residual>& eos, const Number c, std::span<const Number, N> x,
                             const Number T, const std::span<const Number, N> rho_i, std::span<Number, N> log_fug_coeff)
{
    const Number invT = Number{1} / T;
    const Number Z = Number{1} + detail::calc_lambda<0, 1>(eos.residual(), c, x.data(), invT);
    const Number lnZ = std::log(Z);

    // Enzyme reverse mode accumulates into the output, so zero it first.
    std::fill(log_fug_coeff.begin(), log_fug_coeff.end(), Number{0});
    detail::calc_dPsi_drhoi<1>(eos.residual(), rho_i.data(), T, log_fug_coeff.data());
    constexpr Number R = ideal_gas_constant<Number>;
    const Number invRT = Number{1} / (R * T);
    for (auto& ln_phi_i : log_fug_coeff) {
        ln_phi_i *= invRT;
        ln_phi_i -= lnZ;
    }
}

///
/// Fugacities :math:`f_i = \rho_i RT\,\exp(\mu_i^{\text{res}}/(RT))`.
///
/// :param eos: The equation of state.
/// :param rho_i: Partial molar concentrations [mol/m^3] (length ``eos.size()``).
/// :param T: Temperature [K].
/// :param fugacity: Output fugacities [Pa] (length ``eos.size()``).
///             Overwritten (zeroed then filled).
/// :precondition: ``rho_i.size() == fugacity.size() == eos.size()``.
/// :precondition: ``T > 0``.
///
/// \ingroup core
template<IdealEoS Ideal, ResidualEoS Residual, std::floating_point Number, std::size_t N>
void calc_fugacity(const EoS<Ideal, Residual>& eos, std::span<const Number, N> rho_i, const Number T,
                   std::span<Number, N> fugacity)
{
    // Enzyme reverse mode accumulates into the output, so zero it first. After
    // the call `fugacity[idx]` holds the residual chemical potential mu_i^res.
    std::fill(fugacity.begin(), fugacity.end(), Number{0});
    detail::calc_dPsi_drhoi<1>(eos.residual(), rho_i.data(), T, fugacity.data());
    constexpr Number R = ideal_gas_constant<Number>;
    const Number RT = R * T;
    const Number invRT = Number{1} / RT;
    for (std::size_t idx = 0; idx < fugacity.size(); ++idx) {
        fugacity[idx] = rho_i[idx] * RT * std::exp(fugacity[idx] * invRT);
    }
}

} // namespace fugacity
