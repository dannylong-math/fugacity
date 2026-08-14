#pragma once
///
/// NASA-9 ideal model.
///

#include "fugacity/core/concepts.hpp"
#include "fugacity/core/eos_base.hpp"
#include "fugacity/core/horner.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/core/xlnx.hpp"

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

namespace fugacity {

///
/// NASA-9 ideal model.
///
/// For each species, the NASA-9 equations are
///
/// .. math::
///
///    \frac{c_{p,i}^\circ}{R}
///      =a_{0,i}T^{-2}+a_{1,i}T^{-1}+a_{2,i}+a_{3,i}T
///       +a_{4,i}T^2+a_{5,i}T^3+a_{6,i}T^4,
///
/// .. math::
///
///    \frac{h_i^\circ}{RT}
///      =-a_{0,i}T^{-2}+a_{1,i}\frac{\ln T}{T}+a_{2,i}
///       +\frac{a_{3,i}}{2}T+\frac{a_{4,i}}{3}T^2
///       +\frac{a_{5,i}}{4}T^3+\frac{a_{6,i}}{5}T^4+\frac{a_{7,i}}{T},
///
/// .. math::
///
///    \frac{s_i^\circ}{R}
///      =-\frac{a_{0,i}}{2}T^{-2}-a_{1,i}T^{-1}+a_{2,i}\ln T
///       +a_{3,i}T+\frac{a_{4,i}}{2}T^2+\frac{a_{5,i}}{3}T^3
///       +\frac{a_{6,i}}{4}T^4+a_{8,i}.
///
/// The molar ideal Helmholtz energy is
///
/// .. math::
///
///    a^\mathrm{ideal}
///    =\sum_i x_i\left[
///      h_i^\circ-Ts_i^\circ
///      +RT\ln\!\left(\frac{x_i cRT}{p_{i,\mathrm{ref}}}\right)
///      \right]-RT.
///
/// Supply a single coefficient range that is valid at the evaluation
/// temperature. The model does not select among temperature ranges.
///
/// .. code-block:: cpp
///
///    using Ideal = fugacity::Nasa9<1>;
///    const std::array<Ideal::SpeciesInput, 1> species{{{
///        .a0 = 2.210371497e4, .a1 = -3.818461820e2,
///        .a2 = 6.082738360, .a3 = -8.530914410e-3,
///        .a4 = 1.384646189e-5, .a5 = -9.625793620e-9,
///        .a6 = 2.519705809e-12, .a7 = 7.108460860e2,
///        .a8 = -1.076003744e1, .T_ref = 298.15, .p_ref = 1.0e5,
///    }}};
///    const Ideal ideal{species};
///    const fugacity::EoS eos{ideal, fugacity::NoResidual<1>{}};
///
///    const std::array<double, 1> x{1.0};
///    const double cp = fugacity::calc_cp(eos, 40.0, x, 350.0);
///
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \ingroup ideal-models
template<std::size_t N = std::dynamic_extent> class Nasa9 : public BaseEoS<N>, public BaseIdealEoS {
public:
    /// NASA-9 coefficients and standard-state data for one species.
    struct SpeciesInput {
        double a0;    ///< Coefficient :math:`a_0` [K^2].
        double a1;    ///< Coefficient :math:`a_1` [K].
        double a2;    ///< Coefficient :math:`a_2` [-].
        double a3;    ///< Coefficient :math:`a_3` [1/K].
        double a4;    ///< Coefficient :math:`a_4` [1/K^2].
        double a5;    ///< Coefficient :math:`a_5` [1/K^3].
        double a6;    ///< Coefficient :math:`a_6` [1/K^4].
        double a7;    ///< Enthalpy integration coefficient :math:`a_7` [K].
        double a8;    ///< Entropy integration coefficient :math:`a_8` [-].
        double T_ref; ///< Temperature paired with the standard pressure [K].
        double p_ref; ///< Standard-state pressure :math:`p_\mathrm{ref}` [Pa].
    };

    ///
    /// Construct a fixed-size model.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// \id fixed-size
    ///
    explicit Nasa9(const std::array<SpeciesInput, N>& inputs)
        requires(N != std::dynamic_extent)
    {
        for (std::size_t i = 0; i < N; ++i) {
            scatter(i, N, inputs[i]);
        }
    }

    ///
    /// Construct a runtime-size model.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species. ``size()`` is
    ///                set to ``inputs.size()``.
    /// \id runtime-size
    ///
    explicit Nasa9(std::span<const SpeciesInput> inputs)
        requires(N == std::dynamic_extent)
        : BaseEoS<N>(inputs.size())
    {
        const std::size_t n = inputs.size();
        data_.resize(num_params * n);
        for (std::size_t i = 0; i < n; ++i) {
            scatter(i, n, inputs[i]);
        }
    }

    ///
    /// Return the molar ideal Helmholtz energy.
    ///
    /// :param c: Molar concentration [mol/m^3].
    /// :param x: Mole-fraction array [-].
    /// :param T: Temperature [K].
    /// :returns: Molar Helmholtz energy [J/mol].
    ///
    template<std::floating_point Number> [[nodiscard]] Number calc_helmholtz(Number c, const Number* x, Number T) const
    {
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        const Number lnC = std::log(c);
        const Number lnT = std::log(T);
        const Number one_over_T = Number{1} / T;

        Number a{0};
        for (std::size_t i = 0; i < n; ++i) {
            // The NASA-9 Helmholtz contribution is a 5th-degree polynomial in T
            // plus 1/T and lnT-dependent constants, folded into coeffs[0].
            std::array<Number, 6> coeffs;
            coeffs[0] = (data_[(col_a1 * n) + i] * (lnT + Number{1})) + data_[(col_a7 * n) + i] -
                        (data_[(col_a0_over_2 * n) + i] * one_over_T);
            coeffs[1] = (data_[(col_a2 * n) + i] * (Number{1} - lnT)) - data_[(col_a8 * n) + i] + lnT + lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a3_over_2 * n) + i];
            coeffs[3] = -data_[(col_a4_over_6 * n) + i];
            coeffs[4] = -data_[(col_a5_over_12 * n) + i];
            coeffs[5] = -data_[(col_a6_over_20 * n) + i];
            a += R * ((x[i] * eval_polynomial<5>(coeffs, T)) + (T * Number{2} * xlnx<0>(x[i])));
        }
        return a;
    }

    ///
    /// Return the ideal Helmholtz energy density
    /// :math:`\Psi^\mathrm{ideal}=c a^\mathrm{ideal}`.
    ///
    /// :param rho_i: Partial molar concentrations [mol/m^3].
    /// :param T: Temperature [K].
    /// :returns: Helmholtz energy density [J/m^3].
    ///
    template<std::floating_point Number>
    [[nodiscard]] Number calc_helmholtz_density(const Number* rho_i, Number T) const
    {
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        Number c{0};
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
        }
        const Number lnC = std::log(c);
        const Number lnT = std::log(T);
        const Number one_over_T = Number{1} / T;

        Number psi{0};
        for (std::size_t i = 0; i < n; ++i) {
            std::array<Number, 6> coeffs;
            coeffs[0] = (data_[(col_a1 * n) + i] * (lnT + Number{1})) + data_[(col_a7 * n) + i] -
                        (data_[(col_a0_over_2 * n) + i] * one_over_T);
            coeffs[1] = (data_[(col_a2 * n) + i] * (Number{1} - lnT)) - data_[(col_a8 * n) + i] + lnT - lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a3_over_2 * n) + i];
            coeffs[3] = -data_[(col_a4_over_6 * n) + i];
            coeffs[4] = -data_[(col_a5_over_12 * n) + i];
            coeffs[5] = -data_[(col_a6_over_20 * n) + i];
            psi += R * ((rho_i[i] * eval_polynomial<5>(coeffs, T)) + (T * Number{2} * xlnx<0>(rho_i[i])));
        }
        return psi;
    }

    ///
    /// Return a per-component decomposition of the ideal Helmholtz energy density.
    ///
    /// :param rho_i: Partial molar concentrations [mol/m^3].
    /// :param T: Temperature [K].
    /// :param out: Per-component Helmholtz energy density [J/m^3]; length ``size()``.
    ///
    template<std::floating_point Number> void calc_partial_helmholtz(const Number* rho_i, Number T, Number* out) const
    {
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        Number c{0};
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
        }
        const Number lnC = std::log(c);
        const Number lnT = std::log(T);
        const Number one_over_T = Number{1} / T;

        for (std::size_t i = 0; i < n; ++i) {
            std::array<Number, 6> coeffs;
            coeffs[0] = (data_[(col_a1 * n) + i] * (lnT + Number{1})) + data_[(col_a7 * n) + i] -
                        (data_[(col_a0_over_2 * n) + i] * one_over_T);
            coeffs[1] = (data_[(col_a2 * n) + i] * (Number{1} - lnT)) - data_[(col_a8 * n) + i] + lnT - lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a3_over_2 * n) + i];
            coeffs[3] = -data_[(col_a4_over_6 * n) + i];
            coeffs[4] = -data_[(col_a5_over_12 * n) + i];
            coeffs[5] = -data_[(col_a6_over_20 * n) + i];
            out[i] = R * ((rho_i[i] * eval_polynomial<5>(coeffs, T)) + (T * Number{2} * xlnx<0>(rho_i[i])));
        }
    }

private:
    // Same deliberate choices as Nasa7 (see the notes there): flat kernels for
    // Enzyme, Structure-of-Arrays storage, coefficients pre-divided by the
    // constants the Helmholtz polynomial needs.

    // Column layout of the SoA storage.
    static constexpr std::size_t col_a0_over_2 = 0;
    static constexpr std::size_t col_a1 = 1;
    static constexpr std::size_t col_a2 = 2;
    static constexpr std::size_t col_a3_over_2 = 3;
    static constexpr std::size_t col_a4_over_6 = 4;
    static constexpr std::size_t col_a5_over_12 = 5;
    static constexpr std::size_t col_a6_over_20 = 6;
    static constexpr std::size_t col_a7 = 7;
    static constexpr std::size_t col_a8 = 8;
    static constexpr std::size_t col_ln_cref_Tref = 9; // ln(c_ref * T_ref)
    static constexpr std::size_t num_params = 10;

    using Storage =
        std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, num_params * N>>;

    // Derive species i's parameters from `in` and write them into the columns.
    void scatter(std::size_t i, std::size_t n, const SpeciesInput& in)
    {
        constexpr double R = ideal_gas_constant<double>;
        const double c_ref = in.p_ref / (R * in.T_ref);
        data_[(col_a0_over_2 * n) + i] = in.a0 / 2.;
        data_[(col_a1 * n) + i] = in.a1;
        data_[(col_a2 * n) + i] = in.a2;
        data_[(col_a3_over_2 * n) + i] = in.a3 / 2.;
        data_[(col_a4_over_6 * n) + i] = in.a4 / 6.;
        data_[(col_a5_over_12 * n) + i] = in.a5 / 12.;
        data_[(col_a6_over_20 * n) + i] = in.a6 / 20.;
        data_[(col_a7 * n) + i] = in.a7;
        data_[(col_a8 * n) + i] = in.a8;
        data_[(col_ln_cref_Tref * n) + i] = std::log(c_ref * in.T_ref);
    }

    Storage data_{}; // Column-major (SoA) parameter storage.
};

static_assert(IdealEoS<Nasa9<2>>, "Nasa9 must satisfy the IdealEoS concept.");
static_assert(IdealEoS<Nasa9<std::dynamic_extent>>, "Nasa9 must satisfy the IdealEoS concept.");
} // namespace fugacity
