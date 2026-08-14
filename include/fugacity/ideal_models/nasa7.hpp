#pragma once
///
/// NASA-7 ideal model.
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
/// NASA-7 ideal model.
///
/// For each species, the NASA-7 equations are
///
/// .. math::
///
///    \frac{c_{p,i}^\circ}{R}
///      =a_{0,i}+a_{1,i}T+a_{2,i}T^2+a_{3,i}T^3+a_{4,i}T^4,
///
/// .. math::
///
///    \frac{h_i^\circ}{RT}
///      =a_{0,i}+\frac{a_{1,i}}{2}T+\frac{a_{2,i}}{3}T^2
///       +\frac{a_{3,i}}{4}T^3+\frac{a_{4,i}}{5}T^4+\frac{a_{5,i}}{T},
///
/// .. math::
///
///    \frac{s_i^\circ}{R}
///      =a_{0,i}\ln T+a_{1,i}T+\frac{a_{2,i}}{2}T^2
///       +\frac{a_{3,i}}{3}T^3+\frac{a_{4,i}}{4}T^4+a_{6,i}.
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
/// temperature. The model does not select between low- and high-temperature
/// coefficient ranges.
///
/// .. code-block:: cpp
///
///    using Ideal = fugacity::Nasa7<1>;
///    const std::array<Ideal::SpeciesInput, 1> species{{{
///        .a0 = 3.53100528, .a1 = -1.23660988e-4,
///        .a2 = -5.02999433e-7, .a3 = 2.43530612e-9,
///        .a4 = -1.40881235e-12, .a5 = -1046.97628,
///        .a6 = 2.96747038, .T_ref = 298.15, .p_ref = 1.0e5,
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
template<std::size_t N = std::dynamic_extent> class Nasa7 : public BaseEoS<N>, public BaseIdealEoS {
public:
    /// NASA-7 coefficients and standard-state data for one species.
    struct SpeciesInput {
        double a0;    ///< Coefficient :math:`a_0` [-].
        double a1;    ///< Coefficient :math:`a_1` [1/K].
        double a2;    ///< Coefficient :math:`a_2` [1/K^2].
        double a3;    ///< Coefficient :math:`a_3` [1/K^3].
        double a4;    ///< Coefficient :math:`a_4` [1/K^4].
        double a5;    ///< Enthalpy integration coefficient :math:`a_5` [K].
        double a6;    ///< Entropy integration coefficient :math:`a_6` [-].
        double T_ref; ///< Temperature paired with the standard pressure [K].
        double p_ref; ///< Standard-state pressure :math:`p_\mathrm{ref}` [Pa].
    };

    ///
    /// Construct a fixed-size model.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// \id fixed-size
    ///
    explicit Nasa7(const std::array<SpeciesInput, N>& inputs)
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
    explicit Nasa7(std::span<const SpeciesInput> inputs)
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

        Number a{0};
        for (std::size_t i = 0; i < n; ++i) {
            // The NASA-7 Helmholtz contribution is a 5th-degree polynomial in T.
            std::array<Number, 6> coeffs;
            coeffs[0] = data_[(col_a5 * n) + i];
            coeffs[1] = (data_[(col_a0 * n) + i] * (Number{1} - lnT)) - data_[(col_a6 * n) + i] + lnT + lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a1_over_2 * n) + i];
            coeffs[3] = -data_[(col_a2_over_6 * n) + i];
            coeffs[4] = -data_[(col_a3_over_12 * n) + i];
            coeffs[5] = -data_[(col_a4_over_20 * n) + i];
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

        Number psi{0};
        for (std::size_t i = 0; i < n; ++i) {
            std::array<Number, 6> coeffs;
            coeffs[0] = data_[(col_a5 * n) + i];
            coeffs[1] = (data_[(col_a0 * n) + i] * (Number{1} - lnT)) - data_[(col_a6 * n) + i] + lnT - lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a1_over_2 * n) + i];
            coeffs[3] = -data_[(col_a2_over_6 * n) + i];
            coeffs[4] = -data_[(col_a3_over_12 * n) + i];
            coeffs[5] = -data_[(col_a4_over_20 * n) + i];
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

        for (std::size_t i = 0; i < n; ++i) {
            std::array<Number, 6> coeffs;
            coeffs[0] = data_[(col_a5 * n) + i];
            coeffs[1] = (data_[(col_a0 * n) + i] * (Number{1} - lnT)) - data_[(col_a6 * n) + i] + lnT - lnC -
                        data_[(col_ln_cref_Tref * n) + i] - Number{1};
            coeffs[2] = -data_[(col_a1_over_2 * n) + i];
            coeffs[3] = -data_[(col_a2_over_6 * n) + i];
            coeffs[4] = -data_[(col_a3_over_12 * n) + i];
            coeffs[5] = -data_[(col_a4_over_20 * n) + i];
            out[i] = R * ((rho_i[i] * eval_polynomial<5>(coeffs, T)) + (T * Number{2} * xlnx<0>(rho_i[i])));
        }
    }

private:
    // Implementation notes (deliberate, benchmarked choices):
    //  * The kernels are written out in full -- no CRTP indirection, no
    //    parameter-struct reflection -- so the expression Enzyme differentiates
    //    stays flat and the autodiff is robust at every optimization level.
    //  * The derived per-species parameters are stored Structure-of-Arrays
    //    (column-major: one contiguous column per parameter across all
    //    species), which benchmarked faster than an Array-of-Structs layout.
    //  * The a1..a4 coefficients are stored pre-divided by the constants the
    //    Helmholtz polynomial needs, so the kernels only negate them.

    // Column layout of the SoA storage.
    static constexpr std::size_t col_a0 = 0;
    static constexpr std::size_t col_a1_over_2 = 1;
    static constexpr std::size_t col_a2_over_6 = 2;
    static constexpr std::size_t col_a3_over_12 = 3;
    static constexpr std::size_t col_a4_over_20 = 4;
    static constexpr std::size_t col_a5 = 5;
    static constexpr std::size_t col_a6 = 6;
    static constexpr std::size_t col_ln_cref_Tref = 7; // ln(c_ref * T_ref)
    static constexpr std::size_t num_params = 8;

    using Storage =
        std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, num_params * N>>;

    // Derive species i's parameters from `in` and write them into the columns.
    void scatter(std::size_t i, std::size_t n, const SpeciesInput& in)
    {
        constexpr double R = ideal_gas_constant<double>;
        const double c_ref = in.p_ref / (R * in.T_ref);
        data_[(col_a0 * n) + i] = in.a0;
        data_[(col_a1_over_2 * n) + i] = in.a1 / 2.;
        data_[(col_a2_over_6 * n) + i] = in.a2 / 6.;
        data_[(col_a3_over_12 * n) + i] = in.a3 / 12.;
        data_[(col_a4_over_20 * n) + i] = in.a4 / 20.;
        data_[(col_a5 * n) + i] = in.a5;
        data_[(col_a6 * n) + i] = in.a6;
        data_[(col_ln_cref_Tref * n) + i] = std::log(c_ref * in.T_ref);
    }

    Storage data_{}; // Column-major (SoA) parameter storage.
};

static_assert(IdealEoS<Nasa7<2>>, "Nasa7 must satisfy the IdealEoS concept.");
static_assert(IdealEoS<Nasa7<std::dynamic_extent>>, "Nasa7 must satisfy the IdealEoS concept.");
} // namespace fugacity
