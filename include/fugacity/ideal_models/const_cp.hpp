#pragma once
///
/// Constant-heat-capacity ideal model.
///

#include "fugacity/core/concepts.hpp"
#include "fugacity/core/eos_base.hpp"
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
/// Constant-heat-capacity ideal model.
///
/// For each species, the standard-state enthalpy and entropy are
///
/// .. math::
///
///    h_i^\circ(T)=h_{i,\mathrm{ref}}+c_{p,i}(T-T_{i,\mathrm{ref}}),
///
/// .. math::
///
///    s_i^\circ(T)=s_{i,\mathrm{ref}}
///       +c_{p,i}\ln\!\left(\frac{T}{T_{i,\mathrm{ref}}}\right).
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
/// Here :math:`c` is molar concentration [mol/m^3], :math:`T` is temperature
/// [K], and :math:`x_i` is mole fraction [-].
///
/// .. code-block:: cpp
///
///    using Ideal = fugacity::ConstantCp<2>;
///    const std::array<Ideal::SpeciesInput, 2> species{{
///        {.T_ref = 298.15, .p_ref = 1.0e5, .c_p = 29.1,
///         .h_ref = 0.0, .s_ref = 191.6},
///        {.T_ref = 298.15, .p_ref = 1.0e5, .c_p = 33.6,
///         .h_ref = 0.0, .s_ref = 205.2},
///    }};
///    const Ideal ideal{species};
///    const fugacity::EoS eos{ideal, fugacity::NoResidual<2>{}};
///
///    const std::array<double, 2> x{0.5, 0.5};
///    const double h = fugacity::calc_enthalpy(eos, 40.0, x, 350.0);
///
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \ingroup ideal-models
template<std::size_t N = std::dynamic_extent> class ConstantCp : public BaseEoS<N>, public BaseIdealEoS {
public:
    ///
    /// Thermodynamic data for one species.
    struct SpeciesInput {
        double T_ref; ///< Reference temperature :math:`T_\mathrm{ref}` [K].
        double p_ref; ///< Reference pressure :math:`p_\mathrm{ref}` [Pa].
        double c_p;   ///< Isobaric molar heat capacity :math:`c_p` [J/(mol K)].
        double h_ref; ///< Reference molar enthalpy :math:`h_\mathrm{ref}` [J/mol].
        double s_ref; ///< Reference molar entropy :math:`s_\mathrm{ref}` [J/(mol K)].
    };

    ///
    /// Construct a fixed-size model.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// \id fixed-size
    ///
    explicit ConstantCp(const std::array<SpeciesInput, N>& inputs)
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
    explicit ConstantCp(std::span<const SpeciesInput> inputs)
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
        // Pre-calculations hoisted out of the component loop.
        const Number TlnT = T * std::log(T);
        const Number RTlnC = R * T * std::log(c);

        Number a{0};
        for (std::size_t i = 0; i < n; ++i) {
            const double T_ref = data_[(col_T_ref * n) + i];
            const double R_ln_c_ref = data_[(col_R_ln_c_ref * n) + i];
            const double c_p = data_[(col_c_p * n) + i];
            const double ln_T_ref = data_[(col_ln_T_ref * n) + i];
            const double h_ref = data_[(col_h_ref * n) + i];
            const double s_ref = data_[(col_s_ref * n) + i];
            a += (Number{2} * R * T * xlnx<0>(x[i])) +
                 (x[i] * (h_ref + (c_p * (T - T_ref)) - (T * s_ref) + (c_p * T * ln_T_ref) + RTlnC - (T * R_ln_c_ref) +
                          ((R - c_p) * TlnT) - (R * T * ln_T_ref) - (R * T)));
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
        // c = sum_i rho_i, then the hoisted pre-calculations.
        Number c{0};
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
        }
        const Number TlnT = T * std::log(T);
        const Number RTlnC = R * T * std::log(c);

        Number psi{0};
        for (std::size_t i = 0; i < n; ++i) {
            const double T_ref = data_[(col_T_ref * n) + i];
            const double R_ln_c_ref = data_[(col_R_ln_c_ref * n) + i];
            const double c_p = data_[(col_c_p * n) + i];
            const double ln_T_ref = data_[(col_ln_T_ref * n) + i];
            const double h_ref = data_[(col_h_ref * n) + i];
            const double s_ref = data_[(col_s_ref * n) + i];
            psi += (Number{2} * R * T * xlnx(rho_i[i])) +
                   (rho_i[i] * (h_ref + (c_p * (T - T_ref)) - (T * s_ref) + (c_p * T * ln_T_ref) - RTlnC -
                                (T * R_ln_c_ref) + ((R - c_p) * TlnT) - (R * T * ln_T_ref) - (R * T)));
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
        const Number TlnT = T * std::log(T);
        const Number RTlnC = R * T * std::log(c);

        for (std::size_t i = 0; i < n; ++i) {
            const double T_ref = data_[(col_T_ref * n) + i];
            const double R_ln_c_ref = data_[(col_R_ln_c_ref * n) + i];
            const double c_p = data_[(col_c_p * n) + i];
            const double ln_T_ref = data_[(col_ln_T_ref * n) + i];
            const double h_ref = data_[(col_h_ref * n) + i];
            const double s_ref = data_[(col_s_ref * n) + i];
            out[i] = (Number{2} * R * T * xlnx(rho_i[i])) +
                     (rho_i[i] * (h_ref + (c_p * (T - T_ref)) - (T * s_ref) + (c_p * T * ln_T_ref) - RTlnC -
                                  (T * R_ln_c_ref) + ((R - c_p) * TlnT) - (R * T * ln_T_ref) - (R * T)));
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

    // Column layout of the SoA storage.
    static constexpr std::size_t col_T_ref = 0;      // T_ref [K]
    static constexpr std::size_t col_R_ln_c_ref = 1; // R ln(c_ref) [J/(mol K)]
    static constexpr std::size_t col_c_p = 2;        // c_p [J/(mol K)]
    static constexpr std::size_t col_ln_T_ref = 3;   // ln(T_ref) [-]
    static constexpr std::size_t col_h_ref = 4;      // h_ref [J/mol]
    static constexpr std::size_t col_s_ref = 5;      // s_ref [J/(mol K)]
    static constexpr std::size_t num_params = 6;

    using Storage =
        std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, num_params * N>>;

    // Derive species i's parameters from `in` and write them into the columns.
    void scatter(std::size_t i, std::size_t n, const SpeciesInput& in)
    {
        constexpr double R = ideal_gas_constant<double>;
        const double c_ref = in.p_ref / (R * in.T_ref);
        data_[(col_T_ref * n) + i] = in.T_ref;
        data_[(col_R_ln_c_ref * n) + i] = R * std::log(c_ref);
        data_[(col_c_p * n) + i] = in.c_p;
        data_[(col_ln_T_ref * n) + i] = std::log(in.T_ref);
        data_[(col_h_ref * n) + i] = in.h_ref;
        data_[(col_s_ref * n) + i] = in.s_ref;
    }

    Storage data_{}; // Column-major (SoA) parameter storage.
};

static_assert(IdealEoS<ConstantCp<2>>, "ConstantCp must satisfy the IdealEoS concept.");
static_assert(IdealEoS<ConstantCp<std::dynamic_extent>>, "ConstantCp must satisfy the IdealEoS concept.");
} // namespace fugacity
