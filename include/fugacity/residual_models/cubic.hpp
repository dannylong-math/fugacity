#pragma once
///
/// Base class for generalized two-parameter cubic residual models.
///

#include "fugacity/core/assertions.hpp"
#include "fugacity/core/eos_base.hpp"
#include "fugacity/core/numbers.hpp"

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

namespace fugacity {

///
/// Implement the common Helmholtz equations and mixture rules for a generalized
/// two-parameter cubic residual model.
///
/// The molar residual Helmholtz energy is
///
/// .. math::
///
///      \frac{a_r}{RT} = \psi_1 - \frac{a_m}{RT}\,\psi_2,
///      \qquad
///      \psi_1 = -\ln(1 - b_m c),
///      \qquad
///      \psi_2 = \frac{\ln\!\left(\dfrac{\Delta_1 b_m c + 1}{\Delta_2 b_m c + 1}\right)}
///                    {b_m (\Delta_1 - \Delta_2)},
///
/// with the one-fluid mixture rules
///
/// .. math::
///
///      a_m = \sum_i \sum_j x_i x_j\, (1 - \bar{k}_{ij}) \sqrt{a_{ii}(T)\, a_{jj}(T)},
///      \qquad
///      b_m = \sum_i x_i b_{ii},\qquad
///      \bar{k}_{ij}=\frac{k_{ij}+k_{ji}}{2},
///
/// with the temperature-dependent pure attractive parameter
///
/// .. math::
///
///      a_{ii}(T) = a_{0,ii} \left[1 + m_{ii}\left(1 - \sqrt{T/T_{c,i}}\right)\right]^2
///                = a_{0,ii}\, \alpha_i(T)^2,
///
/// A derived class must define ``delta1`` and ``delta2``, with
/// ``delta1 != delta2``, and convert its public inputs to one
/// :cpp:class:`PureSpecies` record per component. See
/// :cpp:class:`fugacity::PengRobinson` for a concrete implementation. Use
/// :cpp:class:`fugacity::VanDerWaals` for the degenerate case
/// :math:`\Delta_1=\Delta_2=0`.
///
///
/// :tparam Derived: Concrete cubic model.
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \ingroup residual-models
template<class Derived, std::size_t N = std::dynamic_extent> class BaseCubic : public BaseEoS<N> {
public:
    ///
    /// Molar residual Helmholtz energy :math:`a_r = RT\psi_1 - a_m\psi_2`.
    ///
    /// :param c: Molar concentration [mol/m^3]. Must satisfy :math:`b_m c < 1`.
    /// :param x: Mole-fraction array [-].
    /// :param T: Temperature [K].
    /// :returns: Molar residual Helmholtz energy [J/mol].
    ///
    template<std::floating_point Number> [[nodiscard]] Number calc_helmholtz(Number c, const Number* x, Number T) const
    {
        static_assert(Derived::delta1 != Derived::delta2,
                      "BaseCubic requires delta1 != delta2; use VanDerWaals for delta1 == delta2 == 0.");
        constexpr double d1 = Derived::delta1;
        constexpr double d2 = Derived::delta2;
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        const Number sT = std::sqrt(T);

        Number am{0};
        Number bm{0};
        for (std::size_t i = 0; i < n; ++i) {
            bm += x[i] * b_[i];
            // |alpha_i|: a_ii = a0 alpha^2 >= 0, so sqrt(a_ii a_jj) carries |alpha|.
            const Number ti = x[i] * std::abs(p_[i] - (q_[i] * sT));
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += x[j] * std::abs(p_[j] - (q_[j] * sT)) * a_[(i * n) + j];
            }
            am += ti * row;
        }
        const Number bc = bm * c;
        const Number psi1 = -std::log(Number{1} - bc);
        const Number psi2 = std::log(((d1 * bc) + Number{1}) / ((d2 * bc) + Number{1})) / (bm * (d1 - d2));
        return (R * T * psi1) - (am * psi2);
    }

    ///
    /// Total residual Helmholtz energy density :math:`\Psi = c\,a_r`.
    ///
    /// Evaluated directly in partial concentrations: with
    ///
    /// :math:`B = \sum_i \rho_i b_{ii} = b_m c` and
    /// :math:`A = \sum_{ij} \rho_i \rho_j a_{ij}(T) = a_m c^2`,
    /// :math:`\Psi = -RTc\ln(1-B) - A
    /// \ln\!\left((\Delta_1 B + 1)/(\Delta_2 B + 1)\right)
    /// / \bigl(B(\Delta_1 - \Delta_2)\bigr)`, which avoids forming mole fractions.
    ///
    ///
    /// :param rho_i: Partial molar concentrations [mol/m^3].
    /// :param T: Temperature [K].
    /// :returns: Residual Helmholtz energy density [J/m^3].
    ///
    template<std::floating_point Number>
    [[nodiscard]] Number calc_helmholtz_density(const Number* rho_i, Number T) const
    {
        static_assert(Derived::delta1 != Derived::delta2,
                      "BaseCubic requires delta1 != delta2; use VanDerWaals for delta1 == delta2 == 0.");
        constexpr double d1 = Derived::delta1;
        constexpr double d2 = Derived::delta2;
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        const Number sT = std::sqrt(T);

        Number c{0};
        Number bc{0}; // b_m * c   = sum_i rho_i b_i
        Number ac{0}; // a_m * c^2 = sum_ij rho_i rho_j A_ij |alpha_i| |alpha_j|
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
            bc += rho_i[i] * b_[i];
            const Number ti = rho_i[i] * std::abs(p_[i] - (q_[i] * sT));
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += rho_i[j] * std::abs(p_[j] - (q_[j] * sT)) * a_[(i * n) + j];
            }
            ac += ti * row;
        }
        // Psi = c a_r: the c psi_2 factor folds into ac / (bc (d1 - d2)).
        const Number psi1 = -std::log(Number{1} - bc);
        return (R * T * c * psi1) -
               (ac * std::log(((d1 * bc) + Number{1}) / ((d2 * bc) + Number{1})) / (bc * (d1 - d2)));
    }

    ///
    /// Per-component residual Helmholtz energy density.
    ///
    /// The residual does not decompose naturally per component, so the
    /// mole-fraction-weighted convention :math:`\Psi_i = (\rho_i / c)\,\Psi` is
    /// used; it satisfies :math:`\sum_i \Psi_i = \Psi` by construction.
    ///
    ///
    /// :param rho_i: Partial molar concentrations [mol/m^3].
    /// :param T: Temperature [K].
    /// :param out: Per-component Helmholtz energy density [J/m^3]; length ``size()``.
    ///
    template<std::floating_point Number> void calc_partial_helmholtz(const Number* rho_i, Number T, Number* out) const
    {
        const std::size_t n = this->size();
        Number c{0};
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
        }
        // Mole-fraction-weighted decomposition: out[i] = (rho_i / c) * Psi.
        const Number scale = calc_helmholtz_density(rho_i, T) / c;
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = rho_i[i] * scale;
        }
    }

protected:
    ///
    /// Generalized cubic parameters for one species.
    struct PureSpecies {
        double a0;  ///< Attractive parameter at :math:`T = T_c` [J m^3/mol^2].
        double b;   ///< Covolume :math:`b_{ii}` [m^3/mol].
        double m;   ///< Alpha-function slope coefficient :math:`m_{ii}` [-].
        double T_c; ///< Critical temperature :math:`T_{c,i}` [K].
    };

    ///
    /// Construct the generalized cubic model.
    ///
    /// :param species: One PureSpecies record per species.
    /// :param kij: Full row-major :math:`n \times n` binary-interaction matrix
    ///                :math:`k_{ij}` [-], or an empty span for all zeros. May be
    ///                asymmetric; the model uses its symmetric part. Supply
    ///                exactly ``species.size() * species.size()`` entries.
    ///
    BaseCubic(std::span<const PureSpecies> species, std::span<const double> kij) : BaseEoS<N>(species.size())
    {
        const std::size_t n = species.size();
        FUGACITY_ASSERT(kij.empty() || kij.size() == n * n);
        if constexpr (N == std::dynamic_extent) {
            b_.resize(n);
            p_.resize(n);
            q_.resize(n);
            a_.resize(n * n);
        }
        for (std::size_t i = 0; i < n; ++i) {
            b_[i] = species[i].b;
            p_[i] = 1.0 + species[i].m;
            q_[i] = species[i].m / std::sqrt(species[i].T_c);
        }
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                // Only the symmetric part of k_ij can affect a_m = sum_ij x_i x_j a_ij.
                const double k_sym = kij.empty() ? 0.0 : 0.5 * (kij[(i * n) + j] + kij[(j * n) + i]);
                a_[(i * n) + j] = (1.0 - k_sym) * std::sqrt(species[i].a0 * species[j].a0);
            }
        }
    }

private:
    using Vec = std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, N>>;
    using Mat = std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, N * N>>;

    // alpha_i(T) = 1 + m_i (1 - sqrt(T/T_c,i)) is evaluated as p_i - q_i * sqrt(T),
    // so one sqrt per evaluation covers every species and every pair. alpha_i can be
    // negative at reasonable operating conditions (T well above T_c,i), so the cross
    // term sqrt(a_ii a_jj) = sqrt(a0_ii a0_jj) |alpha_i| |alpha_j| takes |alpha_i|
    // per species; the pair loop itself stays multiply-accumulate only.
    Vec b_{}; // b_ii covolumes [m^3/mol].
    Vec p_{}; // p_i = 1 + m_i [-].
    Vec q_{}; // q_i = m_i / sqrt(T_c,i) [K^-1/2].
    Mat a_{}; // A_ij = (1 - (k_ij + k_ji)/2) * sqrt(a0_ii * a0_jj), row-major [J m^3/mol^2],
              // so a_m = sum_ij x_i x_j A_ij |alpha_i| |alpha_j|.
};

} // namespace fugacity
