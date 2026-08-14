#pragma once
///
/// van der Waals residual model.
///

#include "fugacity/core/assertions.hpp"
#include "fugacity/core/concepts.hpp"
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
/// van der Waals residual model.
///
/// The molar residual Helmholtz energy is
///
/// .. math::
///
///      a_r = -R T \ln(1 - b_m c) - a_m c,
///
/// with the one-fluid mixture rules
///
/// .. math::
///
///      a_m = \sum_i \sum_j x_i x_j\, (1 - \bar{k}_{ij}) \sqrt{a_{0,ii}\, a_{0,jj}},
///      \qquad
///      b_m = \sum_i x_i b_{ii},\qquad
///      \bar{k}_{ij}=\frac{k_{ij}+k_{ji}}{2},
///
/// with the pure-species parameters built from the critical point:
///
/// .. math::
///
///      a_{0,ii} = \frac{27 (R T_c)^2}{64 P_c}, \qquad b_{ii} = \frac{R T_c}{8 P_c}.
///
/// .. code-block:: cpp
///
///    using VdW = fugacity::VanDerWaals<2>;
///    const std::array<VdW::SpeciesInput, 2> species{{
///        {.T_c = 126.192, .P_c = 3.3958e6},
///        {.T_c = 304.1282, .P_c = 7.3773e6},
///    }};
///    const std::array<double, 4> kij{0.0, 0.05,
///                                     0.05, 0.0};
///    const VdW residual{species, kij};
///
///    const std::array<double, 2> x{0.4, 0.6};
///    const double a_res = residual.calc_helmholtz(500.0, x.data(), 300.0);
///
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \ingroup residual-models
template<std::size_t N = std::dynamic_extent> class VanDerWaals : public BaseEoS<N> {
public:
    /// Critical properties for one species.
    struct SpeciesInput {
        double T_c; ///< Critical temperature :math:`T_c` [K].
        double P_c; ///< Critical pressure :math:`P_c` [Pa].
    };

    ///
    /// Construct a fixed-size model from critical properties.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// :param kij: Full row-major :math:`N \times N` binary-interaction matrix
    ///               :math:`k_{ij}` [-], stored as ``kij[i*N + j]``. The default
    ///               matrix is zero. The model uses its symmetric part.
    /// \id fixed-size
    ///
    explicit VanDerWaals(const std::array<SpeciesInput, N>& inputs, const std::array<double, N * N>& kij = {})
        requires(N != std::dynamic_extent)
    {
        init(inputs, kij);
    }

    ///
    /// Construct a runtime-size model from critical properties.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// :param kij: Full row-major :math:`n \times n` binary-interaction matrix
    ///               :math:`k_{ij}` [-], or an empty span for a zero matrix.
    ///               Supply exactly ``inputs.size() * inputs.size()`` entries.
    /// \id runtime-size
    ///
    explicit VanDerWaals(std::span<const SpeciesInput> inputs, std::span<const double> kij = {})
        requires(N == std::dynamic_extent)
        : BaseEoS<N>(inputs.size())
    {
        b_.resize(inputs.size());
        a_.resize(inputs.size() * inputs.size());
        init(inputs, kij);
    }

    ///
    /// Molar residual Helmholtz energy :math:`a_r = -RT\ln(1-b_m c) - a_m c`.
    ///
    /// :param c: Molar concentration [mol/m^3]. Must satisfy :math:`b_m c < 1`.
    /// :param x: Mole-fraction array [-].
    /// :param T: Temperature [K].
    /// :returns: Molar residual Helmholtz energy [J/mol].
    ///
    template<std::floating_point Number> [[nodiscard]] Number calc_helmholtz(Number c, const Number* x, Number T) const
    {
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        Number am{0};
        Number bm{0};
        for (std::size_t i = 0; i < n; ++i) {
            bm += x[i] * b_[i];
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += x[j] * a_[(i * n) + j];
            }
            am += x[i] * row;
        }
        // a_r = R T psi_1 - a_m psi_2 with psi_1 = -ln(1 - b_m c), psi_2 = c (vdW).
        return (-R * T * std::log(Number{1} - (bm * c))) - (am * c);
    }

    ///
    /// Total residual Helmholtz energy density :math:`\Psi = c\,a_r`.
    ///
    /// Evaluated directly in partial concentrations as
    ///
    /// :math:`\Psi = -RTc\ln(1 - \sum_i \rho_i b_{ii}) - \sum_{ij}\rho_i\rho_j a_{ij}`,
    /// which avoids forming mole fractions.
    ///
    ///
    /// :param rho_i: Partial molar concentrations [mol/m^3].
    /// :param T: Temperature [K].
    /// :returns: Residual Helmholtz energy density [J/m^3].
    ///
    template<std::floating_point Number>
    [[nodiscard]] Number calc_helmholtz_density(const Number* rho_i, Number T) const
    {
        const Number R = ideal_gas_constant<Number>;
        const std::size_t n = this->size();
        Number c{0};
        Number bc{0}; // b_m * c = sum_i rho_i b_i
        Number ac{0}; // a_m * c^2 = sum_ij rho_i rho_j a_ij
        for (std::size_t i = 0; i < n; ++i) {
            c += rho_i[i];
            bc += rho_i[i] * b_[i];
            Number row{0};
            for (std::size_t j = 0; j < n; ++j) {
                row += rho_i[j] * a_[(i * n) + j];
            }
            ac += rho_i[i] * row;
        }
        return (-R * T * c * std::log(Number{1} - bc)) - ac;
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

private:
    using Vec = std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, N>>;
    using Mat = std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, N * N>>;

    // Build the stored parameters from critical data. Cold path: the pairwise
    // square roots are paid once here so the kernels never take any.
    void init(std::span<const SpeciesInput> inputs, std::span<const double> kij)
    {
        constexpr double R = ideal_gas_constant<double>;
        const std::size_t n = inputs.size();
        FUGACITY_ASSERT(kij.empty() || kij.size() == n * n);

        for (std::size_t i = 0; i < n; ++i) {
            b_[i] = R * inputs[i].T_c / (8.0 * inputs[i].P_c);
        }
        const auto a0 = [&](std::size_t i) {
            const double RTc = R * inputs[i].T_c;
            return 27.0 * RTc * RTc / (64.0 * inputs[i].P_c);
        };
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                // Only the symmetric part of k_ij can affect a_m = sum_ij x_i x_j a_ij.
                const double k_sym = kij.empty() ? 0.0 : 0.5 * (kij[(i * n) + j] + kij[(j * n) + i]);
                a_[(i * n) + j] = (1.0 - k_sym) * std::sqrt(a0(i) * a0(j));
            }
        }
    }

    Vec b_{}; // b_ii covolumes [m^3/mol].
    Mat a_{}; // a_ij = (1 - (k_ij + k_ji)/2) * sqrt(a0_ii * a0_jj), row-major [J m^3/mol^2].
              // T-independent for vdW (m_ii = 0), so fully precomputed at construction.
};

static_assert(ResidualEoS<VanDerWaals<2>>, "VanDerWaals must satisfy the ResidualEoS concept.");
static_assert(ResidualEoS<VanDerWaals<std::dynamic_extent>>, "VanDerWaals must satisfy the ResidualEoS concept.");

} // namespace fugacity
