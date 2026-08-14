#pragma once
///
/// Peng-Robinson residual model.
///

#include "fugacity/core/concepts.hpp"
#include "fugacity/core/numbers.hpp"
#include "fugacity/residual_models/cubic.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

namespace fugacity {

///
/// Peng-Robinson residual model.
///
/// The molar residual Helmholtz energy is
///
/// .. math::
///
///    a^\mathrm{res}
///    =-RT\ln(1-b_m c)
///    -\frac{a_m}{b_m(\Delta_1-\Delta_2)}
///    \ln\!\left(\frac{1+\Delta_1b_mc}{1+\Delta_2b_mc}\right),
///    \qquad \Delta_{1,2}=1\mathbin{\pm}\sqrt{2}.
///
/// The mixture parameters are
///
/// .. math::
///
///    a_m=\sum_i\sum_jx_ix_j(1-\bar{k}_{ij})
///        \sqrt{a_{ii}(T)a_{jj}(T)},\qquad
///    b_m=\sum_i x_i b_{ii},\qquad
///    \bar{k}_{ij}=\frac{k_{ij}+k_{ji}}{2},
///
/// .. math::
///
///    a_{ii}(T)=a_{0,ii}
///      \left[1+m_{ii}\left(1-\sqrt{T/T_{c,i}}\right)\right]^2.
///
/// Compute the pure-species parameters from critical properties as
///
/// .. math::
///
///      \eta_c = \left(1 + \sqrt[3]{4 - \sqrt 8} + \sqrt[3]{4 + \sqrt 8}\right)^{-1},
///      \qquad
///      \Omega_a = \frac{8 + 40\eta_c}{49 - 37\eta_c},
///      \qquad
///      \Omega_b = \frac{\eta_c}{3 + \eta_c},
///
/// .. math::
///
///      a_{0,ii} = \Omega_a \frac{(R T_c)^2}{P_c},
///      \qquad
///      b_{ii} = \Omega_b \frac{R T_c}{P_c},
///
/// and compute the alpha-function coefficient from the acentric factor:
///
/// .. math::
///
///      m_{ii} = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^2
///      \quad \text{for } \omega \le 0.491,
///
///      m_{ii} = 0.379642 + 1.48503\,\omega - 0.164423\,\omega^2
///      + 0.016666\,\omega^3
///      \quad \text{for } \omega > 0.491.
///
/// .. code-block:: cpp
///
///    using PR = fugacity::PengRobinson<2>;
///    const std::array<PR::SpeciesInput, 2> species{{
///        {.T_c = 190.564, .P_c = 4.5992e6, .omega = 0.011},
///        {.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394},
///    }};
///    const std::array<double, 4> kij{0.0, 0.09,
///                                     0.09, 0.0};
///    const PR residual{species, kij};
///
///    const std::array<double, 2> x{0.4, 0.6};
///    const double a_res = residual.calc_helmholtz(500.0, x.data(), 300.0);
///
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime count.
///
/// \ingroup residual-models
template<std::size_t N = std::dynamic_extent> class PengRobinson : public BaseCubic<PengRobinson<N>, N> {
public:
    /// Critical properties for one species.
    struct SpeciesInput {
        double T_c;   ///< Critical temperature :math:`T_c` [K].
        double P_c;   ///< Critical pressure :math:`P_c` [Pa].
        double omega; ///< Acentric factor :math:`\omega` [-].
    };

    /// Generalized-cubic constant :math:`\Delta_1 = 1 + \sqrt 2`.
    static constexpr double delta1 = 1.0 + std::numbers::sqrt2;
    /// Generalized-cubic constant :math:`\Delta_2 = 1 - \sqrt 2`.
    static constexpr double delta2 = 1.0 - std::numbers::sqrt2;

    ///
    /// Construct a fixed-size model from critical properties.
    ///
    /// :param inputs: One :cpp:class:`SpeciesInput` per species.
    /// :param kij: Full row-major :math:`N \times N` binary-interaction matrix
    ///               :math:`k_{ij}` [-], stored as ``kij[i*N + j]``. The default
    ///               matrix is zero. The model uses its symmetric part.
    /// \id fixed-size
    ///
    explicit PengRobinson(const std::array<SpeciesInput, N>& inputs, const std::array<double, N * N>& kij = {})
        requires(N != std::dynamic_extent)
        : BaseCubic<PengRobinson<N>, N>(to_pure_all(inputs), kij)
    {
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
    explicit PengRobinson(std::span<const SpeciesInput> inputs, std::span<const double> kij = {})
        requires(N == std::dynamic_extent)
        : BaseCubic<PengRobinson<N>, N>(to_pure_all(inputs), kij)
    {
    }

private:
    using Base = BaseCubic<PengRobinson<N>, N>;
    using Pure = typename Base::PureSpecies;

    // Map one species' critical data to the generalized cubic parameters.
    static Pure to_pure(const SpeciesInput& in)
    {
        constexpr double R = ideal_gas_constant<double>;
        const double s8 = std::sqrt(8.0);
        const double eta_c = 1.0 / (1.0 + std::cbrt(4.0 - s8) + std::cbrt(4.0 + s8));
        const double omega_a = (8.0 + (40.0 * eta_c)) / (49.0 - (37.0 * eta_c));
        const double omega_b = eta_c / (3.0 + eta_c);
        const double w = in.omega;
        const double m = w <= 0.491 ? 0.37464 + (1.54226 * w) - (0.26992 * w * w)
                                    : 0.379642 + (1.48503 * w) - (0.164423 * w * w) + (0.016666 * w * w * w);
        const double RTc = R * in.T_c;
        return {.a0 = omega_a * RTc * RTc / in.P_c, .b = omega_b * RTc / in.P_c, .m = m, .T_c = in.T_c};
    }

    static std::array<Pure, N> to_pure_all(const std::array<SpeciesInput, N>& inputs)
        requires(N != std::dynamic_extent)
    {
        std::array<Pure, N> out{};
        for (std::size_t i = 0; i < N; ++i) {
            out[i] = to_pure(inputs[i]);
        }
        return out;
    }

    static std::vector<Pure> to_pure_all(std::span<const SpeciesInput> inputs)
        requires(N == std::dynamic_extent)
    {
        std::vector<Pure> out;
        out.reserve(inputs.size());
        for (const SpeciesInput& in : inputs) {
            out.push_back(to_pure(in));
        }
        return out;
    }
};

static_assert(ResidualEoS<PengRobinson<2>>, "PengRobinson must satisfy the ResidualEoS concept.");
static_assert(ResidualEoS<PengRobinson<std::dynamic_extent>>, "PengRobinson must satisfy the ResidualEoS concept.");

} // namespace fugacity
