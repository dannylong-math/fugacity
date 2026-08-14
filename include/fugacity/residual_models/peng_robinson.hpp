#pragma once
///
/// File ``peng_robinson.hpp``.
/// Residual (departure) model for the Peng-Robinson equation of state.
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
/// Residual Helmholtz contribution of the Peng-Robinson equation of state.
///
/// A fugacity::BaseCubic model with :math:`\Delta_{1,2} = 1 \pm \sqrt{2}` and the
/// pure-species parameters built from each species' critical point and acentric
/// factor:
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
/// with the alpha-function slope from the acentric factor :math:`\omega`:
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
/// The Helmholtz kernels and mixture rules live in BaseCubic; this class only
/// performs the parameter transformation at construction.
///
/// Pair the model with an ideal contribution in a fugacity::EoS and evaluate
/// properties through the free functions in core_calculations.hpp:
///
/// .. code-block:: cpp
///
///    using namespace fugacity;
///
///    // CH4 and CO2 from critical data, with one interaction coefficient.
///    PengRobinson<2> residual(
///        std::array{PengRobinson<2>::SpeciesInput{.T_c = 190.564, .P_c = 4.5992e6, .omega = 0.011},
///                   PengRobinson<2>::SpeciesInput{.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394}},
///        std::array{0.0, 0.09,
///                   0.09, 0.0});
///    EoS eos{some_ideal_model, residual};
///
///    const std::array<double, 2> x{0.4, 0.6}; // mole fractions
///    const double p = calc_pressure(eos, 500.0, std::span<const double, 2>{x}, 300.0);
///
///
/// Use the ``std::dynamic_extent`` default (e.g. ``PengRobinson<>``) when the
/// number of species is only known at run time.
///
///
/// :tparam N: Component count, or ``std::dynamic_extent`` for a runtime size.
///
/// \ingroup residual-models
template<std::size_t N = std::dynamic_extent> class PengRobinson : public BaseCubic<PengRobinson<N>, N> {
public:
    /// Natural per-species input: critical point and acentric factor.
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
    /// Construct a compile-time-sized model from per-species critical data.
    ///
    /// Only available when the component count ``N`` is known at compile time.
    ///
    ///
    /// :param inputs: One SpeciesInput per species.
    /// :param kij: Full row-major :math:`N \times N` binary-interaction matrix
    ///               :math:`k_{ij}` [-] (entry ``kij[i*N + j]``); the diagonal must
    ///               be zero. Defaults to all zeros. May be asymmetric; only the
    ///               symmetric part affects the model.
    /// \id fixed-size
    ///
    explicit PengRobinson(const std::array<SpeciesInput, N>& inputs, const std::array<double, N * N>& kij = {})
        requires(N != std::dynamic_extent)
        : BaseCubic<PengRobinson<N>, N>(to_pure_all(inputs), kij)
    {
    }

    ///
    /// Construct a runtime-sized model from per-species critical data.
    ///
    /// Only available when ``N`` is ``std::dynamic_extent``; ``size()`` becomes
    /// ``inputs.size()``.
    ///
    ///
    /// :param inputs: One SpeciesInput per species.
    /// :param kij: Full row-major :math:`n \times n` binary-interaction matrix
    ///               :math:`k_{ij}` [-], or an empty span for all zeros (the
    ///               default). Size checked via ``assert``.
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
