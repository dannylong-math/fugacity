#pragma once

///
/// File ``numbers.hpp``.
/// Fundamental physical constants used throughout the library.
///

#include <concepts>
namespace fugacity {

///
/// The universal (molar) gas constant :math:`R`.
///
///
/// :tparam Number: Floating-point type the constant is expressed in.
///
/// Value: ``8.31446261815324`` (the full CODATA value).
/// Units: **J / (mol K)** (joule per mole per kelvin).
///
/// \ingroup core
template<std::floating_point Number = double> inline constexpr Number ideal_gas_constant{8.31446261815324};

} // namespace fugacity
