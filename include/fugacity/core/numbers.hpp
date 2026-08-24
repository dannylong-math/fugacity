#pragma once

///
/// Physical constants.
///

#include <concepts>
namespace fugacity {

///
/// Universal molar gas constant \f$R = 8.31446261815324\f$ [J/(mol K)].
///
/// \tparam Number Floating-point type used for the value.
///
/// \ingroup core
template<std::floating_point Number = double> inline constexpr Number ideal_gas_constant{8.31446261815324};

} // namespace fugacity
