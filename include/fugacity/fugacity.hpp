#pragma once

///
/// Include the complete public API.
///
/// .. code-block:: cpp
///
///    #include <fugacity/fugacity.hpp>
///
/// Include individual headers from ``core/``, ``ideal_models/``, or
/// ``residual_models/`` when the umbrella header is not required.
///

#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/xlnx.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/ideal_models/nasa7.hpp"
#include "fugacity/ideal_models/nasa9.hpp"
#include "fugacity/residual_models/no_residual.hpp"
#include "fugacity/residual_models/peng_robinson.hpp"
#include "fugacity/residual_models/van_der_waals.hpp"
