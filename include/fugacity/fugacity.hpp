#pragma once

///
/// File ``fugacity.hpp``.
/// Umbrella header for Fugacity; including it pulls in the public API.
///
/// This is the single header in the root ``include/fugacity`` directory. It
/// aggregates the core machinery (the EoS pair, concepts, and thermodynamic
/// property calculations) together with the bundled models so that a translation
/// unit can access the whole library with one include:
///
/// .. code-block:: cpp
///
///    #include "fugacity/fugacity.hpp"
///
///
/// Individual components remain available under their respective subdirectories
/// (``core/``, ``ideal_models/``, ``residual_models/``) for finer-grained includes.
///

#include "fugacity/core/core_calculations.hpp"
#include "fugacity/core/xlnx.hpp"
#include "fugacity/ideal_models/const_cp.hpp"
#include "fugacity/ideal_models/nasa7.hpp"
#include "fugacity/ideal_models/nasa9.hpp"
#include "fugacity/residual_models/no_residual.hpp"
#include "fugacity/residual_models/peng_robinson.hpp"
#include "fugacity/residual_models/van_der_waals.hpp"
