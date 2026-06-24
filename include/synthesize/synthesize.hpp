#pragma once

/**
 * @file synthesize.hpp
 * @brief Umbrella header for Synthesize; including it pulls in the public API.
 *
 * This is the single header in the root @c include/synthesize directory. It
 * aggregates the core machinery (the EoS pair, concepts, and thermodynamic
 * property calculations) together with the bundled models so that a translation
 * unit can access the whole library with one include:
 *
 * @code
 * #include "synthesize/synthesize.hpp"
 * @endcode
 *
 * Individual components remain available under their respective subdirectories
 * (@c core/, @c ideal_models/, @c residual_models/) for finer-grained includes.
 */

#include "synthesize/core/core_calculations.hpp"
#include "synthesize/core/xlnx.hpp"
#include "synthesize/ideal_models/const_cp.hpp"
#include "synthesize/ideal_models/nasa7.hpp"
#include "synthesize/ideal_models/nasa9.hpp"
#include "synthesize/residual_models/no_residual.hpp"
#include "synthesize/residual_models/peng_robinson.hpp"
#include "synthesize/residual_models/van_der_waals.hpp"
