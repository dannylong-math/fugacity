#pragma once
///
/// File ``assertions.hpp``.
/// Debug-only precondition macro that throws instead of aborting.
///

#include <stdexcept>
#include <string>

///
/// Macro ``FUGACITY_ASSERT``.
/// Check a precondition, throwing ``std::logic_error`` if it does not hold.
///
/// Like the standard ``assert``, this validates a programmer-side precondition and
/// is compiled out entirely in release builds (``NDEBUG`` defined), so it costs
/// nothing on the hot numerical paths. Unlike ``assert``, a failing check in a
/// debug build **throws** a ``std::logic_error`` rather than calling ``abort``.
///
/// Throwing (rather than aborting) is what makes the failure path testable: a
/// test can assert the throw with boost::ut's ``throws`` matcher, and — because
/// the stack unwinds normally instead of the process dying — coverage
/// instrumentation still flushes, so the failure branch is recorded.
///
/// The thrown message carries the stringized condition and the source location.
///
///
/// :param cond: The precondition expression; must be contextually convertible to
///             ``bool``.
///
/// \ingroup core
#ifdef NDEBUG
#define FUGACITY_ASSERT(cond) ((void)0)
#else
#define FUGACITY_ASSERT(cond)                                                                                        \
    ((cond) ? void(0)                                                                                                  \
            : throw std::logic_error(std::string("FUGACITY_ASSERT failed: " #cond " (" __FILE__ ":") +               \
                                     std::to_string(__LINE__) + ")"))
#endif

///
/// Macro ``FUGACITY_REQUIRE_POSITIVE_TEMPERATURE``.
/// Validate that a temperature is strictly positive, throwing
/// ``std::domain_error`` otherwise.
///
/// Unlike FUGACITY_ASSERT this is an **always-on** runtime check (it is present
/// even in release builds): a non-positive absolute temperature is unphysical
/// input data, not a programmer-side logic error, so it is validated
/// unconditionally. The single floating-point comparison is branch-predictable
/// and negligible on the hot path.
///
/// The throw is textually inlined at each call site, which the Enzyme autodiff
/// pass differentiates cleanly (the throw is dead control flow with no derivative
/// contribution).
///
///
/// :param T: The temperature [K]; must be contextually comparable to ``0``.
///
/// \ingroup core
#define FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T)                                                                     \
    ((T) > 0 ? void(0) : throw std::domain_error("fugacity: temperature must be positive (T > 0 K required)"))
