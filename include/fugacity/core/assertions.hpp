#pragma once
///
/// Runtime precondition checks.
///

#include <stdexcept>
#include <string>

///
/// Check a programmer precondition in debug builds.
///
/// Throw ``std::logic_error`` when ``cond`` is false. The check is omitted when
/// ``NDEBUG`` is defined. The exception message includes the condition and source
/// location.
///
/// :param cond: Expression contextually convertible to ``bool``.
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
/// Require a strictly positive absolute temperature.
///
/// Throw ``std::domain_error`` when ``T <= 0``. This check remains enabled in
/// release builds.
///
/// :param T: Temperature [K].
///
/// \ingroup core
#define FUGACITY_REQUIRE_POSITIVE_TEMPERATURE(T)                                                                     \
    ((T) > 0 ? void(0) : throw std::domain_error("fugacity: temperature must be positive (T > 0 K required)"))
