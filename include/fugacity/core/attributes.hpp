#pragma once
///
/// Compiler attributes used by differentiated kernels.
///

///
/// Force Clang to inline a function used inside a differentiated Helmholtz
/// kernel.
///
/// Annotate small helper functions called by a model's Helmholtz calculations.
/// The macro expands to ``[[clang::always_inline]]`` under Clang and to nothing
/// for documentation or unsupported compilers.
///
#if defined(FUGACITY_DOCUMENTATION)
#define FUGACITY_ALWAYS_INLINE
#elif defined(__clang__)
#define FUGACITY_ALWAYS_INLINE [[clang::always_inline]]
#else
#define FUGACITY_ALWAYS_INLINE
#endif
