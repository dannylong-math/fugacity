#pragma once

#include <array>
#include <string_view>

// name, value function, d/dT function, d/dc function, composition-gradient function
#define SYNTHESIZE_TEST_PROPERTY_CATALOG(X)                                                                            \
    X(helmholtz, calc_helmholtz, calc_helmholtz_dT, calc_helmholtz_dc, calc_helmholtz_dx)                              \
    X(pressure, calc_pressure, calc_pressure_dT, calc_pressure_dc, calc_pressure_dx)                                   \
    X(internal_energy, calc_internal_energy, calc_internal_energy_dT, calc_internal_energy_dc,                         \
      calc_internal_energy_dx)                                                                                         \
    X(enthalpy, calc_enthalpy, calc_enthalpy_dT, calc_enthalpy_dc, calc_enthalpy_dx)                                   \
    X(entropy, calc_entropy, calc_entropy_dT, calc_entropy_dc, calc_entropy_dx)                                        \
    X(gibbs, calc_gibbs, calc_gibbs_dT, calc_gibbs_dc, calc_gibbs_dx)                                                  \
    X(dp_dc, calc_dp_dc, calc_dp_dc_dT, calc_dp_dc_dc, calc_dp_dc_dx)                                                  \
    X(dp_dT, calc_dp_dT, calc_dp_dT_dT, calc_dp_dT_dc, calc_dp_dT_dx)                                                  \
    X(cv, calc_cv, calc_cv_dT, calc_cv_dc, calc_cv_dx)                                                                 \
    X(cp, calc_cp, calc_cp_dT, calc_cp_dc, calc_cp_dx)

// Same shape as the scalar catalogue, but each function also takes an
// effective molar mass immediately after temperature.
#define SYNTHESIZE_TEST_MOLAR_MASS_PROPERTY_CATALOG(X)                                                                 \
    X(sound_speed_squared, calc_sound_speed_squared, calc_sound_speed_squared_dT, calc_sound_speed_squared_dc,         \
      calc_sound_speed_squared_dx)

namespace synthesize_test {

struct property_descriptor {
    std::string_view name;
};

#define SYNTHESIZE_TEST_DESCRIPTOR(NAME, VALUE, DT, DC, DX) property_descriptor{#NAME},
inline constexpr auto property_catalog = std::array{SYNTHESIZE_TEST_PROPERTY_CATALOG(SYNTHESIZE_TEST_DESCRIPTOR)};
#undef SYNTHESIZE_TEST_DESCRIPTOR

static_assert(property_catalog.size() == 10, "Update generated contract coverage when the public catalogue changes.");

#define SYNTHESIZE_TEST_SPECIAL_DESCRIPTOR(NAME, VALUE, DT, DC, DX) property_descriptor{#NAME},
inline constexpr auto molar_mass_property_catalog =
    std::array{SYNTHESIZE_TEST_MOLAR_MASS_PROPERTY_CATALOG(SYNTHESIZE_TEST_SPECIAL_DESCRIPTOR)};
#undef SYNTHESIZE_TEST_SPECIAL_DESCRIPTOR

static_assert(molar_mass_property_catalog.size() == 1);

} // namespace synthesize_test
