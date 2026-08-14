# Fugacity

[![CI test](https://github.com/dannylong-math/fugacity/actions/workflows/ci.yml/badge.svg)](https://github.com/dannylong-math/fugacity/actions/workflows/ci.yml)
[![Documentation](https://github.com/dannylong-math/fugacity/actions/workflows/docs.yml/badge.svg)](https://dannylong-math.github.io/fugacity/)
[![codecov](https://codecov.io/gh/dannylong-math/fugacity/graph/badge.svg)](https://codecov.io/gh/dannylong-math/fugacity)

Fugacity is a header-only C++23 library for calculating thermodynamic
properties from Helmholtz-energy equations of state. Pair an ideal contribution
with a residual contribution; the library uses
[Enzyme](https://enzyme.mit.edu) automatic differentiation to calculate
pressure, caloric properties, chemical potentials, fugacities, and property
derivatives.

The bundled models are:

- ideal: constant heat capacity, NASA-7, and NASA-9;
- residual: zero residual, van der Waals, and Peng-Robinson.

All public thermodynamic inputs and outputs use SI units. The density variable is
molar concentration `c` [mol/m^3].

## Requirements

- C++23
- Clang
- Enzyme built for the same LLVM version as Clang
- CMake 3.21 or newer when using the supplied build configuration

On Linux or macOS, Homebrew can install a matching LLVM and Enzyme toolchain:

```sh
brew install llvm enzyme lld
export PATH="$(brew --prefix llvm)/bin:$PATH"
```

## Minimal example

```cpp
#include <fugacity/fugacity.hpp>

#include <array>
#include <iostream>

int main()
{
    using Ideal = fugacity::ConstantCp<1>;
    const std::array<Ideal::SpeciesInput, 1> species{{{
        .T_ref = 298.15, // K
        .p_ref = 1.0e5,  // Pa
        .c_p = 29.12,    // J/(mol K)
        .h_ref = 0.0,    // J/mol
        .s_ref = 191.61, // J/(mol K)
    }}};

    const fugacity::EoS eos{
        Ideal{species}, fugacity::NoResidual<1>{}};

    const double T = 350.0;           // K
    const double c = 40.0;            // mol/m^3
    const std::array<double, 1> x{1.0}; // mole fraction

    const double p = fugacity::calc_pressure(eos, c, x, T); // Pa
    const double h = fugacity::calc_enthalpy(eos, c, x, T); // J/mol

    std::cout << "p = " << p << " Pa\n"
              << "h = " << h << " J/mol\n";
}
```

Link the interface target from CMake:

```cmake
add_subdirectory(/path/to/fugacity fugacity-build)
target_link_libraries(your_target PRIVATE Fugacity::Fugacity)
```

Configure the consuming project with Clang and the Enzyme package location:

```sh
cmake -S . -B build \
  -D CMAKE_CXX_COMPILER=clang++ \
  -D Enzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
cmake --build build
```

## Documentation

- [Getting started](https://dannylong-math.github.io/fugacity/getting_started.html)
- [Model and property tutorial](https://dannylong-math.github.io/fugacity/tutorial.html)
- [Implementing a new model](https://dannylong-math.github.io/fugacity/implementing_a_new_eos.html)
- [C++ API reference](https://dannylong-math.github.io/fugacity/api/index.html)

Build the documentation locally with:

```sh
python3 -m venv .dependencies/docs
.dependencies/docs/bin/pip install -r docs/requirements.txt
cmake --preset docs
cmake --build --preset docs
```

## Building and testing the repository

```sh
cmake --preset debug \
  -D Enzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
cmake --build --preset debug
ctest --preset debug
```

The repository also provides `release`, `release-max`, `coverage`, and
`debug-tidy` presets.

## License

See [LICENSE](LICENSE).
