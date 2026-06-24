# Synthesize

[![CI test](https://github.com/Miscibility/synthesize/actions/workflows/ci.yml/badge.svg)](https://github.com/Miscibility/synthesize/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/Miscibility/synthesize/graph/badge.svg)](https://codecov.io/gh/Miscibility/synthesize)

A C++ library implementing equations of state in an efficient way.

Synthesize is a header-only library. Thermodynamic properties are derived from
the reduced molar Helmholtz energy and its derivatives, which are obtained by
automatic differentiation with [Enzyme](https://enzyme.mit.edu). Because Enzyme
is an LLVM/Clang plugin, the project **must be built with Clang**.

## Requirements

- A Clang compiler matching the LLVM version Enzyme was built against
- [Enzyme](https://enzyme.mit.edu)
- CMake ≥ 3.21

### Installing Clang and Enzyme (the easy way)

On linux or Mac, the easiest way to install [Enzyme](https://enzyme.mit.edu) is
with [Homebrew](https://brew.sh):

```sh
brew install llvm enzyme lld
# Put the Homebrew Clang ahead of the system compiler:
export PATH="$(brew --prefix llvm)/bin:$PATH"
```

You can install LLVM/Enzyme other ways. Just ensure that CMake can find `clang++`.

## Building and testing

Configure with one of the CMake presets, passing the location of Enzyme's CMake
package:

```sh
# With the Brew installation: -D Enzyme_DIR="$(brew --prefix enzyme)/lib/cmake/Enzyme"
cmake --preset release -D Enzyme_DIR=/path/to/Enzyme/install/lib/cmake/Enzyme
cmake --build build/release
ctest --preset release
```

Available presets:

| Preset        | Purpose                                                        |
| ------------- | ------------------------------------------------------------- |
| `debug`       | `-O1 -g` with AddressSanitizer (Enzyme needs ≥ `-O1`)         |
| `release`     | `-O3 -march=native`                                           |
| `release-max` | `release` plus aggressive fast-math flags                     |
| `coverage`    | Clang source-based code coverage (see below)                  |
| `debug-tidy`  | `debug` with clang-tidy                                        |


## License

See [LICENSE](LICENSE).
