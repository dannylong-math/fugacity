Getting started
===============

Requirements
------------

Fugacity requires:

* C++23;
* Clang;
* Enzyme built for the same LLVM version as Clang; and
* CMake 3.21 or newer when using the supplied build configuration.

Enzyme is a Clang plugin. Matching the Clang and LLVM versions is therefore a
runtime requirement of the compiler toolchain, not only a configuration detail.

Including Fugacity
------------------

Include the umbrella header to use the complete public API:

.. code-block:: cpp

   #include <fugacity/fugacity.hpp>

The CMake configuration exposes the ``Fugacity::Fugacity`` interface target.
Linking this target supplies the include directory, the C++23 requirement, and
the Enzyme compiler-plugin flags.

To use a source checkout from another CMake project:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.21)
   project(thermo_example LANGUAGES CXX)

   add_subdirectory(/path/to/fugacity fugacity-build)

   add_executable(thermo_example main.cpp)
   target_link_libraries(thermo_example PRIVATE Fugacity::Fugacity)

Configure the consuming project with Clang and tell CMake where to find Enzyme:

.. code-block:: console

   $ cmake -S . -B build \
       -D CMAKE_CXX_COMPILER=clang++ \
       -D Enzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
   $ cmake --build build

Building and testing Fugacity
-----------------------------

Configure a development build from the repository root:

.. code-block:: console

   $ cmake --preset debug \
       -D Enzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
   $ cmake --build --preset debug
   $ ctest --preset debug

The ``release`` and ``release-max`` presets build optimized versions. Use
``release-max`` only when the application guarantees finite inputs; its
floating-point options permit transformations that do not preserve NaN and
infinity behavior.

Next step
---------

Continue with :doc:`tutorial` to assemble an ideal and residual model, define a
thermodynamic state, and calculate properties.

Building the documentation
--------------------------

Install the Python dependencies and use the documentation preset:

.. code-block:: console

   $ python3 -m venv .dependencies/docs
   $ .dependencies/docs/bin/pip install -r docs/requirements.txt
   $ cmake --preset docs
   $ cmake --build --preset docs

The generated site is written to ``build/docs/html/index.html``.
