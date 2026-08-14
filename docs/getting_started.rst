Getting started
===============

Fugacity is a header-only library. Add the ``include`` directory to your
compiler's search path and include the umbrella header:

.. code-block:: cpp

   #include <fugacity/fugacity.hpp>

The library requires C++23, Clang, and an Enzyme installation built for the
same LLVM version as Clang. The repository's CMake configuration exposes the
``Fugacity::Fugacity`` interface target and provides presets for development
and testing.

Building the documentation
--------------------------

Install the Python documentation dependencies into the conventional local
environment, then use the documentation preset:

.. code-block:: console

   $ python3 -m venv .dependencies/docs
   $ .dependencies/docs/bin/pip install -r docs/requirements.txt
   $ cmake --preset docs
   $ cmake --build --preset docs

The generated site is written to ``build/docs/html/index.html``.
