C++ API reference
=================

The API operates on molar concentration ``c`` [mol/m^3], mole fractions ``x``
[-], partial molar concentrations ``rho_i`` [mol/m^3], and temperature ``T``
[K]. Scalar property functions accept either a contiguous container or a
pointer to the composition array. Vector-valued functions use ``std::span``.

Functions ending in ``_dT``, ``_dc``, and ``_dx`` return derivatives with
respect to temperature, molar concentration, and mole fractions. Composition
gradients treat the mole fractions as independent coordinates; see
:ref:`composition-derivatives` for derivatives constrained to the composition
simplex.

Core API
--------

.. toctree::
   :maxdepth: 1

   concepts

.. cpp-apigen-group:: core

Ideal-gas models
----------------

.. cpp-apigen-group:: ideal-models

Residual models
---------------

.. cpp-apigen-group:: residual-models
