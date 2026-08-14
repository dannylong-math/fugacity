Equation-of-state concepts
==========================

The concepts below describe the common model interface and distinguish ideal
contributions from residual contributions. Sphinx-Immaterial's current
libclang generator does not emit C++ concept declarations, so they are included
explicitly alongside the generated reference.

.. cpp:concept:: template<class E> EquationOfState

   Requires a model to provide ``size()``, molar Helmholtz energy,
   total Helmholtz-energy density, and per-component Helmholtz-energy density.
   The calculation members are templated by each concrete model so Fugacity can
   instantiate them with the number type used by automatic differentiation.

.. cpp:concept:: template<class E> IdealEoS

   An :cpp:concept:`EquationOfState` that publicly derives from
   :cpp:class:`BaseIdealEoS`.

.. cpp:concept:: template<class E> ResidualEoS

   An :cpp:concept:`EquationOfState` that does not derive from
   :cpp:class:`BaseIdealEoS`.
