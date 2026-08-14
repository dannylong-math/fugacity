Implementing a new equation of state
====================================

A complete equation of state pairs an ideal contribution with a residual
(departure) contribution using :cpp:class:`fugacity::EoS`. Property functions
then differentiate the reduced molar Helmholtz energy

.. math::

   \alpha(c, \boldsymbol{x}, T) = \frac{a(c, \boldsymbol{x}, T)}{RT}

with Enzyme. A model supplies the Helmholtz energy; it does not implement those
derivatives itself.

Model interface
---------------

Every model must satisfy :cpp:concept:`EquationOfState`. It provides
three templated calculations and a component count:

.. code-block:: cpp

   template<std::floating_point Number>
   Number calc_helmholtz(Number c, const Number* x, Number T) const;

   template<std::floating_point Number>
   Number calc_helmholtz_density(const Number* rho_i, Number T) const;

   template<std::floating_point Number>
   void calc_partial_helmholtz(const Number* rho_i, Number T, Number* out) const;

   std::size_t size() const;

``calc_helmholtz`` returns molar energy in J/mol. The density functions return
energy per volume in J/m³ and must remain consistent with

.. math::

   \Psi(\boldsymbol{\rho}, T)
   = c\,a(c, \boldsymbol{x}, T), \qquad
   c = \sum_i \rho_i, \qquad
   x_i = \rho_i / c.

Use :cpp:class:`fugacity::BaseEoS` for compile-time or runtime component-count
storage. Ideal models additionally derive from
:cpp:class:`fugacity::BaseIdealEoS`; residual models do not.

Adding a model
--------------

#. Choose whether the model is ideal or residual and derive from the
   appropriate base classes.
#. Define a ``SpeciesInput`` record for the model's natural per-species data.
#. Provide constructors for fixed-size and, where useful, runtime-size models.
#. Implement the three Helmholtz-energy functions above using a shared
   formulation so their units and values agree.
#. Include the new public header from ``fugacity/fugacity.hpp``.
#. Add focused tests for construction, the equation-of-state concepts, energy
   consistency, and representative thermodynamic properties.

The existing :cpp:class:`fugacity::ConstantCp` and
:cpp:class:`fugacity::NoResidual` classes are compact examples. The
:cpp:class:`fugacity::PengRobinson` and :cpp:class:`fugacity::VanDerWaals`
classes demonstrate mixture models with binary-interaction parameters.

Automatic differentiation constraint
-------------------------------------

Fugacity is compiled with the Enzyme Clang plugin. Keep differentiated kernels
visible to the compiler, avoid opaque calls that Enzyme cannot differentiate,
and use ``FUGACITY_ALWAYS_INLINE`` on small helper kernels where the surrounding
implementation requires inlining. Run the full tests with a matching Clang and
Enzyme pair after adding or changing a model.
