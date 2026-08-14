Equation-of-state concepts
==========================

The model concepts define the Helmholtz interface and distinguish ideal and
residual contributions. The calculation functions are templates over the
floating-point number type so Enzyme can differentiate the concrete model.

.. cpp:concept:: template<class E> EquationOfState

   Require ``size()`` and three mutually consistent Helmholtz calculations:

   .. code-block:: cpp

      template<std::floating_point Number>
      Number calc_helmholtz(Number c, const Number* x, Number T) const;

      template<std::floating_point Number>
      Number calc_helmholtz_density(const Number* rho_i, Number T) const;

      template<std::floating_point Number>
      void calc_partial_helmholtz(
          const Number* rho_i, Number T, Number* out) const;

   ``calc_helmholtz`` returns molar Helmholtz energy [J/mol]. The density
   functions return total or per-component Helmholtz energy density [J/m^3].
   They must satisfy

   .. math::

      \Psi(\boldsymbol{\rho},T)=c\,a(c,\boldsymbol{x},T),\qquad
      c=\sum_i\rho_i,\qquad x_i=\rho_i/c,

   and

   .. math::

      \Psi=\sum_i\Psi_i.

.. cpp:concept:: template<class E> IdealEoS

   Match an :cpp:concept:`EquationOfState` that publicly derives from
   :cpp:class:`fugacity::BaseIdealEoS`.

.. cpp:concept:: template<class E> ResidualEoS

   Match an :cpp:concept:`EquationOfState` that does not derive from
   :cpp:class:`fugacity::BaseIdealEoS`.
