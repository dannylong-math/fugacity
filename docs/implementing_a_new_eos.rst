Implementing a new equation-of-state model
==========================================

This guide implements a readable NASA-7 ideal model and tests it against the
NASA-7 equations. The bundled :cpp:class:`fugacity::Nasa7` follows the same
interface but uses precomputed, structure-of-arrays storage in its evaluation
kernels.

In Fugacity, a model is one contribution to a complete equation of state. An
ideal model and a residual model are combined with
:cpp:class:`fugacity::EoS`:

.. math::

   a(c,\boldsymbol{x},T)
   =a^\mathrm{ideal}(c,\boldsymbol{x},T)
   +a^\mathrm{res}(c,\boldsymbol{x},T).

You implement the Helmholtz contribution. The public property functions obtain
the required derivatives with Enzyme.

The model contract
------------------

Every model must satisfy :cpp:concept:`EquationOfState` and provide:

.. code-block:: cpp

   template<std::floating_point Number>
   Number calc_helmholtz(
       Number c, const Number* x, Number T) const;

   template<std::floating_point Number>
   Number calc_helmholtz_density(
       const Number* rho_i, Number T) const;

   template<std::floating_point Number>
   void calc_partial_helmholtz(
       const Number* rho_i, Number T, Number* out) const;

   std::size_t size() const;

The quantities and units are:

.. list-table:: Model interface
   :header-rows: 1
   :widths: 30 20 20 30

   * - Quantity
     - Symbol
     - Units
     - Array length
   * - Total molar concentration
     - :math:`c`
     - mol/m^3
     - scalar
   * - Mole fraction
     - :math:`x_i`
     - -
     - ``size()``
   * - Partial molar concentration
     - :math:`\rho_i`
     - mol/m^3
     - ``size()``
   * - Temperature
     - :math:`T`
     - K
     - scalar
   * - Molar Helmholtz contribution
     - :math:`a`
     - J/mol
     - scalar
   * - Helmholtz energy density
     - :math:`\Psi`
     - J/m^3
     - scalar
   * - Per-component density contribution
     - :math:`\Psi_i`
     - J/m^3
     - ``size()``

The three calculations must satisfy

.. math::

   \Psi(\boldsymbol{\rho},T)=c\,a(c,\boldsymbol{x},T),\qquad
   c=\sum_i\rho_i,\qquad x_i=\rho_i/c,

and

.. math::

   \Psi=\sum_i\Psi_i.

``calc_helmholtz_density`` must return the scalar density directly. Fugacity
differentiates this function with respect to every :math:`\rho_i` to calculate
chemical potentials.

Step 1: classify the model
--------------------------

Derive every model from :cpp:class:`fugacity::BaseEoS`. Supply the component
count as its template argument.

An ideal model must also derive publicly from
:cpp:class:`fugacity::BaseIdealEoS`:

.. code-block:: cpp

   template<std::size_t N = std::dynamic_extent>
   class ExampleNasa7
       : public fugacity::BaseEoS<N>,
         public fugacity::BaseIdealEoS {
       // ...
   };

A residual model derives only from ``BaseEoS<N>``. This distinction determines
whether the type satisfies :cpp:concept:`IdealEoS` or
:cpp:concept:`ResidualEoS`.

Step 2: state the physical equations
------------------------------------

Write the equations and units before writing the kernels. For NASA-7, the
standard-state functions for species ``i`` are

.. math::

   \frac{c_{p,i}^\circ}{R}
   =a_{0,i}+a_{1,i}T+a_{2,i}T^2+a_{3,i}T^3+a_{4,i}T^4,

.. math::

   \frac{h_i^\circ}{RT}
   =a_{0,i}+\frac{a_{1,i}}{2}T+\frac{a_{2,i}}{3}T^2
    +\frac{a_{3,i}}{4}T^3+\frac{a_{4,i}}{5}T^4
    +\frac{a_{5,i}}{T},

.. math::

   \frac{s_i^\circ}{R}
   =a_{0,i}\ln T+a_{1,i}T+\frac{a_{2,i}}{2}T^2
    +\frac{a_{3,i}}{3}T^3+\frac{a_{4,i}}{4}T^4+a_{6,i}.

The ideal-mixture Helmholtz equation is

.. math::

   a^\mathrm{ideal}
   =\sum_i x_i\left[
     h_i^\circ-Ts_i^\circ
     +RT\ln\!\left(\frac{x_i cRT}{p_{i,\mathrm{ref}}}\right)
     \right]-RT.

The argument of every logarithm is dimensionless. In particular,
:math:`x_i cRT/p_{i,\mathrm{ref}}` is the ratio of the species partial pressure
to its standard-state pressure.

Step 3: define input data with units
------------------------------------

Use a public ``SpeciesInput`` structure when the model takes one record per
species. Give every dimensional coefficient an explicit unit in the API
comment.

.. code-block:: cpp

   struct SpeciesInput {
       double a0;    // [-]
       double a1;    // [1/K]
       double a2;    // [1/K^2]
       double a3;    // [1/K^3]
       double a4;    // [1/K^4]
       double a5;    // [K], enthalpy integration coefficient
       double a6;    // [-], entropy integration coefficient
       double T_ref; // [K]
       double p_ref; // [Pa]
   };

NASA coefficient files normally contain multiple temperature ranges. A model
that stores one range per species should document that the caller must select a
range valid at ``T``. Do not silently extrapolate or switch ranges unless that
behavior is part of the model definition.

``T_ref`` and ``p_ref`` match the input interface of the bundled ``Nasa7``
model. The implementation forms
``c_ref = p_ref / (R T_ref)`` [mol/m\ :sup:`3`]. Its Helmholtz expression uses
the product ``c_ref T_ref = p_ref / R`` [mol K/m\ :sup:`3`], so ``T_ref``
cancels. It neither sets the valid temperature range nor changes the result
when ``p_ref`` is fixed.

Step 4: support fixed and runtime component counts
--------------------------------------------------

Use ``std::array`` for a compile-time component count and ``std::vector`` for
``std::dynamic_extent``:

.. code-block:: cpp

   private:
       using Storage = std::conditional_t<
           N == std::dynamic_extent,
           std::vector<SpeciesInput>,
           std::array<SpeciesInput, N>>;

       Storage data_{};

The constructors mirror the bundled models:

.. code-block:: cpp

   public:
       explicit ExampleNasa7(
           const std::array<SpeciesInput, N>& inputs)
           requires(N != std::dynamic_extent)
           : data_{inputs}
       {
       }

       explicit ExampleNasa7(std::span<const SpeciesInput> inputs)
           requires(N == std::dynamic_extent)
           : fugacity::BaseEoS<N>(inputs.size()),
             data_(inputs.begin(), inputs.end())
       {
       }

The fixed-size base stores no runtime count. The dynamic specialization stores
``inputs.size()`` and returns it from ``size()``.

Step 5: implement the standard-state functions
----------------------------------------------

Keep every calculation templated over ``Number``. Model parameters may remain
``double``, but expressions involving ``T``, ``c``, ``x``, or ``rho_i`` must
preserve the active number type used by automatic differentiation.

The following helpers implement the NASA-7 equations without a derivation:

.. code-block:: cpp

   private:
       template<std::floating_point Number>
       [[nodiscard]] FUGACITY_ALWAYS_INLINE Number
       standard_enthalpy(std::size_t i, Number T) const
       {
           const auto& a = data_[i];
           const Number T2 = T * T;
           const Number T3 = T2 * T;
           const Number T4 = T3 * T;
           const Number h_over_RT =
               a.a0 + (a.a1 * T / Number{2})
               + (a.a2 * T2 / Number{3})
               + (a.a3 * T3 / Number{4})
               + (a.a4 * T4 / Number{5}) + (a.a5 / T);
           return fugacity::ideal_gas_constant<Number> * T * h_over_RT;
       }

       template<std::floating_point Number>
       [[nodiscard]] FUGACITY_ALWAYS_INLINE Number
       standard_entropy(std::size_t i, Number T) const
       {
           const auto& a = data_[i];
           const Number T2 = T * T;
           const Number T3 = T2 * T;
           const Number T4 = T3 * T;
           const Number s_over_R =
               (a.a0 * std::log(T)) + (a.a1 * T)
               + (a.a2 * T2 / Number{2})
               + (a.a3 * T3 / Number{3})
               + (a.a4 * T4 / Number{4}) + a.a6;
           return fugacity::ideal_gas_constant<Number> * s_over_R;
       }

Both functions return molar quantities: ``standard_enthalpy`` returns [J/mol]
and ``standard_entropy`` returns [J/(mol K)]. ``FUGACITY_ALWAYS_INLINE`` keeps
the helper bodies visible to Enzyme at the optimization levels used by the test
presets.

Step 6: implement molar Helmholtz energy
----------------------------------------

Expand the logarithm so that a zero mole fraction is handled by
:cpp:func:`fugacity::xlnx`:

.. code-block:: cpp

   public:
       template<std::floating_point Number>
       [[nodiscard]] Number calc_helmholtz(
           Number c, const Number* x, Number T) const
       {
           const Number R = fugacity::ideal_gas_constant<Number>;
           const Number RT = R * T; // [J/mol]
           Number a = -RT;          // [J/mol]

           for (std::size_t i = 0; i < this->size(); ++i) {
               const Number h0 = standard_enthalpy(i, T); // [J/mol]
               const Number s0 = standard_entropy(i, T);  // [J/(mol K)]
               const Number pressure_ratio =
                   c * RT / data_[i].p_ref; // [-]

               a += x[i] *
                    (h0 - (T * s0) + (RT * std::log(pressure_ratio)));
               a += RT * fugacity::xlnx(x[i]);
           }
           return a;
       }

Preconditions are ``c > 0`` [mol/m^3], ``T > 0`` [K], and a mole-fraction
array of length ``size()``. Physical callers should supply
:math:`x_i\ge0` and :math:`\sum_i x_i=1`.

Step 7: implement the density kernels
-------------------------------------

For the NASA-7 ideal model, a convenient per-component density equation is

.. math::

   \Psi_i
   =\rho_i\left(h_i^\circ-Ts_i^\circ-RT\right)
    +p_{i,\mathrm{ref}}
       \left(y_i\ln y_i\right),
   \qquad
   y_i=\frac{\rho_iRT}{p_{i,\mathrm{ref}}}.

Here :math:`\rho_i` is in mol/m^3, :math:`y_i` is dimensionless, and
:math:`\Psi_i` is in J/m^3. The pressure multiplying :math:`y_i\ln y_i` has
units Pa, equivalent to J/m^3.

Implement the equation once and reuse it in both public density functions:

.. code-block:: cpp

   private:
       template<std::floating_point Number>
       [[nodiscard]] FUGACITY_ALWAYS_INLINE Number density_component(
           std::size_t i, Number rho, Number T) const
       {
           const Number R = fugacity::ideal_gas_constant<Number>;
           const Number RT = R * T; // [J/mol]
           const Number h0 = standard_enthalpy(i, T); // [J/mol]
           const Number s0 = standard_entropy(i, T);  // [J/(mol K)]
           const Number y = rho * RT / data_[i].p_ref; // [-]

           return rho * (h0 - (T * s0) - RT)
                + data_[i].p_ref * fugacity::xlnx(y); // [J/m^3]
       }

   public:
       template<std::floating_point Number>
       [[nodiscard]] Number calc_helmholtz_density(
           const Number* rho_i, Number T) const
       {
           Number psi{0}; // [J/m^3]
           for (std::size_t i = 0; i < this->size(); ++i) {
               psi += density_component(i, rho_i[i], T);
           }
           return psi;
       }

       template<std::floating_point Number>
       void calc_partial_helmholtz(
           const Number* rho_i, Number T, Number* out) const
       {
           for (std::size_t i = 0; i < this->size(); ++i) {
               out[i] = density_component(i, rho_i[i], T); // [J/m^3]
           }
       }

For a residual model, a unique physical per-component decomposition may not
exist. Any documented decomposition is acceptable if its sum equals the scalar
``calc_helmholtz_density`` result. The bundled cubic models use
:math:`\Psi_i=(\rho_i/c)\Psi`.

Step 8: check the concepts and expose the header
------------------------------------------------

Add compile-time checks for both component-count forms:

.. code-block:: cpp

   static_assert(fugacity::IdealEoS<ExampleNasa7<2>>);
   static_assert(
       fugacity::IdealEoS<ExampleNasa7<std::dynamic_extent>>);

For a residual model, use ``ResidualEoS`` instead. Include the new public header
from ``include/fugacity/fugacity.hpp`` so users of the umbrella header can reach
the model.

Step 9: write model-specific unit tests
---------------------------------------

A new model needs tests of its physical equations in addition to generic
interface tests. Generic consistency tests cannot detect an error copied into
all three Helmholtz kernels.

Test the defining equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use an independent implementation of the NASA-7 relations. The following
example checks a one-component model at its reference state:

.. code-block:: cpp

   #include "fugacity/core/core_calculations.hpp"
   #include "fugacity/core/eos_pair.hpp"
   #include "fugacity/core/numbers.hpp"
   #include "fugacity/residual_models/no_residual.hpp"
   #include "example_nasa7.hpp"

   #include <array>
   #include <boost/ut.hpp>
   #include <cmath>

   using Input = ExampleNasa7<1>::SpeciesInput;

   constexpr Input n2{
       .a0 = 3.53100528,
       .a1 = -1.23660988e-4,
       .a2 = -5.02999433e-7,
       .a3 = 2.43530612e-9,
       .a4 = -1.40881235e-12,
       .a5 = -1046.97628,
       .a6 = 2.96747038,
       .T_ref = 298.15, // K
       .p_ref = 1.0e5,  // Pa
   };

   double nasa7_cp(const Input& a, double T) // [J/(mol K)]
   {
       const double R = fugacity::ideal_gas_constant<double>;
       return R * (a.a0 + a.a1*T + a.a2*T*T
                   + a.a3*T*T*T + a.a4*T*T*T*T);
   }

   double nasa7_h(const Input& a, double T) // [J/mol]
   {
       const double R = fugacity::ideal_gas_constant<double>;
       return R*T * (a.a0 + a.a1*T/2.0 + a.a2*T*T/3.0
                     + a.a3*T*T*T/4.0 + a.a4*T*T*T*T/5.0
                     + a.a5/T);
   }

   double nasa7_s(const Input& a, double T) // [J/(mol K)]
   {
       const double R = fugacity::ideal_gas_constant<double>;
       return R * (a.a0*std::log(T) + a.a1*T + a.a2*T*T/2.0
                   + a.a3*T*T*T/3.0 + a.a4*T*T*T*T/4.0
                   + a.a6);
   }

   int main()
   {
       using namespace boost::ut;

       "NASA-7 reference properties"_test = [] {
           const std::array<Input, 1> inputs{{n2}};
           const fugacity::EoS eos{
               ExampleNasa7<1>{inputs}, fugacity::NoResidual<1>{}};

           const double R = fugacity::ideal_gas_constant<double>;
           const double T = n2.T_ref;               // K
           const double c = n2.p_ref / (R * T);     // mol/m^3
           const std::array<double, 1> x{1.0};       // [-]

           expect(std::abs(fugacity::calc_cp(eos, c, x, T)
                           - nasa7_cp(n2, T)) < 1.0e-9);
           expect(std::abs(fugacity::calc_enthalpy(eos, c, x, T)
                           - nasa7_h(n2, T)) < 1.0e-7);
           expect(std::abs(fugacity::calc_entropy(eos, c, x, T)
                           - nasa7_s(n2, T)) < 1.0e-9);
       };

       return cfg<>.run();
   }

Use tolerances justified by the expected magnitude and floating-point path.
Check more than one temperature within the coefficient range; a single
reference temperature can miss errors in higher-order coefficients.

Test molar-density consistency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a representative mixture, verify both interface identities:

.. code-block:: cpp

   const double T = 350.0; // K
   const std::array<double, 2> rho_i{30.0, 70.0}; // mol/m^3
   const double c = rho_i[0] + rho_i[1];          // mol/m^3
   const std::array<double, 2> x{
       rho_i[0] / c, rho_i[1] / c};               // [-]

   const double psi = model.calc_helmholtz_density(rho_i.data(), T);
   const double a = model.calc_helmholtz(c, x.data(), T);

   std::array<double, 2> psi_i{};
   model.calc_partial_helmholtz(rho_i.data(), T, psi_i.data());

   expect(std::abs(psi - c*a) < tolerance); // [J/m^3]
   expect(std::abs(psi - (psi_i[0] + psi_i[1])) < tolerance);

Repeat this test for fixed and dynamic component counts, multiple compositions,
and each intended temperature regime.

Use the repository contract tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The test support in ``tests/support/eos_test_suite.hpp`` checks:

* the ``EquationOfState`` density identities;
* public scalar-property identities;
* container and pointer overloads;
* temperature and size preconditions;
* every public first derivative against adaptive finite differences;
* deterministic sampling over a declared valid domain; and
* ideal-gas or residual dilute-limit identities.

Register those checks from the model's test file:

.. code-block:: cpp

   const auto eos = fugacity::EoS{
       ExampleNasa7<2>{binary_inputs},
       fugacity::NoResidual<2>{}};

   const fugacity_test::eos_test_fixture fixture{
       .contribution = eos.ideal(),
       .eos = eos,
       .states = {
           {.c = 40.0, .x = {0.5, 0.5}, .T = 280.0,
            .effective_molar_mass = 0.029, .label = "cool"},
           {.c = 200.0, .x = {0.7, 0.3}, .T = 500.0,
            .effective_molar_mass = 0.029, .label = "warm"},
       },
       .domain = {
           .c_min = 20.0,  // mol/m^3
           .c_max = 250.0, // mol/m^3
           .T_min = 250.0, // K
           .T_max = 600.0, // K
           .minimum_mole_fraction = 0.01,
           .seed = 0xC0FFEE,
           .random_samples = 50,
       },
   };

   fugacity_test::register_eos_contract_tests(fixture);
   fugacity_test::register_ideal_gas_contract_tests(fixture);

For a residual model, call ``register_residual_contract_tests`` with a dilute
concentration that lies inside the declared domain. Also compare fixed-size and
runtime-size implementations with
``register_static_dynamic_equivalence_tests``.

Step 10: add and run the test
-----------------------------

Place the test in ``tests/test_<model_name>.cpp``. The test CMake configuration
discovers ``.cpp`` files automatically. Reconfigure after adding a file, then
build and run the test suite:

.. code-block:: console

   $ cmake --preset debug \
       -D Enzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
   $ cmake --build --preset debug
   $ ctest --preset debug

Run an optimized preset as well. Enzyme is part of the compiler pipeline, so a
test that passes only at one optimization level is not sufficient.

Automatic-differentiation requirements
---------------------------------------

Follow these rules for functions that contribute to the differentiated
Helmholtz expression:

* Define calculation templates in a public header so Clang can instantiate and
  differentiate them.
* Preserve ``Number`` for active arithmetic. Do not convert ``c``, ``x[i]``,
  ``rho_i[i]``, or ``T`` to ``double``.
* Annotate small helper functions with ``FUGACITY_ALWAYS_INLINE``.
* Prefer standard arithmetic and math functions that Enzyme supports. Avoid
  opaque library calls inside a differentiated kernel.
* Keep the scalar ``calc_helmholtz_density`` path direct. Reverse-mode
  differentiation uses it to obtain :math:`\partial\Psi/\partial\rho_i` [J/mol].
* Test derivatives over the full intended state domain, including dilute,
  dense, low-temperature, and high-temperature states that remain physically
  valid for the model.

Implementation checklist
------------------------

Before submitting a model, verify that:

#. the API comment states the model name and defining Helmholtz equations;
#. every input, output, and coefficient has physical units;
#. the class supports the intended fixed or dynamic component counts;
#. all three Helmholtz calculations satisfy the density identities;
#. the ideal or residual concept ``static_assert`` passes;
#. the umbrella header includes the model;
#. tests compare against independent physical equations;
#. automatic derivatives agree with numerical derivatives over a valid domain;
#. fixed and dynamic implementations agree, when both are provided; and
#. the documentation build completes without warnings.
