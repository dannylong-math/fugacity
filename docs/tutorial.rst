Tutorial: assemble a model and calculate properties
===================================================

This tutorial constructs a two-component equation of state, evaluates its
thermodynamic properties, and calculates fugacities and derivatives. The
example uses a constant-heat-capacity ideal contribution and a Peng-Robinson
residual contribution.

What Fugacity expects
---------------------

Fugacity represents a complete equation of state as

.. math::

   a(c,\boldsymbol{x},T)
   =a^\mathrm{ideal}(c,\boldsymbol{x},T)
   +a^\mathrm{res}(c,\boldsymbol{x},T),

where ``a`` is molar Helmholtz energy [J/mol]. The class
:cpp:class:`fugacity::EoS` stores the two contributions. Property functions
differentiate their sum.

The independent state variables are:

.. list-table:: State variables
   :header-rows: 1
   :widths: 20 20 20 40

   * - C++ name
     - Symbol
     - Units
     - Meaning
   * - ``c``
     - :math:`c`
     - mol/m^3
     - Total molar concentration, :math:`n/V`.
   * - ``x[i]``
     - :math:`x_i`
     - -
     - Mole fraction of component ``i``.
   * - ``T``
     - :math:`T`
     - K
     - Absolute temperature.
   * - ``rho_i[i]``
     - :math:`\rho_i`
     - mol/m^3
     - Partial molar concentration, :math:`\rho_i=c x_i`.

Fugacity calculates pressure from ``c``, ``x``, and ``T``. It does not solve
for molar concentration from a specified pressure and temperature, and it does
not perform phase-equilibrium calculations.

Choose the contributions
------------------------

The bundled ideal models differ in their standard-state caloric equations:

* :cpp:class:`fugacity::ConstantCp` uses one constant molar isobaric heat
  capacity per species.
* :cpp:class:`fugacity::Nasa7` uses one NASA-7 coefficient range per species.
* :cpp:class:`fugacity::Nasa9` uses one NASA-9 coefficient range per species.

The bundled residual models are:

* :cpp:class:`fugacity::NoResidual`, for an ideal gas;
* :cpp:class:`fugacity::VanDerWaals`; and
* :cpp:class:`fugacity::PengRobinson`.

Each model's API entry gives its defining Helmholtz equation and a minimal
example. Select the two contributions independently, but keep their component
counts and component order identical.

Define the species data
-----------------------

This example uses nitrogen followed by carbon dioxide. Keep that order in the
ideal data, residual data, binary-interaction matrix, mole-fraction array,
partial-concentration array, and molar-mass calculation.

The constant-heat-capacity input for each species is:

* ``T_ref``: reference temperature [K];
* ``p_ref``: reference pressure [Pa];
* ``c_p``: molar isobaric heat capacity [J/(mol K)];
* ``h_ref``: molar enthalpy at the reference state [J/mol]; and
* ``s_ref``: molar entropy at the reference state [J/(mol K)].

.. code-block:: cpp

   constexpr std::size_t component_count = 2;
   using Ideal = fugacity::ConstantCp<component_count>;

   const std::array<Ideal::SpeciesInput, component_count> ideal_data{{
       // N2
       {.T_ref = 298.15,
        .p_ref = 1.0e5,
        .c_p = 29.12,
        .h_ref = 0.0,
        .s_ref = 191.61},
       // CO2
       {.T_ref = 298.15,
        .p_ref = 1.0e5,
        .c_p = 37.14,
        .h_ref = -393522.0,
        .s_ref = 213.79},
   }};

   const Ideal ideal{ideal_data};

These rounded values make the example concrete. Use a mutually consistent and
validated thermodynamic dataset for scientific calculations. In particular, a
constant ``c_p`` model is only appropriate over a temperature interval where
the constant approximation is acceptable.

The Peng-Robinson input for each species is:

* ``T_c``: critical temperature [K];
* ``P_c``: critical pressure [Pa]; and
* ``omega``: acentric factor [-].

.. code-block:: cpp

   using Residual = fugacity::PengRobinson<component_count>;

   const std::array<Residual::SpeciesInput, component_count> residual_data{{
       // N2
       {.T_c = 126.192, .P_c = 3.3958e6, .omega = 0.0372},
       // CO2
       {.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394},
   }};

   const std::array<double, component_count * component_count> kij{
       0.0, 0.0,
       0.0, 0.0,
   };

   const Residual residual{residual_data, kij};

``kij`` is a dimensionless, row-major matrix. Entry ``kij[i*N + j]`` is
:math:`k_{ij}`. This example uses zero interaction parameters; supply values
appropriate to the chosen parameter set when they are available.

Pair the models
---------------

Construct a complete equation of state by value:

.. code-block:: cpp

   const fugacity::EoS eos{ideal, residual};

Both contributions must report the same component count. A mismatch throws
``std::logic_error`` in a debug build.

Define a state
--------------

Set temperature, total molar concentration, and composition:

.. code-block:: cpp

   const double T = 320.0; // K
   const double c = 150.0; // mol/m^3
   const std::array<double, component_count> x{0.70, 0.30}; // mole fractions

   std::array<double, component_count> rho_i{}; // mol/m^3
   for (std::size_t i = 0; i < component_count; ++i) {
       rho_i[i] = c * x[i];
   }

The mole fractions must be nonnegative and sum to one. Fugacity does not
renormalize composition arrays. The scalar property functions accept
``std::array``, ``std::vector``, ``std::span``, or another contiguous container
with ``size()`` and ``data()``.

Calculate scalar properties
---------------------------

.. code-block:: cpp

   const double a = fugacity::calc_helmholtz(eos, c, x, T);       // J/mol
   const double p = fugacity::calc_pressure(eos, c, x, T);        // Pa
   const double u = fugacity::calc_internal_energy(eos, c, x, T); // J/mol
   const double h = fugacity::calc_enthalpy(eos, c, x, T);        // J/mol
   const double s = fugacity::calc_entropy(eos, c, x, T);         // J/(mol K)
   const double g = fugacity::calc_gibbs(eos, c, x, T);           // J/mol
   const double cv = fugacity::calc_cv(eos, c, x, T);             // J/(mol K)
   const double cp = fugacity::calc_cp(eos, c, x, T);             // J/(mol K)

The calculated properties use the enthalpy and entropy reference states carried
by the ideal model. Changing those references shifts absolute caloric properties
and chemical potentials but does not change pressure or residual fugacity
coefficients.

Calculate the speed of sound
----------------------------

The speed-of-sound function requires the mixture molar mass ``M`` [kg/mol].
For the composition above:

.. code-block:: cpp

   const std::array<double, component_count> molar_mass{
       0.0280134, // N2  [kg/mol]
       0.0440095, // CO2 [kg/mol]
   };

   double M = 0.0; // kg/mol
   for (std::size_t i = 0; i < component_count; ++i) {
       M += x[i] * molar_mass[i];
   }

   const double w_squared =
       fugacity::calc_sound_speed_squared(eos, c, x, T, M); // m^2/s^2
   const double w = std::sqrt(w_squared);                    // m/s

The function holds the supplied ``M`` constant when calculating its temperature,
concentration, and composition derivatives. If composition changes, calculate
the corresponding molar-mass contribution separately.

Calculate chemical potentials and fugacities
---------------------------------------------

The vector-valued functions use ``std::span``. Create fixed-extent spans from
the state arrays:

.. code-block:: cpp

   const std::span<const double, component_count> x_span{x};
   const std::span<const double, component_count> rho_span{rho_i};

   std::array<double, component_count> mu{};      // J/mol
   std::array<double, component_count> ln_phi{};  // dimensionless
   std::array<double, component_count> f{};       // Pa

   fugacity::calc_chemical_potential(
       eos, rho_span, T, std::span<double, component_count>{mu});

   fugacity::calc_log_fugacity_coeff(
       eos, c, x_span, T, rho_span,
       std::span<double, component_count>{ln_phi});

   fugacity::calc_fugacity(
       eos, rho_span, T, std::span<double, component_count>{f});

The functions evaluate

.. math::

   \mu_i=\left(\frac{\partial\Psi}{\partial\rho_i}\right)_{T,\rho_{j\ne i}}
   \quad [\mathrm{J/mol}],

.. math::

   \ln\varphi_i=\frac{\mu_i^\mathrm{res}}{RT}-\ln Z,
   \qquad Z=\frac{p}{cRT},

and

.. math::

   f_i=\rho_iRT\exp\!\left(\frac{\mu_i^\mathrm{res}}{RT}\right)
   \quad [\mathrm{Pa}].

When calling ``calc_log_fugacity_coeff``, supply mutually consistent arrays:
``rho_i[i] == c*x[i]``. The function does not reconstruct or validate one form
from the other.

Calculate derivatives
---------------------

Scalar properties have temperature, concentration, and composition derivative
functions. For pressure:

.. code-block:: cpp

   const double dp_dT = fugacity::calc_dp_dT(eos, c, x, T); // Pa/K
   const double dp_dc =
       fugacity::calc_dp_dc(eos, c, x, T); // Pa m^3/mol

   std::array<double, component_count> dp_dx{}; // Pa
   fugacity::calc_pressure_dx(eos, c, x, T, dp_dx);

``calc_dp_dT`` and ``calc_dp_dc`` expose the thermodynamic pressure partials
directly. ``calc_pressure_dT`` and ``calc_pressure_dc`` provide the same
partials using the naming pattern shared by the other scalar properties.

.. _composition-derivatives:

Composition derivatives
~~~~~~~~~~~~~~~~~~~~~~~

The ``*_dx`` functions treat :math:`x_0,\ldots,x_{N-1}` as independent
variables. A physical composition perturbation must also satisfy
:math:`\sum_i dx_i=0`.

For a binary mixture, take component 1 as the dependent component. The
directional derivative obtained by increasing :math:`x_0` while decreasing
:math:`x_1` is

.. math::

   \left.\frac{dF}{dx_0}\right|_{x_1=1-x_0}
   =\frac{\partial F}{\partial x_0}
    -\frac{\partial F}{\partial x_1}.

In code:

.. code-block:: cpp

   const double dp_dx0_on_simplex = dp_dx[0] - dp_dx[1]; // Pa

Apply the same subtraction to the gradient returned for any scalar property.

Use a runtime component count
-----------------------------

Omit the template count, store the inputs dynamically, and construct the models
from spans when the component count is not known at compile time:

.. code-block:: cpp

   using DynamicIdeal = fugacity::ConstantCp<>;
   using DynamicResidual = fugacity::PengRobinson<>;

   std::vector<DynamicIdeal::SpeciesInput> ideal_data = load_ideal_data();
   std::vector<DynamicResidual::SpeciesInput> residual_data =
       load_residual_data();
   std::vector<double> kij = load_binary_interactions();

   const DynamicIdeal ideal{
       std::span<const DynamicIdeal::SpeciesInput>{ideal_data}};
   const DynamicResidual residual{
       std::span<const DynamicResidual::SpeciesInput>{residual_data},
       std::span<const double>{kij}};
   const fugacity::EoS eos{ideal, residual};

   std::vector<double> x = load_mole_fractions();
   const double p = fugacity::calc_pressure(eos, c, x, T); // Pa

For :cpp:class:`fugacity::NoResidual` with a runtime component count, construct
``NoResidual<std::dynamic_extent>{component_count}``.

Complete example
----------------

The following program combines the fixed-size steps above:

.. code-block:: cpp

   #include <fugacity/fugacity.hpp>

   #include <array>
   #include <cmath>
   #include <cstddef>
   #include <iostream>
   #include <span>

   int main()
   {
       constexpr std::size_t N = 2;
       using Ideal = fugacity::ConstantCp<N>;
       using Residual = fugacity::PengRobinson<N>;

       const std::array<Ideal::SpeciesInput, N> ideal_data{{
           {.T_ref = 298.15, .p_ref = 1.0e5, .c_p = 29.12,
            .h_ref = 0.0, .s_ref = 191.61},
           {.T_ref = 298.15, .p_ref = 1.0e5, .c_p = 37.14,
            .h_ref = -393522.0, .s_ref = 213.79},
       }};
       const std::array<Residual::SpeciesInput, N> residual_data{{
           {.T_c = 126.192, .P_c = 3.3958e6, .omega = 0.0372},
           {.T_c = 304.1282, .P_c = 7.3773e6, .omega = 0.22394},
       }};
       const std::array<double, N * N> kij{0.0, 0.0,
                                           0.0, 0.0};

       const fugacity::EoS eos{
           Ideal{ideal_data}, Residual{residual_data, kij}};

       const double T = 320.0; // K
       const double c = 150.0; // mol/m^3
       const std::array<double, N> x{0.70, 0.30};
       const std::array<double, N> rho_i{c * x[0], c * x[1]}; // mol/m^3

       const double p = fugacity::calc_pressure(eos, c, x, T); // Pa
       const double h = fugacity::calc_enthalpy(eos, c, x, T); // J/mol
       const double s = fugacity::calc_entropy(eos, c, x, T);  // J/(mol K)
       const double cp = fugacity::calc_cp(eos, c, x, T);      // J/(mol K)

       std::array<double, N> ln_phi{};
       std::array<double, N> f{};
       const std::span<const double, N> x_span{x};
       const std::span<const double, N> rho_span{rho_i};
       fugacity::calc_log_fugacity_coeff(
           eos, c, x_span, T, rho_span, std::span<double, N>{ln_phi});
       fugacity::calc_fugacity(
           eos, rho_span, T, std::span<double, N>{f});

       std::cout << "p [Pa]       = " << p << '\n'
                 << "h [J/mol]    = " << h << '\n'
                 << "s [J/mol/K]  = " << s << '\n'
                 << "cp [J/mol/K] = " << cp << '\n';
       for (std::size_t i = 0; i < N; ++i) {
           std::cout << "component " << i
                     << ": ln(phi) = " << ln_phi[i]
                     << ", f [Pa] = " << f[i] << '\n';
       }
   }

Check the state before evaluation
---------------------------------

Before calculating properties, verify the following conditions in the calling
application:

* ``T > 0`` [K];
* ``c > 0`` [mol/m^3] for models containing logarithms of concentration;
* ``x.size() == eos.size()``;
* :math:`x_i\ge0` and :math:`\sum_i x_i=1`;
* ``rho_i[i] == c*x[i]`` [mol/m^3] when both representations are supplied;
* NASA coefficients are valid at ``T``;
* the species order is identical in every model and state array; and
* cubic models remain inside their logarithm domain, including
  :math:`b_m c<1`.

The library checks positive temperature in the public property functions and
checks selected size preconditions in debug builds. It does not validate every
model-specific physical domain condition.

