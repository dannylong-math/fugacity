Fugacity
========

Fugacity is a header-only C++23 library for calculating thermodynamic properties
from Helmholtz-energy equations of state. Supply an ideal contribution and a
residual contribution; Fugacity uses automatic differentiation to calculate
pressure, caloric properties, chemical potentials, fugacities, and property
derivatives.

All public thermodynamic inputs and outputs use SI units. Molar concentration,
not mass density, is the independent density variable.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   tutorial
   implementing_a_new_eos
   api/index
