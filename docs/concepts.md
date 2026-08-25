---
title: Equation-of-state concepts
description: The model concepts and mutually consistent Helmholtz-energy interface.
---

# Equation-of-state concepts

The model concepts define the Helmholtz interface and distinguish ideal and
residual contributions. Calculation functions are templates over the
floating-point number type so Enzyme can differentiate the concrete model.

## `EquationOfState`

An equation-of-state contribution provides `size()` and three mutually
consistent Helmholtz calculations:

```cpp
template<std::floating_point Number>
Number calc_helmholtz(Number c, const Number* x, Number T) const;

template<std::floating_point Number>
Number calc_helmholtz_density(const Number* rho_i, Number T) const;

template<std::floating_point Number>
void calc_partial_helmholtz(
    const Number* rho_i, Number T, Number* out) const;
```

`calc_helmholtz` returns molar Helmholtz energy [J/mol]. The density functions
return total or per-component Helmholtz energy density [J/m^3]. They satisfy

$$
\Psi(\boldsymbol{\rho},T)=c\,a(c,\boldsymbol{x},T),\qquad
c=\sum_i\rho_i,\qquad x_i=\rho_i/c,
$$

and

$$
\Psi=\sum_i\Psi_i.
$$

See the generated [`EquationOfState` API entry](/api/) for the complete
compile-time requirements.

## `IdealEoS`

`IdealEoS` matches an `EquationOfState` that publicly derives from
`fugacity::BaseIdealEoS`.

## `ResidualEoS`

`ResidualEoS` matches an `EquationOfState` that does not derive from
`fugacity::BaseIdealEoS`.
