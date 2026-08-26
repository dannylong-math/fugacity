# EOS kernel performance V&V report

## Context of use and claims

The primary context is optimized Release/O3 scalar CPU evaluation using the repository's
Clang/Enzyme toolchain. The claimed change is an algebraically equivalent scheduling of the
dynamic-double Peng–Robinson molar attractive mixing term for 10 or more components. It does
not change the equation, coefficients, mixing rule, units, valid-state domain, or
density/reverse/fugacity expression path.

Code and integration verification apply. Calculation/solution convergence does not apply
because the library evaluates pointwise closed-form expressions rather than a discretized or
iterative method. Physical model validation and uncertainty quantification do not apply
because the scientific model and coefficient data are unchanged.

## Independent protocol

The V&V audit evaluated exact integration commit
`3f5fae4ca6da9cb5c633e00d2dcbc36a6c157a4e` using:

- 50-digit independent Peng–Robinson expressions;
- refined finite differences for complete public Release derivatives;
- exact base-revision versus candidate comparisons;
- molar/density/partial consistency and permutation checks; and
- valid states spanning thresholds, block tails, dilute/high-density regimes, and
  temperature-dependent alpha branches away from their nondifferentiable roots.

The installed package identifies the plugin as Enzyme 0.0.79 compatible with LLVM 22.1.7;
Clang 22.1.8 was used. Results apply to this exact installed plugin rather than every Enzyme
release.

## Release/O3 results

| Evidence | Comparisons | Result |
|---|---:|---|
| Multiprecision primal states | 208 | PASS |
| Multiprecision lambda values | 162 | PASS |
| Reconstructed pressure/cp/sound | 36 | PASS |
| Multiprecision fugacity gradients | 27 | PASS |
| Public dT finite differences | 264 | PASS |
| Public dc finite differences | 264 | PASS |
| Representative public dx finite differences | 1254 | PASS |

There were zero scientific-tolerance failures. The largest fractions of an acceptance bound
were `1.28e-5` for dT, `2.01e-4` for dc, and `3.46e-4` for dx. Exact Release base/candidate
comparisons retained bitwise molar Helmholtz, density, partial decomposition, pressure, and
all 4736 sampled fugacity outputs.

Raw ULP parity with the legacy Enzyme graph is not a global invariant. Examples include
134 ULP for internal lambda(2,0), 30 ULP for one sound-speed value in a broader translation
unit, and larger ULP counts for small higher public derivatives. The worst sampled relative
public-derivative difference was `7.94e-12`, and every such result independently passed its
scientific oracle. These are compiler/AD graph-rounding differences rather than equation
errors.

## Debug/O1 limitation

The existing Debug/ASan suite passes, with leak detection disabled because LSan aborts under
the ptrace-based harness. Expanded large-N nested composition-gradient checks are not
numerically reliable under O1: finite differences show order-one and sign-changing errors
for some cp, sound-speed, and related composition derivatives. An isolated N16 `calc_cp_dx`
case has normalized error `1.212e-3` under O1 versus `4.681e-11` under O3.

The user explicitly selected optimized Release as the primary context because Enzyme is
optimization-sensitive. Debug remains useful for diagnostics and ASan, but large-N nested
composition derivatives are not qualified scientific outputs in that profile.

## Verdict and limitations

PASS for the supported Release/O3 context. The optimization is not numerically qualified
for the expanded Debug/O1 nested-composition-AD routes. Results are specific to the installed
Clang/Enzyme toolchain and should be repeated after a material compiler or Enzyme upgrade.

Independent report artifact:
`/tmp/fugacity-eos-vv-3f5fae4/VV_REPORT.md`, SHA-256
`a21e9661383b6ca9c4120986becec820de148e8de5f7f3df5c0787d0338c8e33`.

