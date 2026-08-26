# EOS kernel performance verification

## Integrated change

Integration branch `feature/eos-kernel-performance` is based on
`522c5879bf207be6c3679f27a7e080111dc20f70`. The production change is limited to
`BaseCubic::calc_helmholtz` in `include/fugacity/residual_models/cubic.hpp`:

- runtime-sized `double` models with N >= 10 use an allocation-free eight-row block;
- small dynamic models, static extents, `float`, and `long double` retain the old path; and
- density, partial-density, chemical-potential, and fugacity expression paths are unchanged.

No public API, coefficient, formula, data layout, dependency, CMake target, or relaxed
floating-point flag changed.

## Test-first and regression evidence

Because the baseline formula was correct, new characterization tests were first added and
shown to pass before production code changed. The only pre-change failure criterion was the
measured performance deficit.

The permanent Peng–Robinson suite now includes:

- direct molar, density, and partial Helmholtz values at N=7, 8, 9, 10, 15, 16, 17, 23,
  24, 25, 31, 32, and 33;
- deterministic nonuniform parameters/compositions with seed `0xB10C5EED`;
- threshold, complete-block, partial-tail, static, float, and long-double controls;
- all six lambda orders and representative pressure/cp/sound/fugacity routes;
- species-reversal covariance; and
- a strict-Release N16 `calc_cp_dx` five-point finite-difference regression over all 16
  composition directions, seed `0xF1D1FFB10C5EED11`, scaled tolerance `2e-7`.

The Release-only oracle is compiled when `NDEBUG` is set and neither `__FAST_MATH__` nor
`__NO_MATH_ERRNO__` is defined. Binary inspection confirmed it is present in Release and
absent in Debug and release-max.

## Implementer gate results

| Gate | Result |
|---|---|
| Debug/O1 + ASan full CTest | PASS, 15/15; `ASAN_OPTIONS=detect_leaks=0` |
| Release/O3 full CTest | PASS, 15/15 |
| release-max full CTest | PASS, 15/15 |
| Focused debug-tidy build/test | PASS, 1/1, zero first-party diagnostics |
| Strict Release N16 cp composition FD oracle | PASS, all 16 directions |
| Formatting dry-run | PASS |
| `git diff --check` | PASS |
| Coverage source-site union | PASS, 100% |

Leak detection is disabled because LSan aborts under the ptrace-based test harness; ASan
memory instrumentation otherwise ran. release-max is a compatibility profile with approved
reassociation, not the numerical baseline.

Coverage after the change is:

| Metric | Raw | Excluded | Policy-adjusted |
|---|---:|---:|---:|
| Lines | 1081/1081 | 0 | 100% |
| Source functions | 151/151 | 0 | 100% |
| Branch outcomes | 266/266 | 0 | 100% |

The reviewed denominator guard and AGENTS guidance were updated from 1031/150/244 to
1081/151/266. No exclusion was requested or approved.

## Independent test audit

The skeptic independently compared exact base and candidate builds over thresholds and
tails, exercised nonuniform/asymmetric states, checked comparator behavior, used an
independent long-double oracle, and killed three substantive manual mutants. A threshold
10-to-11 mutant survives numerically because it changes only when the optimization starts;
the paired performance matrix, not a correctness test, enforces that policy.

Release `calc_cp_dx` passed a separate five-point finite-difference campaign with maximum
normalized error `1.21e-10`, comparable to the stencil's own refinement uncertainty. The
skeptic recommended the N16 permanent regression subsequently added in `bef70e2`.

Audit artifact: `/tmp/fugacity-eos-perf-skeptic-3f5fae4/AUDIT.md`, SHA-256
`771e1db2d57507ad49aac250809b3f3fc02f0b606c8adaf8cf2dab1e5ce16d47`.

## Numerical V&V

Independent Release/O3 multiprecision and finite-difference evidence has zero scientific-
tolerance failures across 208 primal states, 162 lambda comparisons, 36 reconstructed
pressure/cp/sound values, 27 fugacity-gradient comparisons, and 1782 complete public
derivative comparisons. See `vv.md` for the context-of-use decision and exact limitations.

Debug/O1 expanded large-N nested composition derivatives are not numerically qualified:
the installed Enzyme can produce order-one/sign-changing differences there despite the
existing Debug suite passing. The user explicitly selected optimized Release as the primary
scientific context. This limitation must remain visible in reviewer/user documentation.

## Performance verification

The paired/bracketed 15-run production A/B demonstrates material improvements outside the
base-control noise:

- N10: pressure 1.48x, cp 1.41x, sound speed 1.38x;
- N50: pressure 2.87x, cp 2.12x, sound speed 2.03x; and
- raw molar Helmholtz: 1.50x at N10 and 2.91x at N50.

Strict-FP assembly shows no heap calls and eight independent scalar FMA accumulation chains.
The speedup comes from reusing each computed column factor across eight rows and exposing
instruction-level parallelism; it does not rely on packed SIMD, fast-math, or reassociation.

## Tooling gaps and residual risk

- No repository-configured property-testing, fuzzing, or mutation-testing workflow exists.
  Deterministic property loops and manual mutants were used; no dependency was installed.
- Hardware counters were unavailable because `perf_event_paranoid=4`; source, compiler
  remarks, assembly, and timings provide the causal evidence.
- Exact numerical/code-generation behavior is compiler and Enzyme sensitive. Repeat the
  Release V&V and performance matrix after a material toolchain change.
- The installed Enzyme package reports version 0.0.79/LLVM 22.1.7 compatibility, but package
  metadata may not fully identify the source revision. Evidence is tied to the plugin used.

## Final independent quality audit

The quality-gate auditor rebuilt exact integration commit `a976ed7`, found only two Markdown
EOF hygiene defects, and reported no production, numerical, performance, test, sanitizer,
coverage, static-analysis, API, or documentation-pipeline defect. After the documentation
owner removed those blank lines in `1689422`, the narrow re-audit confirmed:

- a clean worktree and exactly two documentation-line deletions from the previously audited
  tree;
- `git diff --check 522c587..HEAD` passes;
- the documentation pipeline passes with 29 pages, 34 HTML files, and 2119 checked links;
  and
- all prior unchanged gate evidence remains applicable.

Final verdict: **READY FOR HUMAN REVIEW**.
