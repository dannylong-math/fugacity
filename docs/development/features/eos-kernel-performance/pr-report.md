# Dynamic Peng–Robinson kernel performance refinement

## Summary

This change accelerates the runtime-sized `double` Peng–Robinson molar Helmholtz kernel for
mixtures with at least 10 components. An eight-row, stack-only schedule reuses each
temperature/composition column factor across eight independent matrix rows. The public API,
equation, coefficient interpretation, density/fugacity graph, and dependency set are
unchanged.

The initial audit left ConstantCp, NASA7, NASA9, NoResidual, van der Waals, and static
Peng–Robinson unchanged because no clear evidence-backed restructuring opportunity was
found. The accepted dynamic Peng–Robinson change improves measured N50 pressure by 2.87x,
cp by 2.12x, and squared sound speed by 2.03x.

## Scientific and numerical rationale

- **Context of use:** optimized Release/O3 scalar CPU thermodynamic-property evaluation
  with the installed Clang/Enzyme toolchain.
- **Problem/model:** the dynamic generalized-cubic attractive mixing double sum repeatedly
  recalculated a factor that depends on the column species but not the matrix row.
- **Chosen method:** process eight rows together, compute one column factor per block, and
  update eight independent row accumulators.
- **Key assumptions:** valid finite thermodynamic states; Release/O3 is the production
  numerical profile; each row and final-row reduction retain their original order.
- **Evidence:** source profiling, paired production A/B timing, strict-FP assembly, exact
  base comparisons, independent 50-digit Peng–Robinson expressions, and refined finite
  differences.
- **Accuracy and limitations:** Release public values and derivatives pass all existing
  scientific tolerances. Raw legacy-graph ULP parity is not global because Enzyme is
  optimization-sensitive. Expanded large-N Debug/O1 nested composition derivatives are not
  scientifically qualified, as documented in `vv.md`.

No literature review or physical-data validation was required because neither the equation
of state nor its parameters changed.

## Architecture and compatibility

- **Selected design:** private K=8 blocked helper, selected only for
  `N == std::dynamic_extent`, `Number == double`, and runtime N >= 10.
- **Data flow:** `calc_helmholtz` retains its existing scalar b-mixing accumulation, then
  evaluates the attractive quadratic in row blocks. Each row accumulates columns in the
  original order, and completed rows are added in the original order.
- **Unchanged paths:** dynamic N<10, static extents, `float`, `long double`, density,
  partial-density, chemical-potential, and fugacity calculations.
- **Rejected alternatives:** heap-allocated full cache, mutable model-owned cache,
  caller-provided workspace/API, float blocking, matrix-symmetry summation, and duplicated
  derivative caching.
- **Source API impact:** none.
- **ABI/distribution impact:** Fugacity is header-only and ships no standalone binary ABI.
  Object layout and public signatures are unchanged; consumers rebuild generated code.
- **Data/serialization impact:** none.
- **Decision records:** `plan.md`, `performance.md`, and `vv.md`.

## Implementation

- `include/fugacity/residual_models/cubic.hpp` adds one private static block size, one private
  helper, and a dynamic-double molar dispatch.
- `tests/test_peng_robinson.cpp` adds deterministic threshold/tail characterization,
  legacy-loop comparison, covariance/type canaries, differentiated-property checks, and a
  strict-Release N16 `calc_cp_dx` five-point finite-difference regression.
- `scripts/check_coverage.py` and `AGENTS.md` update the reviewed first-party coverage
  denominator from 1031/150/244 to 1081/151/266.
- No dependency, build-target, compiler-flag, fast-math, SIMD, GPU, or alternative-AD change
  was introduced.
- Target environment remains the repository-supported workstation CPU with Clang/Enzyme.

Why the new code is faster: the old runtime loop computed
`x[j] * abs(p[j] - q[j] * sqrt(T))` once for every `(i,j)` interaction. The new loop
computes it once for each column and eight-row block, reducing redundant arithmetic by
approximately a factor of eight inside full blocks. Clang emits eight independent scalar
FMA dependency chains, allowing instruction-level parallelism. It does not need packed SIMD,
reassociation, or relaxed floating-point arithmetic.

## Tests and adversarial audit

- **Test-first behavior:** the baseline formula was correct, so new characterization tests
  passed before implementation. The failing pre-change criterion was the measured
  performance deficit; no artificial functional failure was manufactured.
- **Boundary coverage:** N=7, 8, 9, 10, 15, 16, 17, 23, 24, 25, 31, 32, and 33, including
  threshold, complete-block, and partial-tail behavior.
- **Inputs:** deterministic nonuniform parameters/compositions and asymmetric nonzero
  interactions, seed `0xB10C5EED`.
- **Independent oracles:** long-double/multiprecision formulas, base-revision differential
  comparison, permutation covariance, molar/density identities, and five-point finite
  differences.
- **Release composition-gradient regression:** N16, all 16 directions, seed
  `0xF1D1FFB10C5EED11`, scaled error limit `2e-7`; passes.
- **Manual mutants:** tail omission, shifted column, and wrong second-block row factor were
  killed. An N10-to-N11 dispatch mutant survives numerical tests because it changes only
  performance policy; the benchmark evidence enforces the threshold.
- **Skeptic verdict:** READY for the primary Release context; no Release production defect.
- **Tooling gaps:** no configured property/shrinking framework, fuzz target/corpus, or
  mutation framework. Manual deterministic loops and mutants were used; nothing was
  installed.
- **Residual risk:** finite tests cannot exclude benchmark-specific special cases, though
  source inspection found none.

## Coverage

| Metric | Raw | Approved exclusions | Policy-adjusted | Requirement | Result |
|---|---:|---:|---:|---:|---|
| Lines | 1081/1081 | 0 | 100% | 100% | PASS |
| Functions | 151/151 | 0 | 100% | 100% | PASS |
| Branch outcomes | 266/266 | 0 | 100% | 100% | PASS |

### Coverage exclusions

None. No exclusion, suppression, or approval record is required.

## Verification and validation

- **Context and claim:** Release/O3 evaluation of the unchanged implemented equation and
  public derivatives within the repository's scientific tolerances.
- **Code verification:** 208 multiprecision primal states, 162 multiprecision lambda
  comparisons, 36 reconstructed pressure/cp/sound values, and 27 multiprecision fugacity
  gradients all pass.
- **Public derivative verification:** 264 dT, 264 dc, and 1254 representative dx finite-
  difference comparisons pass with zero failures. The largest fractions of an acceptance
  bound are `1.28e-5`, `2.01e-4`, and `3.46e-4`, respectively.
- **Exact preservation:** Release molar Helmholtz, density, partial decomposition, pressure,
  and 4736 sampled fugacity outputs are bitwise identical to the base.
- **Calculation/solution verification:** NOT APPLICABLE; this is a pointwise closed-form
  evaluation with no discretization or iterative solver.
- **Model validation:** NOT APPLICABLE; model and parameters are unchanged and no physical-
  accuracy claim is made.
- **Uncertainty:** multiprecision and stencil-refinement disagreement bound oracle error;
  no uncertain-input inference applies.
- **Strongest defensible conclusion:** PASS for Release/O3 on the installed Clang/Enzyme
  toolchain. Expanded Debug/O1 large-N nested composition gradients are not qualified.
- **Report:** `vv.md` and `/tmp/fugacity-eos-vv-3f5fae4/VV_REPORT.md` (SHA-256
  `a21e9661383b6ca9c4120986becec820de148e8de5f7f3df5c0787d0338c8e33`).

## Performance

- **Target:** scalar dynamic-double Peng–Robinson on Intel Core i7-1370P, Clang 22.1.8,
  installed Enzyme, `-O3 -march=native -DNDEBUG`, one pinned P-core.
- **Protocol:** real production classes, 15 alternating bracketed process pairs, at least
  0.1 s calibration per case, deterministic bootstrap interval with seed `0xC0FFEE`.
- **Environment limitation:** powersave governor, turbo/frequency scaling, and variable load;
  paired within-round speedups are the decision statistic. Hardware counters were blocked by
  `perf_event_paranoid=4`.

| Calculation | N10 speedup (95% interval) | N50 speedup (95% interval) |
|---|---:|---:|
| Raw molar Helmholtz | 1.50x [1.43, 1.57] | 2.91x [2.63, 2.96] |
| Pressure | 1.48x [1.34, 1.51] | 2.87x [2.74, 2.95] |
| cp | 1.41x [1.31, 1.44] | 2.12x [1.96, 2.18] |
| Sound speed squared | 1.38x [1.18, 1.46] | 2.03x [1.88, 2.12] |

Density and fugacity remain on the unchanged kernel and show no accepted speedup claim.
Static, small-N, and non-double canaries remain on the old path. The helper has no calls or
heap allocation, a 216-byte stack frame, and eight scalar FMA accumulator chains. Binary text
grew by 8856 bytes (1.55%); a sampled nested Enzyme derivative frame grew from 152 to 1080
bytes, but the public differentiated calculations remain materially faster.

Full results, raw hashes, and commands are in `performance.md` and
`/tmp/fugacity-eos-perf-ab-522c587/PERF2_REPORT.md` (SHA-256
`a8e5d8ec5847b6f84fdc3a54747e5c4bfa9722890cdfdb8b19dab89f2a3bc33b`).

## Quality gates

| Gate | Configuration/command | Result | Evidence |
|---|---|---|---|
| Unit tests | Debug, Release, release-max full CTest | PASS | 15/15 in each profile |
| Coverage | `coverage-check` | PASS | 1081/1081, 151/151, 266/266 |
| Coverage exclusions | Source/register audit | PASS | None |
| clang-tidy/static analysis | Full debug-tidy build/test | PASS | 15/15; zero first-party diagnostics |
| Sanitizers | Debug ASan/LSan outside ptrace sandbox | PASS | 15/15; leak detection enabled in final audit |
| ABI compatibility | Header-only/no binary ABI | NOT APPLICABLE | Public signatures/layout unchanged |
| Adversarial tooling | Deterministic properties/manual mutants | PASS with tooling gaps | Skeptic report |
| Integration/V&V | MP and complete Release derivative campaign | PASS | `vv.md` |
| Documentation | Doxygen + Sourcey | PASS | 29 pages, 34 HTML, 2119 links |
| Formatting/diff hygiene | clang-format + `git diff --check` | PASS | Exact final audit |
| Performance regression | Paired production A/B | PASS | Material N10/N50 speedups |

Final independent auditor verdict: **READY FOR HUMAN REVIEW**.

## Local Git state

- Base revision: `522c5879bf207be6c3679f27a7e080111dc20f70`
- Integration branch: `feature/eos-kernel-performance`
- Audited integration revision: `1689422d5be10b21d45c7f92b229af88d30e379a`
- PR-ready local branch: `feature/eos-kernel-performance-pr`
- Final squash commit: this report's containing squash commit; exact SHA is reported in the
  local handoff because a commit cannot contain its own hash.
- Squashed-tree equivalence check: required after the local squash and reported in the
  handoff.
- Remote status: not pushed; no PR opened.

## Reproduce

```sh
cmake --preset debug -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug --parallel
ASAN_OPTIONS=detect_leaks=0 ctest --preset debug --output-on-failure

cmake --preset release -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release --parallel
ctest --preset release --output-on-failure

cmake --preset release-max -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release-max --parallel
ctest --preset release-max --output-on-failure

cmake --preset debug-tidy -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug-tidy --parallel
ctest --preset debug-tidy --output-on-failure

cmake --preset coverage -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset coverage --parallel
cmake --build build/coverage --target coverage-check

cmake --preset docs
cmake --build --preset docs

find include tests benchmarks -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
  -exec clang-format --dry-run --Werror {} +
git diff --check 522c587..HEAD
```

The throwaway benchmark and independent V&V harnesses are intentionally not committed.
Their exact compile/run commands and artifact hashes are recorded in `performance.md`,
`vv.md`, and the `/tmp` reports cited above.

## Documentation

- Scope, architecture, decisions, and task state: `plan.md`
- Baseline, A/B timings, causal code generation, and raw hashes: `performance.md`
- Numerical context, independent oracles, and limitations: `vv.md`
- Test, coverage, audit, and gate evidence: `verification.md`

## Risks, limitations, tooling gaps, and follow-up

- Expanded large-N Debug/O1 nested composition gradients are not qualified; use Release/O3
  for scientific results.
- Repeat Release V&V and performance evidence after a material Clang or Enzyme upgrade.
- The installed Enzyme package reports 0.0.79/LLVM 22.1.7 compatibility, but package metadata
  may not fully identify the plugin source revision.
- `/tmp` benchmark/audit/V&V artifacts are ephemeral; archive them separately if long-term
  reproduction beyond committed summaries is required.
- Property/shrinking, fuzz, and mutation frameworks are not configured. A future tooling
  project should evaluate a Clang-22-compatible Mull setup and a deterministic property or
  libFuzzer workflow before installation.
- Results are specific to the sampled i7-1370P/Clang/Enzyme environment; other CPUs and
  toolchains require their own A/B evidence.

## Reviewer guide

- [ ] Confirm the specialization guard and old-path fallbacks in
  `include/fugacity/residual_models/cubic.hpp`.
- [ ] Check that row, column, and final accumulation order are preserved.
- [ ] Review threshold/tail fixtures and the strict-Release N16 composition-gradient oracle
  in `tests/test_peng_robinson.cpp`.
- [ ] Review the Release/O3 context and Debug/O1 limitation in `vv.md`.
- [ ] Confirm all coverage metrics are 100% with no exclusions in `verification.md`.
- [ ] Review paired speedups, noise controls, stack/code-size costs, and unchanged density/
  fugacity scope in `performance.md`.
- [ ] Inspect the final squashed diff against base
  `522c5879bf207be6c3679f27a7e080111dc20f70`.
