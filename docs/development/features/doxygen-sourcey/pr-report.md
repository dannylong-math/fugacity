# Migrate documentation and close repository quality gates

## Summary

This feature replaces Sphinx/Sphinx-Immaterial with Doxygen XML and Sourcey
3.6.5, renders mathematics locally with KaTeX 0.18.4, and publishes the checked
static site from `docs/dist`. Its expanded scope also removes the inherited
clang-tidy and clang-format diagnostic baselines, enables tidy warnings as
errors, introduces a correct multi-executable coverage gate, repairs a
multiprecision test-oracle lifetime defect, and makes the ideal-model benchmark
registry consistent with its default sweep.

The independent technical audit passed every configured gate. This report and
the authoritative records were then refreshed for the focused record/clean-
state re-audit that precedes the local squash.

## Scientific and numerical rationale

- Context of use and decision supported: maintain accurate user/API
  documentation and a trustworthy zero-warning/complete-coverage development
  baseline for a header-only Helmholtz-energy equation-of-state library.
- Problem/model: documentation and software-quality infrastructure; no
  thermodynamic model or production numerical method changed.
- Chosen method: Doxygen XML → filtered public XML → Sourcey → parser-aware
  KaTeX → exact static-artifact checks. Static-analysis findings are either
  fixed or covered by the user-approved risk-tiered policy in [plan.md](plan.md).
- Key assumptions: public documentation excludes `fugacity::detail` and Enzyme
  implementation declarations; normal consumers do not require docs tools;
  source API, numerical expression order, layouts, Enzyme spelling, tolerances,
  and deterministic seeds remain protected.
- Literature/evidence: no scientific literature, coefficient data, or
  validation evidence changed.
- Expected accuracy, uncertainty, and limitations: numerical accuracy and
  physical-model uncertainty are unchanged. Existing scientific tests and the
  compatibility probe are regression evidence, not new model validation.

## Architecture and compatibility

- Selected design: CMake owns the docs, tidy, benchmark-registry, and coverage
  entry points. Coverage runs all 15 test executables independently and unions
  source-site outcomes rather than merging incompatible header-only profiles.
- Important interfaces and data flow: `docs` preset/target remains stable;
  `coverage-check` is the 100% completion gate; `debug-tidy` treats enabled
  diagnostics as errors; both benchmark translation units are explicitly
  analyzed.
- Alternatives rejected: shell-only docs orchestration; blanket warning
  suppression; bounds/IR rewrites solely to satisfy generic checks; the native
  15-binary merged branch report, which emits 270 mismatched-function warnings.
- Source API impact: none. Public declarations are unchanged.
- ABI/distribution mode and impact: no separately compiled ABI exists. A
  cleanup-base/current consumer nevertheless produced identical sizes and
  alignments.
- Data/serialization impact: none.
- Decision records: [plan.md](plan.md) and [verification.md](verification.md).

## Implementation

- Major components: Doxygen/Sourcey/KaTeX pipeline; strict math, XML, API,
  redirect, and link checks; warning-clean tidy configuration; exact qualified
  suppressions; authoritative formatting; per-executable coverage checker;
  focused precondition tests; owning multiprecision helpers; canonical ideal
  benchmark registry and registry oracle.
- Generated API/output inventory: 50 public core overloads,
  `ideal_gas_constant`, three concepts, four legacy URLs, and no
  `fugacity::detail` or Enzyme API leakage; 2,119 exact local targets and
  fragments pass.
- Static-analysis disposition: 426 original findings = 293 covered by the two
  approved global disables
  (`cppcoreguidelines-pro-bounds-avoid-unchecked-container-access` and
  `portability-avoid-pragma-once`) + 62 exact qualified suppressions + 71
  actionable fixes.
- Numerical behavior: production calculations and floating-point expression
  order are unchanged. Six hexadecimal thermodynamic outputs matched the
  cleanup base exactly.
- Dependencies/build changes: exact `sourcey@3.6.5` and `katex@0.18.4` in npm
  lockfile version 3, with Node.js 22.12 or newer; `npm audit` reports zero
  known vulnerabilities. No cleanup dependency was added. Coverage and
  registry checkers use Python's standard library.
- Target environments: Linux workstation CPU with Clang/LLVM/clang-tidy
  22.1.8, Enzyme, CMake 3.28.3, Doxygen 1.9.8, Node.js 24.19.0, npm 11.17.0,
  and Python 3.12.3.

## Tests and adversarial audit

- Unit tests and observed test-first failures: initial documentation checks
  exposed 879 broken links, incomplete API output, and math-renderer edge
  cases. Tidy/format started at 426 unique and 89 diagnostics. Coverage exposed
  15 reachable precondition branches. The registry oracle caught 28 unmatched
  default benchmark filters.
- Independent oracles: exact deployed-artifact checker; owning-return static
  assertions; tidy scope/error canaries; external consumer layout and
  hexadecimal-output comparison; coverage uncovered-branch mutant; 336-name
  benchmark registry enumeration.
- Property/metamorphic/differential/fuzz/mutation evidence: existing
  deterministic scientific tests and a focused coverage negative-control
  mutant were used. No repository mutation, fuzz, or shrinking property
  framework is configured.
- Skeptic findings resolved: URL layout/API completeness, math parsing and
  token collisions, invalid tidy configuration, over-broad analysis risks,
  multiprecision expression lifetime, coverage aggregation, and benchmark
  registry/default mismatch.
- Tooling gaps: no approved mutation/fuzz/property framework, workflow linter,
  documentation-example compiler, or complete external-link reachability gate.
- Residual test risk: coverage cannot discover future never-instantiated header
  templates; the function identity assumes one source function declaration per
  line. Browser QA sampled representative pages rather than every viewport.

## Coverage

| Metric | Raw | Approved exclusions | Policy-adjusted | Requirement | Result |
|---|---:|---:|---:|---:|---|
| Lines | 1031/1031 (100%) | 0 | 1031/1031 (100%) | 100% | PASS |
| Functions | 150/150 (100%) | 0 | 150/150 (100%) | 100% | PASS |
| Branch outcomes | 244/244 (100%) | 0 | 244/244 (100%) | 100% | PASS |

### Coverage exclusions

No exclusion, suppression, or unreachable-code approval exists. Raw and
policy-adjusted metrics are identical. A temporary uncovered branch made the
gate fail at 1032/1034 lines and 245/246 branch outcomes; the mutant was removed.

The native merged report remains informational: 1013/1013 lines, 150/150
functions, and 242/244 branch outcomes with 270 mismatched-function warnings.

## Verification and validation

- Context of use and claims: software/documentation correctness and regression
  preservation; no new physical claim.
- Integration/canonical problems: Debug ASan, debug-tidy, Release, and
  release-max each pass 15/15 existing tests.
- Code verification: full existing analytic, derivative, teqp, and contract
  tests pass; tidy scope and coverage gates have independent negative controls.
- Calculation/solution verification: not applicable; no solver/discretization
  changed.
- Observed convergence or invariant evidence: existing derivative-oracle and
  thermodynamic identity tests are unchanged regression evidence.
- Experimental/model validation: not applicable; no experimental dataset or
  physical-model claim changed.
- Uncertainty treatment: not applicable.
- Reproduction report: [verification.md](verification.md).
- Limitations and strongest defensible conclusion: documentation and quality
  infrastructure meet their declared gates without changing the public source
  API or observed production numerical results. This is not new validation of
  the equations of state.

## Performance

- Scope and target: performance optimization is NOT APPLICABLE.
- Benchmark/profile environment: registry and timed smoke only; no controlled
  benchmark baseline was collected.
- Baseline and result: NOT CONFIGURED; no stored baseline or threshold.
- Numerical-equivalence evidence: benchmark calculation/input expressions are
  unchanged; production compatibility probe outputs are bitwise identical.
- API/ABI impact: none.
- Results by target environment: ideal static/dynamic and cubic comparison
  smoke filters executed successfully on the local CPU.
- Remaining bottlenecks: not assessed.

## Quality gates

| Gate | Configuration/command | Result | Evidence |
|---|---|---|---|
| Unit tests | Debug, debug-tidy, Release, release-max | PASS: 15/15 each | [verification.md](verification.md) |
| Line/function/branch coverage | `coverage-check` | PASS: 1031/1031, 150/150, 244/244 | [verification.md](verification.md) |
| Coverage-exclusion audit | feature tree | PASS: zero exclusions | Coverage section above |
| clang-tidy/static analysis | 15 test + 2 benchmark TUs, warnings as errors | PASS: zero diagnostics | [verification.md](verification.md) |
| Sanitizers | Debug ASan/LSan | PASS: 15/15 | [verification.md](verification.md) |
| ABI compatibility | header-only; layout probe | NOT APPLICABLE as ABI gate; layout PASS | [verification.md](verification.md) |
| Adversarial tooling | committed tests, canaries, coverage mutant, skeptic | PASS | [verification.md](verification.md) |
| Integration/V&V | existing numerical regressions and compatibility probe | PASS as regression; new physical V&V NOT APPLICABLE | [verification.md](verification.md) |
| Documentation | clean docs build and exact output checker | PASS: 29 pages, 34 HTML, 2,119 links | [verification.md](verification.md) |
| Formatting | 41-file dry run and idempotence | PASS: zero diagnostics | [verification.md](verification.md) |
| Performance regression | no optimization or stored baseline | NOT APPLICABLE / NOT CONFIGURED | Performance section above |

## Local Git state

- Base revision: `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`.
- Integration branch: `sourcey`.
- Independently audited technical revision:
  `345f94a2cf041672552eb9a19b18158fde9e9f3e`.
- Record-only corrections: current integration-branch HEAD, subject to focused
  record/clean-state re-audit.
- PR-ready local branch: created only after the focused re-audit.
- Final squash commit: created only after the focused re-audit; exact SHA is
  recorded in the final handoff because a commit cannot contain its own SHA.
- Squashed-tree equivalence check: required after squash.
- Remote status: not pushed; no PR opened.

## Reproduce

```sh
cmake --preset docs
cmake --build --preset docs

clang-tidy-22 --verify-config -config-file=.clang-tidy
cmake --preset debug-tidy -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug-tidy --parallel
ctest --preset debug-tidy --output-on-failure

cmake --preset debug -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug --parallel
ctest --preset debug --output-on-failure

cmake --preset release -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release --parallel
ctest --preset release --output-on-failure

cmake --preset release-max -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release-max --parallel
ctest --preset release-max --output-on-failure

cmake --preset coverage -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset coverage --parallel
cmake --build build/coverage --target coverage-check

cmake -S . -B build/benchmark-tidy \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CLANG_TIDY=ON \
  -DBUILD_BENCHMARKS=ON \
  -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build build/benchmark-tidy --parallel
cmake --build build/benchmark-tidy --target bench-ideal-registry-check

find include tests benchmarks -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
  -exec clang-format --dry-run --Werror {} +

git diff --check 9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc..sourcey
```

## Documentation

- API/user docs: `docs/*.md`, `Doxyfile`, `docs/api-groups.dox`, and generated
  `docs/dist`.
- Architecture: [plan.md](plan.md).
- V&V/reproduction: [verification.md](verification.md).
- Benchmarks: registry contract in the benchmark source and checker; no
  performance report because no optimization claim was made.

## Risks, limitations, tooling gaps, and follow-up

- Enzyme emits noisy untagged `unknown tbaa` messages that reduce log signal but
  are not compiler or tidy diagnostics.
- `test_horner` and `test_xlnx` retain documented `std::random_device` usage.
- Coverage cannot detect a future never-instantiated template definition and
  assumes one source function declaration per line.
- No stored benchmark baseline/regression threshold exists.
- TSan/MSan, mutation, fuzzing, property shrinking, actionlint,
  documentation-example compilation, and complete external-link checking are
  not configured. UBSan remains prohibited with Enzyme.

## Reviewer guide

- [ ] Review Doxygen/Sourcey configuration and public XML filtering.
- [ ] Inspect strict artifact/math/API checks and compatibility redirects.
- [ ] Review the approved tidy policy, 426-finding disposition, and every
  check-qualified suppression.
- [ ] Inspect `scripts/check_coverage.py`, the 15 precondition tests, and the
  zero-exclusion coverage evidence.
- [ ] Review the multiprecision owning-return fix and lifetime assertions.
- [ ] Confirm benchmark registry/default consistency and that calculation/input
  expressions are unchanged.
- [ ] Review public declaration/layout/hexadecimal-output/Enzyme compatibility
  evidence in [verification.md](verification.md).
- [ ] Reproduce the principal quality gates with the commands above.
- [ ] Inspect the final squashed diff against
  `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`.
