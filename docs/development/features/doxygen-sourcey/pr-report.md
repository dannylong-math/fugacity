# Migrate documentation from Sphinx to Doxygen and Sourcey

## Summary

This change replaces the Sphinx/Sphinx-Immaterial documentation pipeline with
Doxygen XML and Sourcey 3.6.5 while preserving the public `docs` CMake preset
and target. Narrative pages are Markdown, public C++ API pages come from
Doxygen XML, and KaTeX 0.18.4 renders equations locally. The generated Pages
artifact is `docs/dist`.

The migration-specific implementation and adversarial gates pass. The branch
is **not PR-ready under the current repository policy** because the unchanged
feature base has 426 unique first-party clang-tidy warning lines and 89
clang-format diagnostics. Resolving those gates requires production-source or
repository-policy work outside this documentation-only feature.

## Scientific and numerical rationale

- Context of use and decision supported: provide maintainable user and API
  documentation for a header-only Helmholtz-energy equation-of-state library.
- Problem/model: documentation tooling only; thermodynamic models and property
  calculations are unchanged.
- Chosen method: Doxygen emits public C++ XML, Sourcey builds the site, and a
  parser-aware postprocessor renders TeX with KaTeX.
- Key assumptions: public documentation excludes `fugacity::detail` and Enzyme
  implementation declarations; normal builds must not require documentation
  tools.
- Literature/evidence: no scientific literature or model data changed.
- Expected accuracy, uncertainty, and limitations: numerical accuracy and
  uncertainty are unchanged. The site checks exact local links and fragments,
  but complete external-link reachability and documentation-example compilation
  are not configured.

## Architecture and compatibility

- Selected design: CMake-owned Doxygen XML → filtered public XML → Sourcey →
  KaTeX → exact artifact/API checker.
- Important interfaces and data flow: `cmake --preset docs` and
  `cmake --build --preset docs` remain stable; Doxygen XML is generated under
  `build/doxygen/xml`; Pages publishes `docs/dist`.
- Alternatives rejected: a shell-only pipeline would bypass the repository's
  stable CMake entry point and the approved Rift-derived structure.
- Source API impact: none. Clang token comparison passed for all 17 headers
  whose documentation comments changed.
- ABI/distribution mode and impact: not applicable; Fugacity is header-only and
  has no separately compiled ABI. Source tokens are unchanged.
- Data/serialization impact: none.
- Decision records: [plan.md](plan.md) and [verification.md](verification.md).

## Implementation

- Major components: fatal-warning `Doxyfile`, Doxygen API groups and XML
  filtering, Sourcey configuration, Markdown guides, parser-aware KaTeX
  rendering, compatibility redirects, exact link/fragment/API checks, CMake and
  GitHub Pages integration.
- Numerical behavior: unchanged.
- Dependencies/build changes: exact direct pins `sourcey@3.6.5` and
  `katex@0.18.4`, locked with npm lockfile version 3; Node.js 22.12 or newer;
  Doxygen 1.9.8 or newer.
- Target environments: documentation and C++ regressions were verified on the
  local Linux/Clang/Enzyme environment recorded in [verification.md](verification.md).

## Tests and adversarial audit

- Unit tests and observed test-first failures: committed Node tests reject
  invalid/unmatched TeX and orphan Doxygen references. The first strict audit
  exposed 879 broken static targets and missing functions, constants, and
  concepts; the second exposed escaped-dollar and sentinel-collision cases.
- Independent oracles: the final skeptic re-audit passed with no actionable
  findings. The output checker validates 2,119 exact `href`/`src` targets and
  fragments, 50 core function overloads, `ideal_gas_constant`, three concepts,
  four legacy URLs, and internal-API exclusion.
- Property/metamorphic/differential/fuzz/mutation evidence: not configured and
  not required to distinguish this documentation-only change. Focused negative
  and collision fixtures cover the custom JavaScript transformations.
- Skeptic findings resolved: static URL layout, API completeness, legacy
  aliases, invalid edit links, attribute/code math protection, unmatched TeX,
  escaped-dollar TeX, placeholder collisions, orphan rejection, and strict
  deployed-artifact checking.
- Tooling gaps: no JavaScript coverage/mutation gate, documentation-example
  compiler, or complete external-link checker.
- Residual test risk: browser QA sampled the tutorial, Core API, and
  Peng-Robinson pages rather than every page at every viewport.

## Coverage

| Metric | Raw | Approved exclusions | Policy-adjusted | Requirement | Result |
|---|---:|---:|---:|---:|---|
| Lines | N/A | 0 | N/A | 100% | NOT APPLICABLE |
| Functions | N/A | 0 | N/A | 100% | NOT APPLICABLE |
| Branches | N/A | 0 | N/A | 100% | NOT APPLICABLE |

### Coverage exclusions

No coverage exclusion or suppression was added. Documentation is outside the
first-party denominator, and normalized production token streams are unchanged.

## Verification and validation

- Context of use and claims: verify a documentation migration, not a physical
  model or numerical-method claim.
- Integration/canonical problems: Debug, Release, and release-max each passed
  all 15 existing tests.
- Code verification: 17/17 changed headers are token-identical to the feature
  base; `git diff --check` passes.
- Calculation/solution verification: not applicable.
- Observed convergence or invariant evidence: not applicable.
- Experimental/model validation: not applicable.
- Uncertainty treatment: not applicable.
- Reproduction report: [verification.md](verification.md).
- Limitations and strongest defensible conclusion: the documentation tooling
  and generated site meet their feature-specific acceptance checks without a
  production-token change; this does not add scientific validation evidence.

## Performance

- Scope and target: NOT APPLICABLE; no performance work was requested.
- Benchmark/profile environment: not run.
- Baseline and result: not applicable.
- Numerical-equivalence evidence: normalized C++ tokens are identical and all
  existing regression suites pass.
- API/ABI impact: none.
- Results by target environment: not applicable.
- Remaining bottlenecks: not assessed.

## Quality gates

| Gate | Configuration/command | Result | Evidence |
|---|---|---|---|
| Unit tests | Debug, Release, release-max CTest | PASS: 15/15 each | [verification.md](verification.md) |
| Line/function/branch coverage | docs excluded; source tokens unchanged | NOT APPLICABLE | Coverage table above |
| Coverage-exclusion audit | feature diff | PASS: no exclusions | Coverage table above |
| clang-tidy/static analysis | `debug-tidy`, clang-tidy 22.1.8 | FAIL: inherited 426 unique first-party warning lines plus 15 Enzyme compatibility warnings | [verification.md](verification.md) |
| Sanitizers | Debug ASan/LSan | PASS: 15/15 | [verification.md](verification.md) |
| ABI compatibility | header-only, token-identical | NOT APPLICABLE | [verification.md](verification.md) |
| Adversarial tooling | `npm test --prefix docs`; independent skeptic | PASS | [verification.md](verification.md) |
| Integration/V&V | existing numerical regression suites | PASS as regression; new V&V NOT APPLICABLE | [verification.md](verification.md) |
| Documentation | clean CMake docs build and exact output checker | PASS | [verification.md](verification.md) |
| Formatting | repository clang-format dry run | FAIL: inherited 89 diagnostics in 14 files | [verification.md](verification.md) |
| Performance regression | no executable-token change | NOT APPLICABLE | Performance section above |

## Local Git state

- Base revision: `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`
- Integration branch: `sourcey`
- Audited integration revision: `261f4bbabfe88f59d4de6b25feeab01ab6609da8`
- PR-ready local branch: NOT CREATED; final policy gates failed.
- Final squash commit: NOT CREATED.
- Squashed-tree equivalence check: NOT APPLICABLE.
- Remote status: not pushed; no PR opened.

## Reproduce

```sh
cmake --preset docs
cmake --build --preset docs

npm test --prefix docs
node docs/check-output.mjs

cmake --preset debug -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug --parallel
ctest --preset debug --output-on-failure

cmake --preset release -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release --parallel
ctest --preset release --output-on-failure

cmake --preset release-max -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release-max --parallel
ctest --preset release-max --output-on-failure

cmake --preset debug-tidy \
  -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme \
  -DCLANG_TIDY_EXE=/usr/bin/clang-tidy-22
cmake --build --preset debug-tidy --parallel

find include tests benchmarks -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
  -exec clang-format --dry-run --Werror {} +

git diff --check 9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc..sourcey
```

## Documentation

- API/user docs: `docs/*.md`, `docs/api-groups.dox`, and generated `docs/dist`.
- Architecture: [plan.md](plan.md).
- V&V: not applicable; regression evidence is in [verification.md](verification.md).
- Benchmarks: not applicable.

## Risks, limitations, tooling gaps, and follow-up

- Resolve the inherited clang-tidy and clang-format failures, or explicitly
  revise the repository policy, before creating a PR-ready squash.
- Add warning-as-error enforcement after the warning baseline is clean.
- Consider a documentation-example compilation gate, a complete external-link
  checker, and JavaScript coverage/mutation testing.
- TSan, MSan, property-based testing, fuzzing, and mutation-testing policies
  remain unconfigured repository-wide; UBSan remains prohibited with Enzyme.

## Reviewer guide

- [ ] Review the Doxygen XML filtering and Sourcey data flow in `docs/CMakeLists.txt`,
  `docs/filter-doxygen-xml.mjs`, and `docs/sourcey.config.ts`.
- [ ] Inspect the exact output/API checks in `docs/check-output.mjs` and the
  adversarial fixtures in `docs/test-render-math.mjs` and
  `docs/test-filter-doxygen-xml.mjs`.
- [ ] Review KaTeX postprocessing and collision protections in
  `docs/render-math.mjs`.
- [ ] Confirm legacy URL behavior in `docs/write-legacy-redirects.mjs`.
- [ ] Review the 17 header comment-only changes with the token-equivalence
  evidence in [verification.md](verification.md).
- [ ] Decide whether to resolve the inherited clang-tidy/clang-format debt or
  revise the zero-warning policy before a PR-ready squash is created.
- [ ] Inspect the full feature diff against
  `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`.
