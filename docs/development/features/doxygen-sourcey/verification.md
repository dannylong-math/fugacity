# Doxygen and Sourcey migration verification

## Environment

- Feature base: `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`
- CMake 3.28.3
- Clang 22.1.8
- Doxygen 1.9.8
- Node.js 24.19.0 (project minimum: 22.12)
- npm 11.17.0
- Sourcey 3.6.5, exact lockfile pin
- KaTeX 0.18.4, exact lockfile pin
- Enzyme CMake package: `/opt/enzyme/lib/cmake/Enzyme`

## Replacement baseline and documentation build

Before implementation, the retained Sphinx wrapper was exercised with:

```console
cmake --preset docs
cmake --build --preset docs
```

Configuration stopped because `sphinx-build` was not installed. This captured
the old pipeline baseline before replacement.

The direct replacement pipeline was then run from a clean generated-output
state:

```console
cmake -E make_directory build/doxygen
doxygen Doxyfile
npm ci --prefix docs
npm run --prefix docs build
```

Result: pass. Doxygen emitted no warnings; the XML filter excluded 55 private
implementation members and removed 608 ambiguous member references with zero
remaining orphan references. It also reassigned 51 duplicate namespace member
references to their public groups and attached the three public concepts.
Sourcey built 29 content pages; four compatibility redirects plus the generated
index produced 34 HTML artifacts. KaTeX rendered 43 display and 266 inline
equations. `npm audit` reported zero vulnerabilities during the original
migration build.

The stable project wrapper was also rebuilt after removing `build/docs`,
`build/doxygen`, and `docs/dist`:

```console
cmake --preset docs
cmake --build --preset docs
```

Result: pass. Doxygen 1.9.8 and Node.js 24.19.0 were detected, and the same page,
API, link, redirect, and equation inventory was generated under `docs/dist`.
`WARN_AS_ERROR=FAIL_ON_WARNINGS` gates every enabled Doxygen warning.
`WARN_NO_PARAMDOC` remains narrowly disabled because Doxygen 1.9.8 reports
false positives for documented constrained C++23 functions; parameter
documentation remains a review rule rather than a reliable 1.9.8 build gate.

The Pages workflow no longer installs npm dependencies separately: the CMake
`docs` target owns the single `npm ci` and then runs the checked npm build.

## Generated-site checks

The initial post-migration checker was permissive: it treated a missing
`api/name.html` target as successful when `api/name/index.html` existed. The
skeptic's exact-artifact oracle correctly found 879 missing targets among 1,673
links. This was a false pass in the initial verification, not a deployment-safe
site.

The correction uses Sourcey's flat `.html` output, supplies its missing
`api/index.html` compatibility target, and commits an exact checker that runs in
every npm/CMake/CI documentation build:

```console
node docs/check-output.mjs
```

Result: pass, 2,119 exact local `href`/`src` targets and HTML fragments checked
across 34 HTML pages, with no fallback from one URL layout to another. The same
gate verified all 50 public core function overloads, `ideal_gas_constant`, the
`EquationOfState`, `IdealEoS`, and `ResidualEoS` concept pages, the four legacy
URLs, absence of Edit-this-page links, and absence of `fugacity::detail` or
Enzyme declarations from generated API pages.

The XML filter's orphan-reference behavior and the parser-aware math renderer
have committed negative tests:

```console
npm test --prefix docs
```

Result: pass. An unmatched inline delimiter (`$x^2`), unmatched display
delimiter, invalid TeX, and a compound XML orphan all fail. Dollar signs inside
HTML attributes, `pre`, and `code` remain untouched. Valid escaped dollars in
`$\\$5$` and `$\\text{price \\$5}$` render successfully: rendered KaTeX is held
behind a temporary token until unmatched source delimiters have been checked.
Temporary tag, block, and rendered-math tokens use a prefix proven absent from
each input, and regression fixtures preserve the former literal
`FUGACITYPROTECTEDTAG0END` and `FUGACITYPROTECTEDBLOCK0END` sentinels. Formula
markers use unambiguous non-base64 framing, and KaTeX CSS/fonts are copied
locally.

Original browser inspection used a temporary localhost server. The tutorial and
Peng-Robinson API pages were inspected through the rendered DOM and screenshots.
Tables, code blocks, navigation, full public include paths, inheritance links,
inline equations, and display equations rendered. The Peng-Robinson page had
16 KaTeX equations (six display equations), no literal `$$` delimiters, and no
detail/Enzyme API text. Long equations remained within the content column.

A corrective browser run was completed on the integrated `sourcey` tree after
the exact-link, API, and math fixes. The tutorial rendered 19 KaTeX formulas
(five display formulas), 11 code blocks, and one table without literal `$$`
delimiters or horizontal overflow. The Core API rendered 51 member headings,
including four `calc_pressure` overloads and `ideal_gas_constant`, with 125
KaTeX formulas and no internal API leakage. The Peng-Robinson page rendered 16
KaTeX formulas (six display formulas), seven tables, its full public include
path, and the `BaseCubic` relationship without horizontal overflow. The browser
reported no warnings or errors on that representative model page.

Tracked-tree searches found no surviving Sphinx configuration, dependency, RST
source, or Sphinx workflow reference outside this feature's historical plan and
verification record.

## C++ regression gates

```console
cmake --preset debug -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset debug --parallel
ctest --preset debug --output-on-failure

cmake --preset release -DEnzyme_DIR=/opt/enzyme/lib/cmake/Enzyme
cmake --build --preset release --parallel
ctest --preset release --output-on-failure
```

Result: both configurations and builds passed; Debug and Release each passed
15/15 tests. The first sandboxed Debug CTest run reported LeakSanitizer's known
ptrace incompatibility after every test's assertions passed. Repeating CTest
outside that sandbox passed 15/15 with ASan/LSan enabled.

The expanded quality cleanup reran coverage with each of the 15 independently
linked test executables profiled and exported separately. The committed
`coverage-check` target unions resolved first-party source sites under
`include/fugacity` and reports 1031/1031 lines, 150/150 source functions, and
244/244 branch outcomes (122 conditions), with no exclusions. A temporary
uncovered branch in `core/horner.hpp` made the gate fail at 1032/1034 lines and
245/246 outcomes, proving fail-under enforcement. The mutant was removed.

The legacy merged native report remains informational. Before the focused
precondition tests it reported 1013/1013 lines, 150/150 functions, and 227/244
branch outcomes. After the tests it reports 1013/1013, 150/150, and 242/244,
while still warning that 270 functions have mismatched data. Merging unrelated
executables drops 18 mapped source lines and conflicting header-template branch
instances, so the merged report is not used as the completion gate.

## Source-token and repository checks

For every modified public header, Clang 22 raw tokens were dumped for the
feature tree and feature base. Comment and whitespace (`unknown`) tokens and
source locations were removed before byte comparison:

```console
clang++ -std=c++23 -Xclang -dump-raw-tokens -fsyntax-only <header>
```

Result: pass for all 17 modified headers. The normalized token streams were
identical, establishing that declarations and executable/source tokens did not
change. The corrective run repeated this with Clang's preprocessed token dump,
normalizing built-in source-location macros (`__FILE__` and `__LINE__`) so
documentation-only line movement could not create false differences.

```console
find include tests benchmarks -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
  -exec clang-format --dry-run --Werror {} +
git diff --check
```

`git diff --check` passed. The initial repository-wide formatting dry check
reported 89 diagnostics in 14 files. After formatting the authoritative
`include`, `tests`, and `benchmarks` trees, the same command reports zero.

The initial independent `debug-tidy` audit reported 426 unique first-party
warning lines plus 15 Enzyme analyzer-plugin compatibility warnings. The exact
disposition ledger is:

| Disposition | Count | Checks |
|---|---:|---|
| Approved global disables | 293 | 277 bounds-index findings; 16 `#pragma once` findings |
| Exact approved suppressions | 62 | 37 Enzyme casts; 15 test-main exception escapes; 6 fully assigned arrays; 2 retained getters; 1 public temperature macro; 1 CRTP constructor |
| Actionable fixes | 71 | 26 include-cleaner; 26 const-correctness; 8 use-auto; 6 redundant-typename; 3 concise-preprocessor; 1 redundant-parentheses; 1 multiprecision lifetime |

The final all-target `debug-tidy` build treats every enabled diagnostic as an
error and reports zero warnings and errors across all 15 test translation
units. Both benchmark translation units were analyzed explicitly and also
report zero. The real compiler retains Enzyme's plugin; only clang-tidy's false
analyzer-subprocess plugin warning is suppressed. Every `NOLINT` annotation is
check-qualified, and there are no bare suppressions.

## Applicability

The quality cleanup preserves the public source API, storage layout, numerical
expression order, Enzyme intrinsic spelling and behavior, exception semantics,
test tolerances, and deterministic seeds. A base/current consumer produced
identical layout values and hexadecimal numerical results. The only intentional
calculation-support change is the test oracle's owning multiprecision return
lifetime; it does not alter production calculations, data formats, model
validation, uncertainty quantification, or performance claims.
