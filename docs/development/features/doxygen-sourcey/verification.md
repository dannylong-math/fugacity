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
remaining orphan references. Sourcey built 26 content pages, and KaTeX rendered
41 display plus 138 inline equations across 27 HTML files. `npm audit` reported
zero vulnerabilities.

The stable project wrapper was also rebuilt after removing `build/docs`,
`build/doxygen`, and `docs/dist`:

```console
cmake --preset docs
cmake --build --preset docs
```

Result: pass. Doxygen 1.9.8 and Node.js 24.19.0 were detected, and the same page
and equation inventory was generated under `docs/dist`.

## Generated-site checks

Representative guide and API outputs were asserted to exist, including
`introduction`, `getting-started`, `tutorial`, `implementing-a-new-eos`,
`concepts`, `fugacity::EoS`, `fugacity::Nasa7`, and
`fugacity::PengRobinson`. Searches of generated HTML and Doxygen compound XML
found no `fugacity::detail` compound or Enzyme declaration exposed as public
API.

An independent local static-link script checked generated `href` targets and
fragments:

```console
node /tmp/check_sourcey_links.mjs docs/dist
```

Result: pass, 1,673 internal links checked across 27 HTML pages. No repository
or generated-file link checker was added to the product.

Browser inspection used a temporary localhost server. The tutorial and
Peng-Robinson API pages were inspected through the rendered DOM and screenshots.
Tables, code blocks, navigation, full public include paths, inheritance links,
inline equations, and display equations rendered. The Peng-Robinson page had
16 KaTeX equations (six display equations), no literal `$$` delimiters, and no
detail/Enzyme API text. Long equations remained within the content column.

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

Coverage was not rerun because executable tokens are unchanged and coverage's
first-party denominator excludes documentation. Existing Debug and Release
tests are regression evidence, not new coverage or physical validation.

## Source-token and repository checks

For every modified public header, Clang 22 raw tokens were dumped for the
feature tree and feature base. Comment and whitespace (`unknown`) tokens and
source locations were removed before byte comparison:

```console
clang++ -std=c++23 -Xclang -dump-raw-tokens -fsyntax-only <header>
```

Result: pass for all 17 modified headers. The normalized token streams were
identical, establishing that declarations and executable/source tokens did not
change.

```console
find include tests benchmarks -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
  -exec clang-format --dry-run --Werror {} +
git diff --check
```

`git diff --check` passed. The repository-wide formatting dry check reported
89 diagnostics in 14 files. Repeating the same command against the feature base
reported the same 89 diagnostics with the same per-file counts; the only
modified header among those files is `assertions.hpp`, whose three diagnostics
are unchanged pre-existing macro-layout findings. No formatting rewrite was
made because this documentation-only task may not alter executable source.

## Applicability

The migration changes documentation tooling, prose, and comments only. ABI,
numerical behavior, scientific data formats, calculation/solution verification,
physical-model validation, uncertainty quantification, and performance are
unchanged. New numerical V&V and benchmarking are therefore not applicable.
