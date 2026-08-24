# Doxygen and Sourcey documentation migration

## Scope

Replace the Sphinx/Sphinx-Immaterial documentation pipeline with Doxygen XML
and Sourcey 3.6.5. Preserve the existing `docs` CMake configure/build preset as
the project-facing entry point. Convert retained narrative documentation from
reStructuredText to Markdown and publish the generated Sourcey site from
`docs/dist`.

The migration must not change public declarations, source API, ABI,
thermodynamic equations, numerical behavior, precision, tolerances,
determinism, or performance.

## Base and compatibility

- Integration branch: `sourcey`.
- Feature base: `9b6933f5c9bbc99e42cfb09758c5f8cb31189ecc`.
- Original project base: `33458c05b1f8abe25c7d85504d40c29386cb9107`.
- Distribution mode: header-only; no separately compiled ABI.
- Stable documentation entry point:
  `cmake --preset docs && cmake --build --preset docs`.
- Toolchain: Doxygen 1.9.8-compatible configuration, Sourcey 3.6.5, and
  Node.js 22.12 or newer.
- The isolated `performance` worktree is outside this feature and must remain
  untouched.

## Architecture decision

Use a CMake-owned two-stage pipeline:

1. Doxygen scans first-party public headers and writes XML to
   `build/doxygen/xml`. Doxygen HTML output is disabled and warnings are fatal.
2. Sourcey consumes that XML plus explicit Markdown navigation and writes the
   static site to `docs/dist`.
3. `docs/CMakeLists.txt` keeps the `docs` target and invokes the reproducible
   lockfile-based Sourcey build.
4. GitHub Actions builds documentation on pull requests and pushes to `main`,
   and publishes `docs/dist` only for pushes.
5. Generated Doxygen XML, `docs/node_modules`, and `docs/dist` are ignored and
   never edited manually.
6. Sourcey emits flat `.html` artifacts so its generated Doxygen links resolve
   on a static host. Build-time compatibility redirects retain the former
   `getting_started.html`, `tutorial.html`, `implementing_a_new_eos.html`, and
   `api/concepts.html` entry points and bridge Sourcey 3.6.5's `api/index.html`
   navigation target to its emitted `api.html` overview.
7. The npm build renders TeX with exactly pinned KaTeX 0.18.4, then fails on
   malformed or unmatched delimiters, orphan XML references, missing public API
   inventory, internal API leakage, and any non-existent exact static target or
   fragment. HTML attributes and code-like elements are excluded from TeX
   parsing.

A standalone shell-only pipeline was rejected because it would bypass the
project's stable CMake preset/target and diverge from the approved Rift
reference workflow.

## Acceptance criteria

- `sourcey` is pinned exactly to 3.6.5 in a committed npm lockfile.
- `katex` is pinned exactly to 0.18.4 in the same lockfile.
- A clean `docs` configure/build produces Doxygen XML and the Sourcey site.
- Doxygen and Sourcey complete without warnings or errors.
- The public API reference excludes `fugacity::detail` and does not expose
  Enzyme implementation declarations as user API.
- Every retained Sphinx prose page has an equivalent navigable Markdown page.
- Representative equations, code blocks, API signatures, cross-links, and
  external links render correctly in a browser.
- README documentation links and local build instructions match the new site.
- No Sphinx dependency, configuration, generated page, workflow step, or RST
  source remains after replacement verification.
- Generated artifacts are ignored by Git.
- The normal Debug and Release library configure/build/test workflows remain
  green, establishing that the documentation-only migration did not affect
  executable behavior.
- C++ formatting and `git diff --check` pass.

## Task plan

| ID | Dependency | Owner | Branch/worktree | Status | Acceptance evidence |
| --- | --- | --- | --- | --- | --- |
| A1 | — | scientific software architect | `sourcey` | Complete | Architecture decision and compatibility risks recorded above |
| I1 | A1 | scientific implementation engineer | `task/doxygen-sourcey/migration` | Complete | Doxyfile, exact pins, Markdown, CMake/CI, strict output checks, legacy redirects, and clean documentation build |
| S1 | I1 | test skeptic | integration tree | Re-review pending | Initial audit found URL, API inventory, math, redirect, edit-link, and checking gaps; corrective implementation is ready for re-review |
| Q1 | S1 | quality gate auditor | integration tree | Pending | Independent documentation, build/test, formatting, and repository audit |
| R1 | Q1 | project manager | PR-ready branch | Pending | Squashed commit and evidence-based PR report |

## Risks and controls

- Existing header comments contain Sphinx roles and directives. Convert only
  documentation text; verify production declarations and executable tokens do
  not change.
- Doxygen may expose internal namespaces or compiler-plugin declarations.
  Configure exclusions and inspect the generated API inventory.
- Sourcey changes page URLs and anchors. Update repository links and preserve
  legacy paths or redirects where practical.
- Doxygen and Sourcey may interpret equations differently. Inspect
  representative inline and display mathematics visually.
- Documentation dependencies must remain optional: normal consumers and
  non-doc builds must not require Node.js, npm, Doxygen, or Sourcey.

## V&V and performance applicability

This feature changes documentation tooling and prose only. Physical model
validation, calculation/solution verification, uncertainty quantification,
and performance optimization are not applicable. Existing unit and numerical
verification tests are regression gates, not new V&V evidence.
