# Fugacity project guidance

## Project overview

- Fugacity is a header-only C++23 library for calculating thermodynamic
  properties from Helmholtz-energy equations of state.
- Pair an ideal Helmholtz contribution with a residual contribution. Enzyme
  automatic differentiation supplies pressure, caloric properties, chemical
  potentials, fugacities, and public property derivatives.
- Bundled ideal models are constant heat capacity, NASA-7, and NASA-9.
  Bundled residual models are zero residual, van der Waals, and Peng-Robinson.
- Public thermodynamic inputs and outputs use SI units. Temperature is K,
  concentration and partial densities are mol/m^3, and composition is
  dimensionless.
- The library evaluates properties at supplied `(c, x, T)` states. It does not
  solve phase equilibrium or invert pressure and temperature, and callers
  remain responsible for model-specific physical-domain validity.
- Supported execution is workstation CPU: scalar reference,
  `std::experimental::simd`, and optional Google Highway. MPI, OpenMP, GPU,
  CUDA, HIP, SYCL, Kokkos, PETSc, and Trilinos are not supported backends.
- Ubuntu CPU with Clang/Enzyme is CI-verified. Linux/macOS Homebrew setup is
  documented, but macOS is not CI-verified.
- Preserve the public headers, `fugacity` namespace, `Fugacity::Fugacity`
  target, SI conventions, and numerical behavior by default. A formal
  source-API stability/versioning policy is TBD.
- Fugacity has no separately compiled library, install/export package, plugin
  ABI, or language binding. Source compatibility is therefore primary and no
  standalone binary ABI baseline currently applies.
- Benchmark JSON/CSV output and the benchmark registry are not persistent
  scientific serialization formats. Data/schema compatibility policy is TBD.

## Build and dependencies

- Run project commands from the repository root.
- CMake minimum version is 3.21 and the verified generator is Unix Makefiles.
- Clang is mandatory. Enzyme must be built for the same LLVM version.
- Configure, build, and test Debug with:

  ```console
  cmake --preset debug -DEnzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
  cmake --build --preset debug --parallel
  ctest --preset debug --output-on-failure
  ```

- Replace `debug` with `release` for optimized verification.
- `debug` deliberately uses `-O1 -g` plus ASan because Enzyme derivatives are
  not reliable at `-O0`.
- `release-max` permits reassociation, reciprocal transformations, and changed
  signed-zero behavior. Do not use it as the numerical correctness baseline.
  Run it only as an additional compatibility check unless relaxed arithmetic
  is explicitly approved for the change.
- Run clang-tidy during compilation with:

  ```console
  cmake --preset debug-tidy -DEnzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
  cmake --build --preset debug-tidy --parallel
  ctest --preset debug-tidy --output-on-failure
  ```

- Clean a configured preset with
  `cmake --build --preset <preset> --target clean`.
- No CMake workflow presets or repository Makefile targets are defined.
- Non-documentation configuration replaces the source-root
  `compile_commands.json` symlink. Do not edit it manually, and avoid
  concurrent configurations that race to repoint it.
- Test configuration may fetch Boost.UT, Boost.Config, and
  Boost.Multiprecision. Highway and Google Benchmark may also be fetched when
  their options are enabled. Configuration can therefore access the network
  and modify build directories; do not run it merely for discovery.

## Tests

- Boost.UT is the unit-test framework. Each `tests/*.cpp` file becomes a CTest
  executable named after its filename stem.
- Run one registered test with:

  ```console
  ctest --preset debug -R '^<exact_test_name>$' --output-on-failure
  ```

- Run all unit and numerical-verification tests with:

  ```console
  ctest --preset debug --output-on-failure
  ctest --preset release --output-on-failure
  ```

- For numerical kernels, Enzyme integration, SIMD, or optimization-sensitive
  changes, also run the applicable `release-max` and Highway-enabled matrices.
- Preserve deterministic test inputs and record seeds. Existing contract
  sweeps use seed `0xC0FFEE`; tests using `std::random_device` are a known
  reproducibility gap and should be converted to recorded deterministic seeds.
- Do not describe analytic identities, derivative oracles, teqp golden values,
  or coverage as physical model validation.

## Adversarial testing

- Property-based framework and command: TBD.
- Fuzzer targets, command, corpus, and reproducer locations: TBD.
- Mutation-testing tool and command: TBD.
- Do not install adversarial tools unilaterally. Report missing tooling and
  propose a version compatible with the required Clang/LLVM/Enzyme toolchain.
- Every adopted property or fuzz workflow must record its seed, input domain,
  corpus location, minimized reproducer convention, and bounded CI command.

## Quality gates

- First-party coverage scope is `include/fugacity/**`. Tests, fetched
  dependencies, documentation, and benchmarks are outside the production
  denominator.
- Required line coverage: 100%.
- Required function coverage: 100%.
- Required branch coverage: 100%.
- Configure and run the current coverage workflow with:

  ```console
  cmake --preset coverage -DEnzyme_DIR=/path/to/Enzyme/lib/cmake/Enzyme
  cmake --build --preset coverage --parallel
  mkdir -p build/coverage/profraw
  LLVM_PROFILE_FILE="$PWD/build/coverage/profraw/%p.profraw" \
    ctest --preset coverage --output-on-failure
  ```

- Use the `llvm-profdata merge` and `llvm-cov export/report` sequence in
  `.github/workflows/ci.yml`. Ensure the matching LLVM toolchain is on `PATH`;
  local installations may expose only versioned names such as `llvm-cov-22`
  and `llvm-profdata-22`.
- Current Codecov checks are informational and are not the completion gate.
  Automated fail-under enforcement remains a tooling gap; record raw
  line/function/branch metrics in each feature verification report.
- No coverage exclusion is approved by default. Record proposed exclusions in
  `docs/development/quality/coverage-exclusions.md` with a stable ID, exact
  source locations, proof of unreachability, tool-specific suppression,
  reviewer, approval date, and reevaluation trigger. Only the user may approve
  an exclusion.
- Coverage is a gate, not evidence that tests contain adequate scientific
  oracles or properties.
- First-party compiler and clang-tidy warnings must be zero. The build currently
  lacks `-Werror` and clang-tidy warnings-as-errors; report this enforcement
  gap until it is closed.
- Debug ASan is the verified sanitizer gate.
- UBSan is prohibited for Enzyme-differentiated targets until compatible
  evidence exists. TSan and MSan policy is TBD.
- Check formatting without modifying files with:

  ```console
  find include tests benchmarks -type f \
    \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
    -exec clang-format --dry-run --Werror {} +
  ```

- Apply formatting with `./format.sh`. It rewrites files and currently omits
  `benchmarks/`; update it when benchmark sources become in-scope.

## Documentation and development records

- The first managed task after this bootstrap is to replace Sphinx with Doxygen
  XML plus Sourcey, using Rift as the local reference implementation.
- Pin Sourcey 3.6.5 in `docs/package.json` and `docs/package-lock.json`; require
  Node.js 22.12 or newer.
- Target documentation build:

  ```console
  cmake -E make_directory build/doxygen
  doxygen Doxyfile
  npm ci --prefix docs
  npm run --prefix docs build
  ```

- Retain `cmake --preset docs` and `cmake --build --preset docs` as the stable
  project-facing wrapper around that build.
- Generate Doxygen XML under `build/doxygen/xml` and the Sourcey site under
  `docs/dist`. Update GitHub Pages to publish `docs/dist`.
- Convert narrative RST pages to Markdown and preserve the current public API
  scope, including exclusion of `fugacity::detail`.
- Make Doxygen warnings fatal. Verify Sourcey output, navigation, API coverage,
  internal links, mathematical notation, code examples, and representative
  rendered pages before removing Sphinx.
- Do not delete the Sphinx configuration or dependencies until the replacement
  passes locally and in the documentation CI job.
- Do not hand-edit generated Doxygen XML, `docs/node_modules/`, `docs/dist/`,
  `docs/api/generated/`, or CMake-generated benchmark headers.
- Store feature records under `docs/development/features/<slug>/`, using only
  applicable files from `plan.md`, `verification.md`, `vv.md`,
  `performance.md`, and `pr-report.md`.
- Architecture/ADR convention: TBD.

## Performance

- Google Benchmark v1.9.4 is the benchmark framework. The
  `bench_calculations` launcher selects ideal/residual models, validates each
  backend against the scalar reference, then benchmarks the public `calc_*`
  functions.
- The launcher configures, builds, may fetch dependencies, and writes results;
  it is not a read-only discovery command.
- Representative smoke command:

  ```console
  ./bench_calculations \
    --ideal=const_cp --residual=none \
    --components=1 --batch=16 \
    --no-highway \
    --benchmark_repetitions=10 \
    --benchmark_report_aggregates_only=true
  ```

- Benchmark inputs are deterministic synthetic states, not physical validation
  mixtures.
- Record CPU model, ISA, core allocation, compiler/LLVM/Enzyme versions, build
  profile, backend, component and batch sizes, repetitions, system load,
  frequency/governor controls, and raw output.
- Stored baseline, regression threshold, and canonical benchmark report
  location: TBD.
- Do not compare `release` and `release-max` results as equivalent numerical
  configurations.

## V&V and literature

- Existing verification assets include analytic thermodynamic identities,
  multiprecision adaptive finite-difference derivative oracles, ideal-gas and
  critical-point identities, static/dynamic and batch-backend parity, NASA
  reference calculations, and teqp 0.23.1 golden values.
- Context-of-use convention: TBD.
- Physical validation datasets, provenance, and uncertainty records: TBD.
- Calculation/solution convergence command: not applicable to the current
  pointwise library; revisit if iterative solvers or discretizations are added.
- V&V report location: `docs/development/features/<slug>/vv.md`.
- Adopted V&V terminology or standard: TBD.
- Bibliography format/file and Zotero collection/tag convention: TBD.
- Record the provenance and reproduction procedure for coefficient datasets
  and external golden values before treating them as validation evidence.

## Local Git workflow

- Base branch: `main`.
- Use `feature/<slug>` for managed integration branches,
  `task/<slug>/<work-package>` for task branches, and `fix/<slug>` for small
  fixes.
- Use Conventional Commits with a meaningful scope.
- WIP commits are allowed only on task branches. The final local squash commit
  must summarize scientific rationale and verification evidence.
- Preserve integration and task history until the user accepts the squashed
  branch.
- Never push, create a remote branch, open or modify a pull request, deploy an
  artifact, trigger a remote workflow, or upload benchmark/coverage results.
  The user owns remote operations.

## Change policy

- Preserve public source API and numerical behavior unless a breaking change is
  explicitly approved.
- Preserve precision, tolerances, stopping criteria, deterministic seeds,
  coefficient interpretation, mixing rules, SI conventions, and normal
  floating-point ordering unless scientifically justified and approved.
- Dependency changes require approval and coordinated version pins, lockfiles,
  README instructions, and CI updates.
- Generated files must be reproducible and must not be edited manually.
- Commands that modify state include configuration, builds, tests that write
  profiles, dependency fetching, `npm ci`, documentation generation,
  benchmarks, formatting, coverage upload, Pages deployment, and remote Git
  operations.
