# EOS kernel performance refinement

## Scope and constraints

- Audit the scalar CPU performance of every currently bundled ideal and residual
  equation-of-state implementation.
- Change production code only where profiling and a representative benchmark identify a
  clear restructuring opportunity.
- Introduce no new library or runtime dependency.
- Preserve the public source API, mathematical formulas, coefficient interpretation, SI
  conventions, compiler floating-point mode, and numerical behavior. Fugacity is
  header-only, so there is no standalone binary ABI to preserve.
- Do not add SIMD, GPU execution, alternative automatic differentiation, fast-math, or a
  new batch/traversal API.
- Use a throwaway benchmark for each accepted optimization and retain the commands and raw
  evidence in `performance.md`, not the benchmark source in the final production diff.
- If the initial audit finds performance reasonable and no clear source-level opportunity,
  leave the corresponding implementation unchanged.

## Context of use

The target is repeated scalar CPU evaluation of thermodynamic states using the current
Clang/Enzyme toolchain. A false optimization could silently change thermodynamic values or
derivatives, so correctness and floating-point behavior take precedence over small or
noise-level timing changes.

## Acceptance rules

An implementation change is accepted only when:

1. a baseline profile identifies the affected expression as material to runtime;
2. the change is an internal restructuring with no public API or dependency change;
3. focused and full numerical tests pass under the repository's standard floating-point
   profile;
4. a paired, repeated, same-toolchain benchmark shows a speedup outside measurement noise;
5. optimized IR/assembly explains the speedup without relying on relaxed arithmetic; and
6. the complete coverage, warning, formatting, sanitizer, and applicable V&V gates pass.

## Architecture decision

The initial audit found one material opportunity: the runtime-sized generalized-cubic
attractive mixing term recomputes a temperature/composition factor inside every matrix
interaction. The accepted design processes eight rows at a time and keeps eight independent
row accumulators on the stack.

- Apply the blocked path only to runtime-sized `double` models with at least 10 components.
- Preserve the existing scalar loop literally for static extents, smaller dynamic models,
  `float`, and `long double`.
- Keep each row's inner accumulation and the final row accumulation in their original order.
- Use one private helper only from the molar kernel; allocate no heap storage and mutate no
  model state. The density kernel retains its original expression graph because an early
  blocked prototype changed reverse-mode fugacity bits.
- Limit production changes to `include/fugacity/residual_models/cubic.hpp` and permanent
  characterization tests to `tests/test_peng_robinson.cpp`.

Alternatives rejected for this pass are a heap-allocated full cache, mutable model-owned
cache, caller workspace/API, changing the static loop, and applying blocking to `float`.
Although the float primal was bitwise identical, sampled Enzyme-derived float sound speed
and fugacity values rounded differently. The dynamic-double-only design therefore preserves
the project's current numerical-compatibility policy without requiring a new tolerance
decision.

## Protected unrelated work

- `docs/development/features/flash-calculations/plan.md` is pre-existing and untracked in
  the primary worktree. Initial SHA-256:
  `9a47b5514c08a272c6ca4ee2f9068018f622f5525f7e8cee69bfad44b75fabc6`.
- The `feature/eos-model-catalog` worktree is outside this feature and must not be edited.
  Initial tree: `c5fc08cc088b067ee67d87ddc600650151f01fa4`.
- The dirty legacy `performance` worktree is outside this feature and must not be edited or
  used as the source of production changes.

## Git state

- Base branch: `main`
- Base revision: `522c5879bf207be6c3679f27a7e080111dc20f70`
- Integration branch: `feature/eos-kernel-performance`
- Audit branch: `task/eos-kernel-performance/audit`
- Audit worktree:
  `/home/dannylong/research/eos_library/fugacity-eos-kernel-performance-audit`

## Tasks

| ID | Dependency | Owner | Branch/worktree | Status | Required evidence |
|---|---|---|---|---|---|
| PERF-1 | none | C++ performance engineer | audit branch/worktree | complete; no commit | Baselines, profiles, model matrix, assembly, and candidate/no-change recommendation |
| ARCH-1 | PERF-1 | software architect | integration branch, read-only source review | complete; no commit | Dynamic-double-only eight-row blocking; no API or allocation |
| IMPL-1 | ARCH-1 | implementation engineer | task branch/worktree | complete at `d15bec8`; follow-up pending | Characterization tests, source change, full regression, coverage |
| TEST-1 | IMPL-1 | test skeptic | review worktree | complete; no commit | Release READY; Debug/O1 nested-AD limitation; Release-only regression requested |
| VV-1 | IMPL-1 | V&V scientist | review worktree | complete; no commit | Release PASS against MP/finite differences; Debug limitation confirmed |
| PERF-2 | IMPL-1 | C++ performance engineer | audited implementation tree | complete; no commit | Paired production A/B and causal code-generation explanation |
| QA-1 | accepted tree | quality-gate auditor | integration branch | pending | Independent final gate matrix and READY/NOT READY verdict |

## Performance decision record

The frozen molar-only candidate materially improves dynamic-double pressure, cp, and sound
speed, but Enzyme's derivative graph is not universally bitwise identical. The user
approved proceeding with bounded roundoff after the paired performance evidence was
presented.

## Numerical-policy decision

On 2026-08-26, the user approved the measured derivative roundoff in exchange for the
material dynamic Peng–Robinson speedup, then clarified that optimized Release is the
primary context of use because Enzyme is optimization-sensitive.

The focused legacy-characterization states retain these local regression checks:

- direct primal Helmholtz values, pressure, density-derived quantities, and fugacity remain
  bitwise identical in the characterization matrix;
- internal differentiated lambda values remain within 32 ULP of the legacy runtime loop;
- cp and sound speed squared remain within 16 ULP of the legacy runtime loop; and
- all independent scientific oracles and existing public tolerances remain unchanged.

The 32/16 ULP limits are local regression checks, not global error bounds. An expanded V&V
matrix exceeded them in raw ULP distance while remaining well inside independent scientific
tolerances. Release acceptance therefore depends on multiprecision and refined
finite-difference oracles for the public API. Debug/O1 must compile and pass the existing
ASan/test suite, but large-N nested composition gradients are not numerically qualified
because the installed Enzyme produces optimization-sensitive results there.
