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
| PERF-1 | none | C++ performance engineer | audit branch/worktree | in progress | Baselines, profiles, model matrix, assembly, and candidate/no-change recommendation |
| ARCH-1 | PERF-1 | software architect | integration branch, read-only source review | pending | Internal design, compatibility and numerical-ordering assessment |
| IMPL-1 | ARCH-1 | implementation engineer | dedicated task branch/worktree | pending if needed | Red equivalence test, minimal source change, focused/regression tests, commit SHA |
| TEST-1 | IMPL-1 | test skeptic | audited implementation tree | pending if needed | Fault hypotheses, adequacy findings, coverage/tooling gaps |
| VV-1 | IMPL-1 | V&V scientist | audited implementation tree | pending if needed | Predeclared code-verification criteria and numerical parity evidence |
| PERF-2 | IMPL-1 | C++ performance engineer | audited implementation tree | pending if needed | Paired before/after benchmark and causal code-generation explanation |
| QA-1 | accepted tree | quality-gate auditor | integration branch | pending | Independent final gate matrix and READY/NOT READY verdict |

