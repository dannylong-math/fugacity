# EOS kernel performance report

## First-pass audit

The audit used exact base revision
`522c5879bf207be6c3679f27a7e080111dc20f70`, Clang 22.1.8, the installed
ClangEnzyme-22 plugin, C++23, and `-O3 -march=native -DNDEBUG` without relaxed
floating-point flags. Runs were pinned to one performance core of an Intel Core i7-1370P.
The governor was `powersave`, energy preference was `balance_performance`, and turbo was
enabled, so relative paired measurements are more reliable than absolute timings.

The broad baseline used nine repetitions with a 20 ms minimum calibration. Focused
candidate measurements used 13 repetitions with a 50 ms minimum calibration. Hardware
performance counters were unavailable because `perf_event_paranoid=4`.

### Raw Helmholtz baseline

Median CPU time in nanoseconds:

| Model | Static N2 | Dynamic N2 | Static N10 | Dynamic N10 | Dynamic N50 |
|---|---:|---:|---:|---:|---:|
| ConstantCp | 13.9 | 16.7 | 46.1 | 54.8 | 267.8 |
| NASA7 | 14.1 | 19.5 | 52.0 | 61.7 | 309.1 |
| NASA9 | 14.4 | 19.0 | 58.2 | 67.4 | 298.6 |
| van der Waals | 4.12 | 6.18 | 15.9 | 26.7 | 944.9 |
| Peng–Robinson | 9.32 | 12.3 | 25.3 | 65.0 | 1624.9 |

Most coefficients of variation were 0.3–3.5 percent. The public dynamic-N50 timings
illustrate the additional differentiated work:

| Residual model | Pressure | cp | Sound speed squared | Fugacity |
|---|---:|---:|---:|---:|
| van der Waals | 0.907 us | 3.100 us | 5.663 us | 1.962 us |
| Peng–Robinson | 1.623 us | 9.871 us | 18.459 us | 3.363 us |

### No-change decisions

- ConstantCp already hoists its temperature/logarithm invariants and scales linearly.
- NASA7 and NASA9 use preprocessed coefficients, hoisted temperature terms, and scalarized
  Horner evaluation.
- NoResidual constant-folds to roughly 0.63–0.94 ns in the sampled public calculations.
- van der Waals performs one necessary row-major quadratic form. Exploiting symmetry would
  change the floating-point summation order.
- Static Peng–Robinson specializes well: Clang can unroll and eliminate repeated work when
  the component count is known.
- Reusing differentiated quantities in sound-speed calculations changed output bits in the
  first experiment and was rejected under the numerical-compatibility rule.

## Accepted candidate: dynamic-double cubic row blocking

The dynamic Peng–Robinson attractive mixing term contains a species factor of the form

```text
weight[j] * abs(p[j] - q[j] * sqrt(T))
```

that depends on `j` but not on the matrix row `i`. The current runtime loop calculates that
factor in every `(i,j)` interaction. Processing eight rows together calculates each `j`
factor once for those eight rows and updates eight independent row accumulators. Each row
still receives terms in increasing `j` order, and completed rows are still added in
increasing `i` order.

This is faster for two related reasons:

1. expensive multiply/absolute-value factor work falls from approximately N-squared
   evaluations to about N-squared divided by eight, plus N row factors; and
2. eight independent accumulators expose instruction-level parallelism. Under the strict
   floating-point profile Clang unrolls them into independent scalar FMA chains and reuses
   the computed column factor; it does not pack-vectorize or reassociate an individual
   reduction.

The design uses two fixed eight-element stack arrays and no heap allocation. It remains an
O(N-squared) generalized-cubic mixing calculation; it reduces the constant work inside the
same algorithm rather than changing the model.

### Preliminary throwaway result

These results use a surrogate copied from the production expression and establish the
candidate, not the final acceptance proof:

| Type | N | Current | Blocked-8 | Speedup |
|---|---:|---:|---:|---:|
| double | 10 | 61.5 ns | 35.0 ns | 1.76x |
| double | 20 | 201 ns | 95.3 ns | 2.11x |
| double | 50 | 1444 ns | 533 ns | 2.71x |
| double | 100 | 6522 ns | 1995 ns | 3.27x |
| float | 10 | 138 ns | 95.5 ns | 1.45x |
| float | 20 | 698 ns | 312 ns | 2.24x |
| float | 50 | 4966 ns | 1660 ns | 2.99x |
| float | 100 | 26519 ns | 6647 ns | 3.99x |

N2 regressed by roughly 8–13 percent. The accepted production dispatch therefore retains
the existing loop below N10. Only dynamic `double` selects blocking; static extents, `float`,
and `long double` remain unchanged. The float speedup is not accepted because sampled
Enzyme-derived float sound-speed and fugacity results rounded differently despite identical
primal values.

## Frozen production-candidate A/B

The candidate was narrowed to the molar Helmholtz kernel after blocking the density kernel
changed reverse-mode fugacity bits. The density kernel and all chemical-potential/fugacity
paths therefore retain the original expression graph.

The final throwaway benchmark compiled one identical harness against the exact base and
the frozen uncommitted candidate. It used 15 alternating, bracketed process pairs on one
core, at least 0.1 seconds calibration per case, and bootstrap confidence intervals for the
paired median speedup.

| Calculation | N10 speedup (95% interval) | N50 speedup (95% interval) |
|---|---:|---:|
| Raw molar Helmholtz | 1.50x [1.43, 1.57] | 2.91x [2.63, 2.96] |
| Pressure | 1.48x [1.34, 1.51] | 2.87x [2.74, 2.95] |
| cp | 1.41x [1.31, 1.44] | 2.11x [1.96, 2.18] |
| Sound speed squared | 1.38x [1.18, 1.46] | 2.03x [1.88, 2.12] |
| Fugacity | 1.03x [0.96, 1.08] | 0.98x [0.92, 0.99] |

Raw density, small-N, static, and float controls stayed within the noisy baseline-control
envelope. The N50 fugacity result shows an approximately two-percent apparent regression,
but that code path is source-identical and the interval is consistent with observed process
drift rather than a causal production change.

Strict-FP assembly contains eight independent scalar FMA accumulator chains. The helper
has a 216-byte stack frame and no allocation. The deepest sampled dynamic Enzyme derivative
frame increased from 152 to 1080 bytes, while the measured differentiated calculations
still improved. Binary text grew by 8856 bytes, about 1.55 percent.

### Numerical decision gate

A deterministic sweep covered 125 states at N=10, 16, 20, 50, and 100, temperatures from
180 to 900 K, concentrations from 5 to 1000 mol/m3, and 6025 scalar outputs.

- Molar/density primal values, pressure, lambda(0,1), lambda(0,2), and all 4900 fugacity
  outputs were bitwise identical.
- lambda(1,1) differed by at most 3 ULP.
- lambda(2,0) differed by at most 20 ULP, with maximum relative difference
  `3.98e-15`.
- cp differed in 1 of 125 states, by 1 ULP.
- sound speed squared differed in 4 of 125 states, by at most 2 ULP.
- The permanent-test states found a wider cp/sound envelope of 8 and 6 ULP under Release.

The candidate therefore passes the performance gate but is not accepted until the user
decides whether this derivative roundoff envelope is compatible with the project's
numerical policy.

## First-pass artifact provenance

| Artifact | SHA-256 |
|---|---|
| `/tmp/fugacity-eos-perf-audit-522c587/models_baseline.json` | `6e6c664b0f77dce77802f38672d24a2c3edc271d0e64e2e1cb1bcc7e0a23b52f` |
| `/tmp/fugacity-eos-perf-audit-522c587/bench_models.cpp` | `7a7b9dde61964c267982e22975c734f5742aa7084432c4175d6043d4f42ab13a` |
| `/tmp/fugacity-eos-perf-audit-522c587/bench_pr_cache.cpp` | `c5c3f1f16052220ab2e798b053caef6ee9b7d228cdd9f52f82089e2ebc2faf63` |
| `/tmp/fugacity-eos-perf-audit-522c587/pr_block_probe.json` | `1beaf404c43425cc7f3768716f1bf24968a5a5def83622c79f9e54a7288b0e1a` |
| `/tmp/fugacity-eos-perf-audit-522c587/pr_block_float_probe.json` | `e4e3e2b2896df3dd0074e0116a514a305e03d69bc3c27d7531bc07775a876b7d` |
| `/tmp/fugacity-eos-perf-audit-522c587/pr_vectorization_remarks.txt` | `195ec95de77484e6aca712658b27a716c543eb27dd4413e9412d443e8ada725e` |
| `/tmp/fugacity-eos-perf-audit-522c587/pr_cache_objdump.txt` | `9b59b653dd6bbc3edf02a2eefce3ceed1264945183dfdd4b430624ac246441ca` |
| `/tmp/fugacity-eos-perf-audit-522c587/enzyme_block_smoke.stdout` | `4613867d3f098b7e3defd67884f9ce45767943baf6205793eda337688b4a2ae4` |
| `/tmp/fugacity-eos-perf-ab-522c587/PERF2_REPORT.md` | `a8e5d8ec5847b6f84fdc3a54747e5c4bfa9722890cdfdb8b19dab89f2a3bc33b` |
| `/tmp/fugacity-eos-perf-ab-522c587/paired_runs.tar` | `1d75cb2530dd3d7b9e8135df134baf4acc9c51e33e9e695be51c97a43c75571e` |
| `/tmp/fugacity-eos-perf-ab-522c587/paired_summary.csv` | `14d20cd75146457ae7a2b07b1931ce19e26bb03f5af1bc7c918edac13a2905fc` |
| `/tmp/fugacity-eos-perf-ab-522c587/paired_raw.csv` | `b2605121278b32b82d4857f98b2e4645a4a9a539819886ea13355dd003245721` |
| `/tmp/fugacity-eos-perf-ab-522c587/ulp_summary.csv` | `fdcce2b573b5befb26796d635a73e6eec95aa07be246daae50f10f0f9a3ac1f6` |
| `/tmp/fugacity-eos-perf-ab-522c587/candidate_cubic.patch` | `c8544cadebf02e069b43f85085539cabbca7001e3b5eae085cc0eaf3322333ec` |
