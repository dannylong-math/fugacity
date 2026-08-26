#!/usr/bin/env python3
"""Require complete first-party LLVM source coverage across test executables.

Header-only templates are instantiated differently in each independently linked
test executable. Merging their profiles before asking llvm-cov to interpret the
mapping loses data as "mismatched functions". This checker exports each binary
against only its own profile, then unions JSON function/branch source sites and
llvm-cov's canonical lcov line records across the per-executable reports.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


Site = tuple[object, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--llvm-profdata", required=True)
    parser.add_argument("--llvm-cov", required=True)
    parser.add_argument("--expected-lines", type=int, default=1081)
    parser.add_argument("--expected-functions", type=int, default=151)
    parser.add_argument("--expected-branch-outcomes", type=int, default=266)
    return parser.parse_args()


def run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[Any]:
    return subprocess.run(command, check=True, **kwargs)


def production_path(filename: str, production_root: Path) -> Path | None:
    path = Path(filename).resolve()
    try:
        path.relative_to(production_root)
    except ValueError:
        return None
    return path


def native_denominator_is_zero(file_record: dict[str, Any]) -> bool:
    summary = file_record["summary"]
    return all(summary[metric]["count"] == 0 for metric in ("lines", "functions", "branches"))


def add_count(store: dict[Site, int], key: Site, count: int) -> None:
    store[key] = store.get(key, 0) + count


def add_branch(store: dict[Site, list[int]], key: Site, branch: list[Any]) -> None:
    totals = store.setdefault(key, [0, 0])
    totals[0] += branch[4]
    totals[1] += branch[5]


def collect_report(
    report: dict[str, Any],
    production_root: Path,
    functions: dict[Site, int],
    branches: dict[Site, list[int]],
) -> None:
    data = report["data"][0]
    for file_record in data["files"]:
        path = production_path(file_record["filename"], production_root)
        if path is None:
            continue
        if native_denominator_is_zero(file_record):
            continue
        for branch in file_record["branches"]:
            key = (str(path), "direct", branch[0], branch[1], branch[2], branch[3])
            add_branch(branches, key, branch)
        for expansion in file_record["expansions"]:
            call_site = expansion["source_region"][0:4]
            for branch_index, branch in enumerate(expansion["branches"]):
                key = (str(path), "expansion", *call_site, branch_index)
                add_branch(branches, key, branch)

    for function in data["functions"]:
        if not function["regions"]:
            continue
        primary = function["regions"][0]
        file_id = primary[5]
        if file_id >= len(function["filenames"]):
            continue
        path = production_path(function["filenames"][file_id], production_root)
        if path is None:
            continue
        # All first-party declarations place one source function per line.
        # This collapses template instantiations and the alternate zero-count
        # constructor mapping that begins at its member-initializer region.
        key = (str(path), primary[0])
        add_count(functions, key, function["count"])


def collect_lcov_lines(report: Path, production_root: Path, lines: dict[Site, int]) -> None:
    current_file: Path | None = None
    with report.open(encoding="utf-8") as source:
        for record in source:
            if record.startswith("SF:"):
                current_file = production_path(record[3:].strip(), production_root)
            elif current_file is not None and record.startswith("DA:"):
                line, count, *_ = record[3:].split(",")
                add_count(lines, (str(current_file), int(line)), int(count))


def ratio(covered: int, total: int) -> str:
    percent = 100.0 if total == 0 else (100.0 * covered / total)
    return f"{covered}/{total} ({percent:.2f}%)"


def main() -> int:
    args = parse_args()
    build_dir = args.build_dir.resolve()
    production_root = (args.source_root.resolve() / "include" / "fugacity").resolve()
    binaries = sorted(
        path
        for path in (build_dir / "tests").glob("test_*")
        if path.is_file() and os.access(path, os.X_OK)
    )
    if not binaries:
        raise RuntimeError(f"no test executables found under {build_dir / 'tests'}")

    profraw_dir = build_dir / "profraw"
    profdata_dir = build_dir / "profdata"
    json_dir = build_dir / "coverage-json"
    lcov_dir = build_dir / "coverage-lcov"
    for directory in (profraw_dir, profdata_dir, json_dir, lcov_dir):
        directory.mkdir(parents=True, exist_ok=True)

    line_counts: dict[Site, int] = {}
    function_counts: dict[Site, int] = {}
    branch_counts: dict[Site, list[int]] = {}

    for binary in binaries:
        raw_profile = profraw_dir / f"{binary.name}.profraw"
        indexed_profile = profdata_dir / f"{binary.name}.profdata"
        json_report = json_dir / f"{binary.name}.json"
        lcov_report = lcov_dir / f"{binary.name}.lcov"
        for artifact in (raw_profile, indexed_profile, json_report, lcov_report):
            artifact.unlink(missing_ok=True)

        environment = os.environ.copy()
        environment["LLVM_PROFILE_FILE"] = str(raw_profile)
        run([str(binary)], cwd=build_dir, env=environment)
        run([args.llvm_profdata, "merge", "-sparse", str(raw_profile), "-o", str(indexed_profile)])
        with json_report.open("w", encoding="utf-8") as output:
            run(
                [
                    args.llvm_cov,
                    "export",
                    str(binary),
                    f"-instr-profile={indexed_profile}",
                    "-ignore-filename-regex=(_deps|/tests/).*",
                ],
                stdout=output,
            )
        with json_report.open(encoding="utf-8") as source:
            report = json.load(source)
        with lcov_report.open("w", encoding="utf-8") as output:
            run(
                [
                    args.llvm_cov,
                    "export",
                    str(binary),
                    f"-instr-profile={indexed_profile}",
                    "-format=lcov",
                    "-ignore-filename-regex=(_deps|/tests/).*",
                ],
                stdout=output,
            )
        totals = report["data"][0]["totals"]
        print(
            f"native {binary.name}: "
            f"lines {ratio(totals['lines']['covered'], totals['lines']['count'])}, "
            f"functions {ratio(totals['functions']['covered'], totals['functions']['count'])}, "
            f"branches {ratio(totals['branches']['covered'], totals['branches']['count'])}"
        )
        collect_lcov_lines(lcov_report, production_root, line_counts)
        collect_report(report, production_root, function_counts, branch_counts)

    covered_line_count = sum(count > 0 for count in line_counts.values())
    covered_function_count = sum(count > 0 for count in function_counts.values())
    covered_branch_outcomes = sum(
        (true_count > 0) + (false_count > 0) for true_count, false_count in branch_counts.values()
    )
    branch_outcomes = 2 * len(branch_counts)
    print("source-site union:")
    print(f"  lines: {ratio(covered_line_count, len(line_counts))}")
    print(f"  functions: {ratio(covered_function_count, len(function_counts))}")
    print(f"  branch outcomes: {ratio(covered_branch_outcomes, branch_outcomes)}")

    failures: list[str] = []
    failures.extend(f"uncovered line: {key[0]}:{key[1]}" for key, count in line_counts.items() if count == 0)
    failures.extend(
        f"uncovered function: {key[0]}:{key[1]}"
        for key, count in function_counts.items()
        if count == 0
    )
    for key, (true_count, false_count) in branch_counts.items():
        if true_count == 0:
            failures.append(f"uncovered true outcome: {key}")
        if false_count == 0:
            failures.append(f"uncovered false outcome: {key}")
    if len(line_counts) != args.expected_lines:
        failures.append(f"line denominator changed: expected {args.expected_lines}, observed {len(line_counts)}")
    if len(function_counts) != args.expected_functions:
        failures.append(
            f"function denominator changed: expected {args.expected_functions}, observed {len(function_counts)}"
        )
    if branch_outcomes != args.expected_branch_outcomes:
        failures.append(
            f"branch denominator changed: expected {args.expected_branch_outcomes}, observed {branch_outcomes}"
        )

    if failures:
        print("coverage failures:", file=sys.stderr)
        for failure in sorted(failures):
            print(f"  {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
