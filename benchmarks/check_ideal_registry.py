#!/usr/bin/env python3
"""Check that the ideal benchmark's grouped default sweep is complete."""

from __future__ import annotations

import argparse
from collections import Counter
import itertools
from pathlib import Path
import re
import subprocess


SIZES = (1, 2, 10, 50, 100, 1000)
FAMILIES = ("ConstantCp", "Nasa7")
PHASES = ("Static", "Dynamic")
CALCULATIONS = (
    "helmholtz",
    "pressure",
    "internal_energy",
    "enthalpy",
    "entropy",
    "gibbs",
    "dp_dc",
    "dp_dT",
    "cv",
    "cp",
    "sound_speed_sq",
    "chemical_potential",
    "log_fugacity_coeff",
    "fugacity",
)
BENCHMARK_NAME = re.compile(r"^(Static|Dynamic)/N([0-9]+)/(ConstantCp|Nasa7)/([a-zA-Z0-9_]+)$")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()

    completed = subprocess.run(
        [str(args.executable.resolve()), "--benchmark_list_tests=true"],
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise RuntimeError(f"benchmark list command failed ({completed.returncode}):\n{output}")
    if "Failed to match any benchmarks" in output:
        raise AssertionError(f"default sweep requested an unregistered benchmark:\n{output}")

    observed = Counter(line for line in output.splitlines() if BENCHMARK_NAME.fullmatch(line))
    expected = Counter(
        f"{phase}/N{size}/{family}/{calculation}"
        for size, calculation, family, phase in itertools.product(SIZES, CALCULATIONS, FAMILIES, PHASES)
    )
    if observed != expected:
        missing = sorted((expected - observed).elements())
        unexpected = sorted((observed - expected).elements())
        raise AssertionError(
            f"ideal benchmark registry mismatch: expected={sum(expected.values())}, "
            f"observed={sum(observed.values())}, missing={missing}, unexpected={unexpected}"
        )

    print(f"Ideal benchmark registry/default sweep: {sum(observed.values())} expected cases, zero unmatched filters.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
