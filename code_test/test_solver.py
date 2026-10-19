#!/usr/bin/env python3
"""Run the two intentionally small code tests.

The suite has exactly two responsibilities: a DG/HDG smoke test and a
backward-Euler, one-step, converged-Newton update-equivalence test. It has no
fixed short time step: both cases select their step from the requested CFL.
"""
import argparse
import csv
import math
from pathlib import Path
import subprocess


def row(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--output", default="code_test/results")
    args = parser.parse_args()

    exe = str(Path(args.exe).resolve())
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=True)

    def run(name, flags):
        prefix = root / name
        result = subprocess.run([exe, "--output", str(prefix),
                                 "--newton-tol", "1e-12", "--linear-tol", "1e-12",
                                 *map(str, flags)],
                                text=True, capture_output=True)
        (root / f"{name}.log").write_text(result.stdout + result.stderr,
                                            encoding="utf-8")
        if result.returncode:
            raise AssertionError(f"{name} failed:\n{result.stdout}\n{result.stderr}")
        return prefix

    smoke = run("smoke", [
        "--nx", 8, "--ny", 8, "--k", 1, "--cfl", .1, "--steps", 1,
        "--form", "both", "--no-save-solution"])
    for form in ("dg", "hdg"):
        summary = row(f"{smoke}_{form}_summary.csv")
        assert summary["integrator"] == "ssp-dirk3"
        assert int(summary["time_steps"]) == 1
        assert int(summary["newton_updates"]) > 0
    smoke_comparison = row(f"{smoke}_comparison.csv")
    assert float(smoke_comparison["max_nodal_difference"]) < 1e-8

    backward_euler = run("be_converged_update", [
        "--nx", 8, "--ny", 8, "--k", 2, "--cfl", .1, "--steps", 1,
        "--integrator", "be", "--form", "both", "--check-update",
        "--newton-tol", 1e-12, "--linear-tol", 1e-12,
        "--no-save-solution"])
    for form in ("dg", "hdg"):
        summary = row(f"{backward_euler}_{form}_summary.csv")
        assert summary["integrator"] == "be"
        assert int(summary["time_steps"]) == 1
        assert int(summary["fixed_newton_steps"]) == -1
        assert int(summary["newton_updates"]) > 0
    for update in rows(f"{backward_euler}_hdg_iterations.csv"):
        check = float(update["update_check_relative"])
        assert math.isfinite(check) and check < 1e-8
    comparison = row(f"{backward_euler}_comparison.csv")
    assert float(comparison["max_nodal_difference"]) < 1e-8
    assert float(comparison["solution_L2_difference"]) < 1e-8
    print("Smoke and backward-Euler converged-update checks passed.")


if __name__ == "__main__":
    main()
