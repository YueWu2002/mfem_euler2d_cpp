#!/usr/bin/env python3
"""Explicitly invoked numerical/CLI regression suite. Python standard library only.

Nothing here runs during configure/build. Use --quick to reduce the degree sweep;
--full-time additionally exercises converged DG/HDG evolution to T=0.2.
"""
import argparse
import csv
import math
from pathlib import Path
import subprocess


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--output", default="results/tests")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--skip-algebra", action="store_true")
    parser.add_argument("--full-time", action="store_true")
    args = parser.parse_args()
    exe = str(Path(args.exe).resolve())
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=True)

    def run(name, flags, success=True):
        prefix = root / name
        command = [exe, "--output", str(prefix), *map(str, flags)]
        result = subprocess.run(command, text=True, capture_output=True)
        (root / (name + ".log")).write_text(
            result.stdout + result.stderr, encoding="utf-8")
        if (result.returncode == 0) != success:
            raise AssertionError(f"{name}: unexpected exit {result.returncode}\n"
                                 f"{result.stdout}\n{result.stderr}")
        return prefix

    if not args.skip_algebra:
        prefix = run("algebra", ["--verify"])
        algebra = rows(str(prefix) + "_verification.csv")
        assert {int(r["k"]) for r in algebra} == set(range(9))
        invariant = rows(str(prefix) + "_weighted_gmres.csv")
        assert len(invariant) == 2
        for row in invariant:
            assert float(row["physical_relative_difference"]) < 1e-8
            assert row["original_iterations"] == row["transformed_iterations"]

    degrees = [0, 1, 2, 8] if args.quick else range(9)
    # Same nonlinear target and update count, differing global actions and PCs.
    for degree in degrees:
        reference = None
        for pc in ["block-jacobi", "block-ilu0"]:
            for assembly in ["partial", "element", "full"]:
                prefix = run(f"k{degree}_{pc}_{assembly}", [
                    "--nx", 3, "--ny", 2, "--k", degree,
                    "--steps", 1, "--dt", 1e-4, "--newton-steps", 2,
                    "--form", "both", "--assembly", assembly, "--precond", pc,
                    "--check-update", "--linear-tol", 1e-12])
                compare = rows(str(prefix) + "_comparison.csv")[0]
                assert float(compare["max_nodal_difference"]) < 1e-8
                current = rows(str(prefix) + "_dg_coefficients.csv")
                values = [float(r[c]) for r in current
                          for c in ["rho", "rhou", "rhov", "energy"]]
                if reference is None:
                    reference = values
                assert max(abs(a-b) for a, b in zip(reference, values)) < 1e-8
                for form in ["dg", "hdg"]:
                    summary = rows(f"{prefix}_{form}_summary.csv")[0]
                    updates = rows(f"{prefix}_{form}_iterations.csv")
                    assert int(summary["newton_updates"]) == len(updates) == 2
                    assert int(summary["time_steps"]) == 1
                    assert math.isclose(float(summary["mean_newton_s"]),
                                        sum(float(r["iteration_s"]) for r in updates)/2,
                                        rel_tol=1e-12)
                    assert float(summary["max_conservation_drift"]) < 1e-8
                    assert all(math.isfinite(float(r["linear_relative_residual"]))
                               for r in updates)

    # Centered geometry with an independently adjustable angle, and selected CFL.
    prefix = run("cfl_angle", ["--nx", 4, "--ny", 3, "--k", 2,
        "--alpha-deg", 30, "--cfl", .05, "--steps", 2,
        "--newton-steps", 2, "--form", "both"])
    for form in ["dg", "hdg"]:
        row = rows(f"{prefix}_{form}_summary.csv")[0]
        assert math.isclose(float(row["actual_max_cfl"]), .05, rel_tol=1e-13)
        assert int(row["newton_updates"]) == 4
    text = Path(str(prefix) + "_dg.mesh").read_text(encoding="utf-8").splitlines()
    index = text.index("vertices")
    count = int(text[index+1])
    assert int(text[index+2]) == 2
    vertices = [list(map(float, s.split())) for s in text[index+3:index+3+count]]
    for axis, half_length in [(0, 1/math.cos(math.pi/6)), (1, 2.)]:
        assert abs(min(v[axis] for v in vertices) + half_length) < 2e-5
        assert abs(max(v[axis] for v in vertices) - half_length) < 2e-5

    # Both optional direct references solve only the Schur system on HDG.
    run("dense", ["--nx", 3, "--ny", 2, "--k", 1, "--steps", 1,
                  "--dt", 1e-4, "--newton-steps", 1,
                  "--solver", "dense", "--form", "both", "--check-update"])
    for name, flags in [
        ("degree", ["--k", 9]), ("angle", ["--alpha-deg", 0]),
        ("cfl", ["--cfl", 0]), ("time", ["--T", 1]),
        ("quadrature", ["--k", 2, "--volume-order", 9]),
        ("direct_mf", ["--solver", "dense", "--matrix-free"]),
        ("preconditioner", ["--precond", "unsupported"]),
    ]:
        run("invalid_" + name, flags, success=False)
    if args.full_time:
        prefix = run("T02", ["--nx", 8, "--ny", 6, "--k", 2,
            "--cfl", .2, "--T", .2, "--form", "both", "--no-save-solution"])
        for form in ["dg", "hdg"]:
            assert math.isclose(float(rows(f"{prefix}_{form}_summary.csv")[0]["T"]), .2)
    print("All explicitly requested regression checks passed.")


if __name__ == "__main__":
    main()
