#!/usr/bin/env python3
"""Single-thread paired timings: warm up, alternate DG/HDG, report medians.

Measures each run's mean Newton time; does not substitute subprocess wall time.
The executable's iteration timer includes residual/Jacobian/condensation/assembly,
preconditioner factorization, linear solve, recovery, validation and line search.
"""
import argparse
import csv
import os
from pathlib import Path
import statistics
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", required=True)
    p.add_argument("--nx", type=int, default=16)
    p.add_argument("--ny", type=int, default=16)
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--cfl", type=float, default=1.)
    p.add_argument("--steps", type=int, default=3)
    p.add_argument("--newton-steps", type=int, default=3)
    p.add_argument("--repeat", type=int, default=5)
    p.add_argument("--precond", choices=["block-jacobi", "block-ilu0"], default="block-jacobi")
    p.add_argument("--assembly", choices=["partial", "element", "full"], default="partial")
    p.add_argument("--output", default="results/benchmark")
    args = p.parse_args()
    if args.repeat < 1:
        p.error("--repeat must be positive")
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    exe = str(Path(args.exe).resolve())
    env = dict(os.environ)
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[name] = "1"
    values = {"dg": [], "hdg": []}
    samples = []
    for repeat in range(-1, args.repeat):
        forms = ["dg", "hdg"] if repeat % 2 == 0 else ["hdg", "dg"]
        for form in forms:
            prefix = root / f"r{repeat}_{form}"
            command = [exe, "--form", form, "--nx", str(args.nx), "--ny", str(args.ny),
                "--k", str(args.k), "--cfl", str(args.cfl), "--steps", str(args.steps),
                "--newton-steps", str(args.newton_steps), "--precond", args.precond,
                "--assembly", args.assembly, "--timing-only", "--output", str(prefix)]
            result = subprocess.run(command, text=True, capture_output=True, env=env)
            (root / f"r{repeat}_{form}.log").write_text(result.stdout + result.stderr, encoding="utf-8")
            if result.returncode:
                raise RuntimeError(f"Benchmark failed: {command}\n{result.stderr}")
            with Path(f"{prefix}_{form}_summary.csv").open(newline="") as stream:
                row = next(csv.DictReader(stream))
            expected = args.steps*args.newton_steps
            if int(row["newton_updates"]) != expected:
                raise AssertionError("Unexpected Newton update count")
            if repeat >= 0:
                value = float(row["mean_newton_s"])
                values[form].append(value)
                samples.append({"repeat": repeat, "form": form, "mean_newton_s": value,
                                "newton_updates": expected})
    with (root / "samples.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=samples[0].keys())
        writer.writeheader(); writer.writerows(samples)
    dg, hdg = (statistics.median(values[f]) for f in ["dg", "hdg"])
    print(f"median mean Newton seconds: DG={dg:.8g}, HDG={hdg:.8g}; DG/HDG={dg/hdg:.5g}")


if __name__ == "__main__":
    main()
