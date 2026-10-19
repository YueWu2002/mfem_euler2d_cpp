#!/usr/bin/env python3
"""Run the requested k=2 DG/HDG SSP-DIRK3 convergence study."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor


def summary(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--grids", nargs="+", type=int, default=[32, 64, 128, 256, 512])
    parser.add_argument("--T", type=float, default=0.2)
    parser.add_argument("--cfl", type=float, default=2.0)
    parser.add_argument("--newton-tol", type=float, default=1e-13)
    parser.add_argument("--linear-tol", type=float, default=1e-13)
    parser.add_argument("--output", default="accuracy_test/convergence_order/convergence_k2")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Independent serial cases in parallel; concurrent timings are not benchmarks")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    if args.newton_tol <= 0 or args.linear_tol <= 0:
        parser.error('--newton-tol and --linear-tol must be positive')
    if sorted(args.grids) != args.grids or len(set(args.grids)) != len(args.grids):
        parser.error("--grids must be strictly increasing")

    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    exe = str(Path(args.exe).resolve())
    env = dict(os.environ)
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[name] = "1"

    def run_case(n):
        prefix = root / f"n{n}"
        command = [exe, "--nx", str(n), "--ny", str(n), "--k", "2",
                   "--T", str(args.T), "--cfl", str(args.cfl), "--form", "both",
                   "--integrator", "ssp-dirk3", "--assembly", "partial",
                   "--precond", "block-jacobi",
                   "--newton-tol", str(args.newton_tol),
                   "--linear-tol", str(args.linear_tol),
                   "--no-save-solution",
                   "--output", str(prefix)]
        # Keep diagnostics on disk even if Slurm interrupts the sweep.
        with (root / f"n{n}.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        if result.returncode:
            raise RuntimeError(f"Convergence run failed: {command}; see {root / f'n{n}.log'}")
        print(f'Completed DG+HDG N={n}', flush=True)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(run_case, args.grids))
    records = []
    components = ["rho", "rhou", "rhov", "energy"]
    previous = {}
    for n in args.grids:
        prefix = root / f"n{n}"
        comparison = summary(Path(f"{prefix}_comparison.csv"))
        for form in ["dg", "hdg"]:
            row = summary(Path(f"{prefix}_{form}_summary.csv"))
            errors = {component: float(row[f"{component}_L2"])
                      for component in components}
            rates = {component: "" for component in components}
            if form in previous:
                old_n, old_errors = previous[form]
                rates = {component:
                         math.log(old_errors[component]/errors[component])/math.log(n/old_n)
                         for component in components}
            previous[form] = (n, errors)
            # Default alpha=pi/4 gives hx=hy=2*sqrt(2)/n.
            records.append({"form": form, "n": n, "h": 2*math.sqrt(2)/n,
                            **{f"{c}_L2": errors[c] for c in components},
                            **{f"{c}_rate": rates[c] for c in components},
                            "cfl": row["cfl"],
                            "newton_tol": args.newton_tol,
                            "linear_tol": args.linear_tol,
                            "max_nodal_dg_hdg_difference": comparison["max_nodal_difference"],
                            "rho_L2_dg_hdg_difference": comparison["rho_L2_difference"],
                            **{f"{c}_L2_dg_hdg_difference": comparison[f"{c}_L2_difference"]
                               for c in ["rhou", "rhov", "energy"]},
                            "solution_L2_dg_hdg_difference": comparison["solution_L2_difference"],
                            "time_steps": row["time_steps"],
                            "newton_updates": row["newton_updates"],
                            "max_conservation_drift": row["max_conservation_drift"],
                            "mean_newton_s": row["mean_newton_s"],
                            "time_loop_s": row["time_loop_s"]})

    path = root / "convergence.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    for row in records:
        rate = "-" if row["rho_rate"] == "" else f"{row['rho_rate']:.4f}"
        print(f"{row['form']:3s} n={row['n']:4d} rho_L2={row['rho_L2']:.8e} "
              f"rate={rate} time_loop_s={float(row['time_loop_s']):.6g}")


if __name__ == "__main__":
    main()
