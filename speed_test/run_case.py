#!/usr/bin/env python3
"""Run one (degree, preconditioner) DG/HDG Newton timing case."""
import argparse
import csv
import os
from pathlib import Path
import platform
import subprocess
from measure import run_measured


def summary(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--k", required=True, type=int)
    parser.add_argument("--precond", required=True,
                        choices=["block-jacobi", "block-ilu0"])
    parser.add_argument("--cfl", type=float, default=5.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if not 0 <= args.k <= 8:
        parser.error("--k must be in [0,8]")

    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    exe = str(Path(args.exe).resolve())
    env = dict(os.environ)
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[name] = "1"

    cfl = args.cfl
    rows = []
    resource_rows = []
    for form in ["dg", "hdg"]:
        prefix = output / f"cfl{cfl:g}_{form}"
        command = [exe, "--nx", "256", "--ny", "256", "--k", str(args.k),
                   "--cfl", str(cfl), "--steps", "1", "--form", form,
                   "--integrator", "ssp-dirk3", "--newton-steps", "3",
                   "--solver", "gmres", "--assembly", "partial",
                   "--precond", args.precond, "--timing-only",
                   "--output", str(prefix)]
        result, memory = run_measured(command, output / f"cfl{cfl:g}_{form}_resources.json", env)
        (output / f"cfl{cfl:g}_{form}.log").write_text(
            result.stdout + result.stderr, encoding="utf-8")
        resource_rows.append({
            "k": args.k, "preconditioner": args.precond, "form": form,
            "nx": 256, "ny": 256, "cfl": cfl,
            "status": memory["status"], "returncode": memory["returncode"],
            "peak_rss_kib": memory["peak_rss_kib"],
            "peak_rss_gib": memory["peak_rss_gib"],
            "process_user_cpu_s": memory["user_cpu_s"],
            "process_system_cpu_s": memory["system_cpu_s"],
            "wall_s": memory["wall_s"], "host": platform.node(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        })
        with (output / "resources.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=resource_rows[0].keys())
            writer.writeheader()
            writer.writerows(resource_rows)
        if result.returncode:
            raise RuntimeError(
                f"{form} failed at k={args.k}, CFL={cfl}: {result.stderr}")
        data = summary(Path(f"{prefix}_{form}_summary.csv"))
        expected = 4 * 3  # four DIRK stages, three fixed updates per stage
        if int(data["newton_updates"]) != expected:
            raise RuntimeError("Unexpected Newton update count")
        rows.append({
            "k": args.k,
            "preconditioner": args.precond,
            "form": form,
            "nx": 256,
            "ny": 256,
            "cfl": cfl,
            "dirk_stages": 4,
            "newton_updates_per_stage": 3,
            "measured_newton_updates": expected,
            "mean_newton_s": data["mean_newton_s"],
            "peak_rss_kib": memory["peak_rss_kib"],
            "peak_rss_gib": memory["peak_rss_gib"],
            "process_user_cpu_s": memory["user_cpu_s"],
            "process_system_cpu_s": memory["system_cpu_s"],
            "time_loop_s": data["time_loop_s"],
            "volume_dofs": data["volume_dofs"],
            "trace_dofs": data["trace_dofs"],
            "host": platform.node(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        })
    dg = float(rows[0]["mean_newton_s"])
    hdg = float(rows[1]["mean_newton_s"])
    ratio = dg/hdg
    for row in rows:
        row["dg_over_hdg"] = ratio
    path = output / "result.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"k={args.k} pc={args.precond} CFL={cfl:g} "
          f"DG={dg:.9g}s HDG={hdg:.9g}s DG/HDG={ratio:.6g}")
    return


if __name__ == "__main__":
    main()
