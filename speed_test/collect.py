#!/usr/bin/env python3
"""Collect completed speed cases into CPU time per Newton step.csv."""
import csv
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--update-results", action="store_true")
parser.add_argument("--root", type=Path,
                    default=Path(__file__).resolve().parent / "results" / "be_newton3")
args = parser.parse_args()
root = args.root.resolve()
root.mkdir(parents=True, exist_ok=True)
rows = []
for path in sorted(root.glob("job_*/k*_*/result.csv")):
    with path.open(newline="", encoding="utf-8") as stream:
        data = list(csv.DictReader(stream))
    if len(data) != 2:
        raise RuntimeError(f"Expected DG and HDG rows in {path}")
    by_form = {row["form"]: row for row in data}
    dg, hdg = by_form["dg"], by_form["hdg"]
    for row in data:
        if (row.get("integrator") != "be"
                or int(row["measured_newton_updates"]) != 3):
            raise RuntimeError(f"Expected backward Euler with 3 updates: {path}")
    rows.append({
        "k": int(dg["k"]),
        "preconditioner": dg["preconditioner"],
        "CFL": dg["cfl"],
        "DG CPU time per Newton step (s)": dg["mean_newton_s"],
        "HDG CPU time per Newton step (s)": hdg["mean_newton_s"],
        "DG/HDG": dg["dg_over_hdg"],
        "DG peak RSS (GiB)": dg.get("peak_rss_gib", ""),
        "HDG peak RSS (GiB)": hdg.get("peak_rss_gib", ""),
        "node": dg["host"],
        "job_id": dg["slurm_job_id"],
    })
rows.sort(key=lambda row: (row["k"], row["preconditioner"]))
output = root / "CPU time per Newton step.csv"
if rows:
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(output)
headers = list(rows[0]) if rows else ["k", "preconditioner", "CFL",
    "DG CPU time per Newton step (s)", "HDG CPU time per Newton step (s)",
    "DG/HDG", "DG peak RSS (GiB)", "HDG peak RSS (GiB)", "node", "job_id"]
lines = ["### CPU time per Newton step", "",
         "Configuration: 256x256, CFL=5.0, one backward-Euler time step,",
         "three fixed Newton updates (3 measured updates), matrix-free partial",
         "assembly, restarted weighted right GMRES (restart dimension 20, "
         "maximum 1500 Krylov iterations), and one process/thread.", "",
         "Peak RSS is measured per solver process over its complete lifetime,",
         "including setup; missing measurements are shown as —, not estimates.", "",
         "| " + " | ".join(headers) + " |",
         "| " + " | ".join(["---"] * len(headers)) + " |"]
for row in rows:
    lines.append("| " + " | ".join(str(row[h]) if row[h] != "" else "—" for h in headers) + " |")
if not rows:
    lines.append("| " + " | ".join(["—"] * len(headers)) + " |")
memory_rows = []
for path in sorted(root.glob("job_*/k*_*/resources.csv")):
    with path.open(newline="", encoding="utf-8") as stream:
        memory_rows.extend(csv.DictReader(stream))
memory_rows.sort(key=lambda row: (int(row["k"]), row["preconditioner"], row["form"]))
lines += ["", "### Measured process peak RSS", "",
          "Failed processes are retained here; their RSS is not a completed-run estimate.", "",
          "| k | preconditioner | form | status | peak RSS (GiB) | wall (s) | node | job ID |",
          "|---:|---|---|---|---:|---:|---|---:|"]
for row in memory_rows:
    lines.append("| {k} | {preconditioner} | {form} | {status} | {peak_rss_gib} | {wall_s} | {host} | {slurm_job_id} |".format(**row))
if not memory_rows:
    lines.append("| — | — | — | no measurements | — | — | — | — |")
markdown = "\n".join(lines) + "\n"
(root / "CPU time per Newton step.md").write_text(markdown)
if args.update_results:
    target = Path(__file__).resolve().parent.parent / "RESULTS.md"
    content = target.read_text()
    start, end = "<!-- speed-test-be-results:start -->", "<!-- speed-test-be-results:end -->"
    if content.count(start) != 1 or content.count(end) != 1:
        raise RuntimeError("Expected exactly one speed-test table marker pair")
    first = content.index(start) + len(start)
    last = content.index(end)
    if last < first:
        raise RuntimeError("Invalid speed-test table markers")
    target.write_text(content[:first] + "\n" + markdown + content[last:])
