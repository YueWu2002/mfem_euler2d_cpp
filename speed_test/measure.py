#!/usr/bin/env python3
"""Measure each solver process separately; do not reuse cumulative child rusage."""
import argparse
import json
from pathlib import Path
import subprocess


def run_measured(command, output, env=None):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = Path(str(output) + ".time")
    result = subprocess.run(
        ["/usr/bin/time", "-o", str(raw), "-f",
         "peak_rss_kib=%M\nuser_cpu_s=%U\nsystem_cpu_s=%S\nwall_s=%e",
         "--", *command], text=True, capture_output=True, env=env)
    values = {}
    for line in raw.read_text().splitlines():
        key, separator, value = line.partition("=")
        if separator and key in {"peak_rss_kib", "user_cpu_s", "system_cpu_s", "wall_s"}:
            values[key] = float(value)
    if "peak_rss_kib" not in values:
        raise RuntimeError(f"Missing peak RSS in {raw}")
    values["peak_rss_gib"] = values["peak_rss_kib"] / 1024**2
    values["returncode"] = result.returncode
    values["status"] = "completed" if result.returncode == 0 else "failed"
    values["command"] = command
    output.write_text(json.dumps(values, indent=2) + "\n")
    return result, values


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("A command is required")
    result, values = run_measured(command, args.output)
    print(result.stdout, end="")
    print(result.stderr, end="", file=__import__("sys").stderr)
    print(json.dumps(values, indent=2))
    raise SystemExit(result.returncode)
