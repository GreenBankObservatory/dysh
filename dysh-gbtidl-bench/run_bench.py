#!/usr/bin/env python
"""
dysh vs GBTIDL benchmark runner.

Usage:
    uv run python run_bench.py [--benchmarks NAME ...] --mode {warm,cold,both}
                               --iterations N [--tmpdir /path/to/tmp]
                               [--output results.json] [--verbose]

Cold cache mode copies data to a unique temp directory before each timed run
and calls posix_fadvise(DONTNEED) on every file to evict them from the OS page
cache.  For a guaranteed cold cache (root required):
    echo 3 | sudo tee /proc/sys/vm/drop_caches

Subtract the 'exit' baseline from other benchmarks to get net data-processing
time (the 'exit' benchmark measures Python/dysh process-startup overhead only).
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path
from statistics import mean, stdev

from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

console = Console()

SCRIPT_MS_RE = re.compile(r"(?P<tool>DYSH|GBTIDL)_BENCH_SCRIPT_MS=\s*(?P<ms>[0-9.dDeE+\-]+)")
STAGE_MS_RE = re.compile(r"(?P<tool>DYSH|GBTIDL)_BENCH_STAGE_MS\[(?P<stage>[^\]]+)\]=\s*(?P<ms>[0-9.dDeE+\-]+)")

# ---------------------------------------------------------------------------
# Benchmark registry
# ---------------------------------------------------------------------------

BENCHMARKS = {
    "argus_vanecal": {
        "dysh_script": "scripts/argus_vanecal/dysh_script.py",
        "gbtidl_script": "scripts/argus_vanecal/gbtidl.pro",
        "data_path": "/home/scratch/ajschmie/training/dysh/datasets/argus/TGBT22A_603_05_vanecal.raw.vegas",
        "has_output": False,
        "verify_script": "scripts/argus_vanecal/verify.py",
        "verify_needs_stdout": True,
    },
    "hi_otf": {
        "dysh_script": "scripts/hi_otf/dysh_script.py",
        "gbtidl_script": "scripts/hi_otf/gbtidl.pro",
        "data_path": "/home/scratch/dfrayer/DATAdemo/TGBT17A_506_11.raw.vegas",
        "has_output": True,
        "verify_script": "scripts/hi_otf/verify.py",
        "verify_needs_stdout": False,
    },
    "hi_survey": {
        "dysh_script": "scripts/hi_survey/dysh_script.py",
        "gbtidl_script": "scripts/hi_survey/gbtidl.pro",
        "data_path": "/home/astro-util/HIsurvey/Session02",
        "has_output": False,
        "verify_script": "scripts/hi_survey/verify.py",
        "verify_needs_stdout": True,
    },
    "nod_kfpa": {
        "dysh_script": "scripts/nod_kfpa/dysh_script.py",
        "gbtidl_script": None,
        "data_path": "/home/dysh/example_data/nod-KFPA/data/TGBT22A_503_02.raw.vegas",
        "has_output": False,
    },
    "summary_lib": {
        "dysh_script": "scripts/summary_lib/dysh_script.py",
        "gbtidl_script": "scripts/summary_lib/gbtidl.pro",
        "data_path": "/home/astro-util/HIsurvey/Session02",
        "has_output": False,
    },
    "summary_cli": {
        "dysh_script": "scripts/summary_cli/dysh_script.py",
        "gbtidl_script": "scripts/summary_cli/gbtidl.pro",
        "data_path": "/home/astro-util/HIsurvey/Session02",
        "has_output": False,
    },
    "exit": {
        "dysh_script": "scripts/exit/dysh_script.py",
        "gbtidl_script": "scripts/exit/gbtidl",
        "data_path": None,
        "has_output": False,
        "script_body_zero": True,
    },
}

# ---------------------------------------------------------------------------
# Page-cache eviction
# ---------------------------------------------------------------------------


def _evict_from_pagecache(path: Path, verbose: bool = False) -> None:
    """Walk all files under *path* and call posix_fadvise(DONTNEED) on each."""
    if not hasattr(os, "posix_fadvise"):
        console.print("[yellow]warn:[/] posix_fadvise not available on this platform — skipping eviction")
        return
    for fpath in sorted(path.rglob("*")):
        if not fpath.is_file():
            continue
        try:
            fd = os.open(str(fpath), os.O_RDONLY)
            try:
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                if verbose:
                    console.print(f"  evicted: {fpath}")
            finally:
                os.close(fd)
        except OSError as exc:
            console.print(f"[yellow]warn:[/] posix_fadvise failed for {fpath}: {exc}")


# ---------------------------------------------------------------------------
# Runner detection
# ---------------------------------------------------------------------------


def _gbtidl_available() -> bool:
    return shutil.which("gbtidl") is not None


def _make_cmd(tool: str, script: str) -> list[str]:
    if tool == "dysh":
        return ["uv", "run", "python", script]
    else:
        return ["gbtidl", script]


# ---------------------------------------------------------------------------
# Single-iteration helpers
# ---------------------------------------------------------------------------


def _parse_script_seconds(output: str) -> float | None:
    match = SCRIPT_MS_RE.search(output)
    if match is None:
        return None
    return float(match.group("ms").replace("D", "E").replace("d", "e")) / 1000.0


def _parse_stage_seconds(output: str) -> dict[str, float]:
    stages: dict[str, float] = {}
    for match in STAGE_MS_RE.finditer(output):
        stages[match.group("stage")] = float(match.group("ms").replace("D", "E").replace("d", "e")) / 1000.0
    return stages


def _read_peak_rss_kb(pid: int) -> int:
    """Read the current high-water RSS for a process from procfs."""
    vm_hwm_kb = 0
    vm_rss_kb = 0
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    vm_hwm_kb = int(line.split()[1])
                elif line.startswith("VmRSS:"):
                    vm_rss_kb = int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return 0
    return max(vm_hwm_kb, vm_rss_kb)


def _drain_stdout(stream, output_lines: list[str], verbose: bool) -> None:
    for line in stream:
        output_lines.append(line)
        if verbose:
            console.print(line, end="")


def _run_script(cmd: list[str], env: dict, verbose: bool = False, require_script_marker: bool = True) -> dict:
    """Run *cmd* with *env*. Exits on failure."""
    if verbose:
        console.print(f"[dim]$ {' '.join(cmd)}[/]")
    t0 = time.perf_counter()
    proc = subprocess.Popen(
        cmd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    output_lines: list[str] = []
    assert proc.stdout is not None
    reader = threading.Thread(target=_drain_stdout, args=(proc.stdout, output_lines, verbose), daemon=True)
    reader.start()
    peak_rss_kb = 0
    while proc.poll() is None:
        peak_rss_kb = max(peak_rss_kb, _read_peak_rss_kb(proc.pid))
        time.sleep(0.01)
    peak_rss_kb = max(peak_rss_kb, _read_peak_rss_kb(proc.pid))
    reader.join()
    result = subprocess.CompletedProcess(cmd, proc.returncode, "".join(output_lines), "")
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        console.print(f"[red bold]FAILED[/] {' '.join(cmd)} exited with code {result.returncode}")
        if not verbose and (result.stdout or result.stderr):
            console.print(result.stdout, end="")
            console.print(result.stderr, end="", style="red")
        sys.exit(1)
    script_body_s = _parse_script_seconds(result.stdout)
    if require_script_marker and script_body_s is None:
        console.print(f"[red bold]FAILED[/] {' '.join(cmd)} did not emit a BENCH_SCRIPT_MS marker")
        if not verbose and result.stdout:
            console.print(result.stdout, end="")
        sys.exit(1)
    return {
        "elapsed_s": elapsed,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "script_body_s": script_body_s,
        "stage_s": _parse_stage_seconds(result.stdout),
        "peak_rss_mb": peak_rss_kb / 1024.0,
    }


def _run_one_cold(
    tool: str,
    script: str,
    data_path: str | None,
    has_output: bool,
    tmpdir: str | None,
    verbose: bool,
    require_script_marker: bool,
) -> float:
    """Copy data to a fresh temp dir, evict from page cache, time the run."""
    tmp_root = Path(tempfile.mkdtemp(prefix=f"dysh_cold_{uuid.uuid4().hex[:8]}_", dir=tmpdir))
    try:
        env = dict(os.environ)
        if data_path:
            tmp_data = tmp_root / Path(data_path).name
            if verbose:
                console.print(f"  copying {data_path} -> {tmp_data}")
            shutil.copytree(data_path, tmp_data)
            _evict_from_pagecache(tmp_data, verbose=verbose)
            env["DYSH_BENCH_DATA_PATH"] = str(tmp_data)
        if has_output:
            env["DYSH_BENCH_OUT_PATH"] = str(tmp_root / "output.fits")
        result = _run_script(_make_cmd(tool, script), env, verbose=verbose, require_script_marker=require_script_marker)
        return result
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)


def _run_iterations(
    tool: str,
    label: str,
    script: str,
    data_path: str | None,
    has_output: bool,
    cache_mode: str,
    n_iterations: int,
    tmpdir: str | None,
    verbose: bool,
    overall_task: int | None,
    progress: Progress | None,
    require_script_marker: bool,
) -> list[dict]:
    env = dict(os.environ)
    if data_path:
        env["DYSH_BENCH_DATA_PATH"] = data_path
    if has_output:
        env["DYSH_BENCH_OUT_PATH"] = "bench_output.fits"

    cmd = _make_cmd(tool, script)

    def _run():
        return _run_script(cmd, env, verbose=verbose, require_script_marker=require_script_marker)

    if cache_mode == "warm":
        if verbose:
            console.print(f"[dim]--- {label} warmup ---[/]")
        elif overall_task is not None:
            progress.update(overall_task, description=f"{label} — warmup")
        _run()
        if overall_task is not None:
            progress.advance(overall_task)

    runs: list[dict] = []
    for i in range(n_iterations):
        if verbose:
            console.print(f"[dim]--- {label} trial {i + 1}/{n_iterations} ---[/]")
        elif overall_task is not None:
            progress.update(overall_task, description=f"{label} — trial")
        if cache_mode == "cold":
            run = _run_one_cold(tool, script, data_path, has_output, tmpdir, verbose, require_script_marker)
        else:
            run = _run()
        runs.append(run)
        if overall_task is not None:
            progress.advance(overall_task)

    return runs


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------


def _run_verify(name: str, cfg: dict, data_path: str | None, has_gbtidl: bool, verbose: bool) -> None:
    """Run the benchmark's verify.py once, capturing stdout as needed."""
    verify_script = cfg.get("verify_script")
    if not verify_script or not Path(verify_script).exists():
        console.print(f"[yellow]verify:[/] no verify.py for [bold]{name}[/], skipping")
        return

    needs_stdout = cfg.get("verify_needs_stdout", False)
    env = dict(os.environ)
    if data_path:
        env["DYSH_BENCH_DATA_PATH"] = data_path

    if needs_stdout:
        if not has_gbtidl:
            console.print(f"[yellow]verify:[/] gbtidl not available — cannot verify [bold]{name}[/]")
            return

        dysh_out = Path(tempfile.mktemp(suffix=f".{name}.dysh.out"))
        gbtidl_out = Path(tempfile.mktemp(suffix=f".{name}.gbtidl.out"))
        try:
            for tool, script, outfile in [
                ("dysh", cfg["dysh_script"], dysh_out),
                ("gbtidl", cfg["gbtidl_script"], gbtidl_out),
            ]:
                if verbose:
                    console.print(f"[dim]verify: capturing {tool} stdout for {name}[/]")
                result = subprocess.run(_make_cmd(tool, script), env=env, capture_output=True, text=True)
                if result.returncode != 0:
                    console.print(f"[red bold]VERIFY FAIL[/] {name} ({tool} exited with code {result.returncode})")
                    if result.stdout:
                        console.print(result.stdout, end="")
                    if result.stderr:
                        console.print(result.stderr, end="", style="red")
                    sys.exit(1)
                outfile.write_text(result.stdout)

            result = subprocess.run(
                ["uv", "run", "python", verify_script, str(dysh_out), str(gbtidl_out)],
                capture_output=False,
                text=True,
            )
            if result.returncode == 0:
                console.print(f"[green bold]VERIFY PASS[/] {name}")
            else:
                console.print(f"[red bold]VERIFY FAIL[/] {name}")
                sys.exit(1)
        finally:
            dysh_out.unlink(missing_ok=True)
            gbtidl_out.unlink(missing_ok=True)
    else:
        # hi_otf: verify.py reads FITS files from its own directory
        script_dir = Path(verify_script).parent
        result = subprocess.run(
            ["uv", "run", "python", Path(verify_script).name],
            cwd=script_dir,
            capture_output=False,
            text=True,
        )
        if result.returncode == 0:
            console.print(f"[green bold]VERIFY PASS[/] {name}")
        else:
            console.print(f"[red bold]VERIFY FAIL[/] {name}")
            sys.exit(1)


# ---------------------------------------------------------------------------
# Stats + table
# ---------------------------------------------------------------------------


def _stats(times: list[float]) -> dict:
    n = len(times)
    m = mean(times)
    s = stdev(times) if n > 1 else 0.0
    return {"n": n, "mean": m, "std": s, "min": min(times), "max": max(times), "times": times}


def _stats_from_runs(runs: list[dict]) -> dict:
    stats = _stats([run["elapsed_s"] for run in runs])
    peak_rss_values = [run["peak_rss_mb"] for run in runs if run.get("peak_rss_mb") is not None]
    if peak_rss_values:
        stats["peak_rss_mb"] = _stats(peak_rss_values)
    script_times = [run["script_body_s"] for run in runs if run["script_body_s"] is not None]
    if len(script_times) == len(runs):
        stats["script_body"] = _stats(script_times)
        stats["startup_overhead"] = _stats([max(0.0, run["elapsed_s"] - run["script_body_s"]) for run in runs])
    all_stage_keys = sorted({key for run in runs for key in run.get("stage_s", {}).keys()})
    if all_stage_keys:
        stage_stats = {}
        for key in all_stage_keys:
            stage_times = [run["stage_s"][key] for run in runs if key in run.get("stage_s", {})]
            if len(stage_times) == len(runs):
                stage_stats[key] = _stats(stage_times)
        if stage_stats:
            stats["stages"] = stage_stats
    return stats


def _apply_zero_script_body(stats: dict) -> dict:
    zeroes = [0.0] * stats["n"]
    stats["script_body"] = _stats(zeroes)
    stats["startup_overhead"] = _stats(list(stats["times"]))
    return stats


def _print_results(results: dict, modes: list[str], all_columns: bool = False) -> None:
    table = Table(show_header=True, header_style="bold")
    table.add_column("Benchmark")
    table.add_column("Mode")
    table.add_column("Tool")
    table.add_column("N", justify="right")
    table.add_column("Total (s)", justify="right")
    table.add_column("Script (s)", justify="right")
    if all_columns:
        table.add_column("Startup (s)", justify="right")
        table.add_column("Peak RSS", justify="right")
        table.add_column("Std (s)", justify="right")
    table.add_column("Total Speedup", justify="right")
    table.add_column("Script x", justify="right")
    if all_columns:
        table.add_column("Startup x", justify="right")

    for name, mode_data in results.items():
        for mode in modes:
            if mode not in mode_data:
                continue
            tool_data = mode_data[mode]
            dysh_proc_mean = tool_data.get("dysh", {}).get("mean")
            gbtidl_proc_mean = tool_data.get("gbtidl", {}).get("mean")
            dysh_script_mean = tool_data.get("dysh", {}).get("script_body", {}).get("mean")
            gbtidl_script_mean = tool_data.get("gbtidl", {}).get("script_body", {}).get("mean")
            dysh_start_mean = tool_data.get("dysh", {}).get("startup_overhead", {}).get("mean")
            gbtidl_start_mean = tool_data.get("gbtidl", {}).get("startup_overhead", {}).get("mean")

            for tool in ("dysh", "gbtidl"):
                if tool not in tool_data:
                    continue
                st = tool_data[tool]
                script_speedup = ""
                startup_speedup = ""
                total_speedup = ""
                if tool == "dysh":
                    if dysh_proc_mean and gbtidl_proc_mean:
                        total_speedup = f"{gbtidl_proc_mean / dysh_proc_mean:.2f}x"
                    if dysh_script_mean and gbtidl_script_mean:
                        script_speedup = f"{gbtidl_script_mean / dysh_script_mean:.2f}x"
                    if dysh_start_mean and gbtidl_start_mean:
                        startup_speedup = f"{gbtidl_start_mean / dysh_start_mean:.2f}x"
                elif tool == "gbtidl" and (gbtidl_proc_mean or gbtidl_script_mean or gbtidl_start_mean):
                    total_speedup = "[dim](ref)[/]" if dysh_proc_mean and gbtidl_proc_mean else ""
                    script_speedup = "[dim](ref)[/]" if dysh_script_mean and gbtidl_script_mean else ""
                    startup_speedup = "[dim](ref)[/]" if dysh_start_mean and gbtidl_start_mean else ""
                table.add_row(
                    name,
                    mode,
                    tool,
                    str(st["n"]),
                    f"{st['mean']:.2f}",
                    f"{st['script_body']['mean']:.2f}" if "script_body" in st else "",
                    *(
                        [
                            f"{st['startup_overhead']['mean']:.2f}" if "startup_overhead" in st else "",
                            f"{st['peak_rss_mb']['max']:.0f} MB" if "peak_rss_mb" in st else "",
                            f"{st['std']:.2f}",
                        ]
                        if all_columns
                        else []
                    ),
                    total_speedup,
                    script_speedup,
                    *([startup_speedup] if all_columns else []),
                )

    console.print(table)

    stage_tables_printed = 0
    for name, mode_data in results.items():
        for mode in modes:
            if mode not in mode_data:
                continue
            tool_data = mode_data[mode]
            dysh = tool_data.get("dysh", {})
            gbtidl = tool_data.get("gbtidl", {})
            dysh_stages = dysh.get("stages", {})
            gbtidl_stages = gbtidl.get("stages", {})
            common = [stage for stage in dysh_stages if stage in gbtidl_stages]

            def _print_single_tool_stage_table(tool: str, stages: dict, script_mean: float | None) -> None:
                nonlocal stage_tables_printed
                if not stages:
                    return
                stage_tables_printed += 1
                stage_table = Table(
                    title=f"{name} {mode} {tool} stage breakdown",
                    show_header=True,
                    header_style="bold",
                )
                stage_table.add_column("Stage")
                stage_table.add_column(f"{tool} (s)", justify="right")
                stage_table.add_column(f"{tool} %", justify="right")
                for stage, stage_stats in stages.items():
                    value = stage_stats["mean"]
                    share = f"{(100.0 * value / script_mean):.1f}%" if script_mean else ""
                    stage_table.add_row(stage, f"{value:.3f}", share)
                console.print()
                console.print(stage_table)

            if common:
                stage_tables_printed += 1
                stage_table = Table(
                    title=f"{name} {mode} stage breakdown",
                    show_header=True,
                    header_style="bold",
                )
                stage_table.add_column("Stage")
                stage_table.add_column("dysh (s)", justify="right")
                stage_table.add_column("gbtidl (s)", justify="right")
                stage_table.add_column("dysh/gbtidl", justify="right")
                stage_table.add_column("dysh %", justify="right")
                dysh_script_mean = dysh.get("script_body", {}).get("mean")
                for stage in common:
                    dy = dysh_stages[stage]["mean"]
                    gb = gbtidl_stages[stage]["mean"]
                    ratio = f"{dy / gb:.2f}x" if gb else ""
                    share = f"{(100.0 * dy / dysh_script_mean):.1f}%" if dysh_script_mean else ""
                    stage_table.add_row(stage, f"{dy:.3f}", f"{gb:.3f}", ratio, share)
                console.print()
                console.print(stage_table)
                continue

            _print_single_tool_stage_table("dysh", dysh_stages, dysh.get("script_body", {}).get("mean"))
            _print_single_tool_stage_table("gbtidl", gbtidl_stages, gbtidl.get("script_body", {}).get("mean"))

    if stage_tables_printed:
        console.print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--benchmarks",
        nargs="+",
        default=list(BENCHMARKS),
        choices=list(BENCHMARKS),
        metavar="NAME",
        help="benchmarks to run (default: all)",
    )
    parser.add_argument(
        "--mode",
        choices=["warm", "cold", "both"],
        default="warm",
        help="cache mode (default: warm)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=3,
        metavar="N",
        help="number of timed iterations per benchmark/mode (default: 3)",
    )
    parser.add_argument(
        "--tmpdir",
        default=None,
        metavar="DIR",
        help="parent directory for cold-mode temp copies (default: system temp)",
    )
    parser.add_argument(
        "--output",
        default=None,
        metavar="FILE",
        help="write JSON results to this file",
    )
    parser.add_argument(
        "--dysh-only",
        action="store_true",
        help="skip gbtidl runs even if gbtidl is on PATH",
    )
    parser.add_argument(
        "--data-path",
        default=None,
        metavar="PATH",
        help="override data path for all selected benchmarks",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="stream script stdout/stderr live and show per-file eviction",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="after timing, run each benchmark's verify.py to check dysh/gbtidl agreement",
    )
    parser.add_argument(
        "--all-columns",
        action="store_true",
        help="show startup, RSS, stddev, and startup speedup columns in the summary table",
    )
    args = parser.parse_args()

    modes = ["warm", "cold"] if args.mode == "both" else [args.mode]
    has_gbtidl = _gbtidl_available() and not args.dysh_only
    if args.dysh_only:
        console.print("[dim]--dysh-only: skipping gbtidl runs[/]\n")
    elif not has_gbtidl:
        console.print("[yellow]warn:[/] gbtidl not found on PATH — skipping GBTIDL runs\n")

    results: dict = {}

    # Compute total steps: each tool run = n_iterations timed + 1 warmup (warm mode only)
    steps_per_run = args.iterations + (1 if args.mode != "cold" else 0)
    n_total_steps = sum(
        steps_per_run * (1 + (1 if has_gbtidl and BENCHMARKS[n]["gbtidl_script"] else 0)) * len(modes)
        for n in args.benchmarks
    )

    def _run_all(overall, progress):
        for name in args.benchmarks:
            cfg = BENCHMARKS[name]
            data_path = args.data_path or cfg["data_path"]
            has_output = cfg.get("has_output", False)
            results[name] = {}

            for mode in modes:
                results[name][mode] = {}
                require_script_marker = not cfg.get("script_body_zero", False)

                runs = _run_iterations(
                    "dysh",
                    f"dysh ({name})",
                    cfg["dysh_script"],
                    data_path,
                    has_output,
                    mode,
                    args.iterations,
                    args.tmpdir,
                    args.verbose,
                    overall,
                    progress,
                    require_script_marker,
                )
                dysh_stats = _stats_from_runs(runs)
                if cfg.get("script_body_zero"):
                    dysh_stats = _apply_zero_script_body(dysh_stats)
                results[name][mode]["dysh"] = dysh_stats

                if has_gbtidl and cfg["gbtidl_script"]:
                    runs = _run_iterations(
                        "gbtidl",
                        f"gbtidl ({name})",
                        cfg["gbtidl_script"],
                        data_path,
                        has_output,
                        mode,
                        args.iterations,
                        args.tmpdir,
                        args.verbose,
                        overall,
                        progress,
                        require_script_marker,
                    )
                    gbtidl_stats = _stats_from_runs(runs)
                    if cfg.get("script_body_zero"):
                        gbtidl_stats = _apply_zero_script_body(gbtidl_stats)
                    results[name][mode]["gbtidl"] = gbtidl_stats

            if args.verify:
                _run_verify(name, cfg, data_path, has_gbtidl, args.verbose)

    if args.verbose:
        _run_all(None, None)
    else:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            MofNCompleteColumn(),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
            transient=True,
        ) as progress:
            overall = progress.add_task("starting", total=n_total_steps)
            _run_all(overall, progress)

    _print_results(results, modes, all_columns=args.all_columns)

    if args.output:
        with open(args.output, "w") as fh:
            json.dump(results, fh, indent=2)
        console.print(f"\nResults written to [bold]{args.output}[/]")


if __name__ == "__main__":
    main()
