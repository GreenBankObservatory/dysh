#!/usr/bin/env python
"""
dysh vs GBTIDL workflow benchmark runner.

Usage (from any directory; script, golden, and verifier paths are relative to this file, and
every child process runs with this directory as its working directory):
    uv run python run_bench.py [--benchmarks NAME ...] --mode {warm,cold,both}
                               --iterations N [--tmpdir /path/to/tmp]
                               [--output results.json] [--verbose]

Data paths resolve in this order: --data-path, $DYSH_BENCH_DATA_PATH, the
canonical GBO path (if it exists on this host), then the dysh_data alias
(which may download the data on first use).

Cold cache mode copies data to a unique temp directory before each timed run
and calls posix_fadvise(DONTNEED) on every file to evict them from the OS page
cache.  For a guaranteed cold cache (root required):
    echo 3 | sudo tee /proc/sys/vm/drop_caches

Subtract the 'exit' baseline from other benchmarks to get net data-processing
time (the 'exit' benchmark measures Python/dysh process-startup overhead only).

Verification: --verify compares dysh output against GBTIDL when gbtidl is on
PATH, otherwise against a committed golden capture (scripts/<name>/golden.txt,
created with --capture-golden on the canonical dataset at a known-good
baseline). Golden values are dataset-dependent: only compare against them when
running on the canonical dataset.
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

HERE = Path(__file__).resolve().parent
"""Directory containing this file: base for relative registry paths and cwd of child processes."""

Script = str | list[str]
"""A benchmark script: a path, or an argv list whose first item is the path and the rest are arguments."""

SCRIPT_MS_RE = re.compile(r"(?P<tool>DYSH|GBTIDL)_BENCH_SCRIPT_MS=\s*(?P<ms>[0-9.dDeE+\-]+)")
STAGE_MS_RE = re.compile(r"(?P<tool>DYSH|GBTIDL)_BENCH_STAGE_MS\[(?P<stage>[^\]]+)\]=\s*(?P<ms>[0-9.dDeE+\-]+)")

# ---------------------------------------------------------------------------
# Benchmark registry
# ---------------------------------------------------------------------------

# ``dysh_script`` / ``gbtidl_script`` are a `Script`: a path, or an argv list such as
# ``["../bench_getps.py", "-t", "-l", "4"]``. Relative paths (the first item only) are resolved
# against HERE, so a script need not live under ``scripts/``. ``verify_script`` is a path.
BENCHMARKS = {
    "argus_vanecal": {
        "dysh_script": "scripts/argus_vanecal/dysh_script.py",
        "gbtidl_script": "scripts/argus_vanecal/gbtidl.pro",
        # Canonical: vane/sky subset of TGBT22A_603_05 (GBO only). The otf4
        # alias is the full session (same vane scans, bigger load stage).
        "data_path": "/home/scratch/ajschmie/training/dysh/datasets/argus/TGBT22A_603_05_vanecal.raw.vegas",
        "data_alias": {"example": "otf4"},
        "has_output": False,
        "verify_script": "scripts/argus_vanecal/verify.py",
        "verify_needs_stdout": True,
    },
    "hi_survey": {
        "dysh_script": "scripts/hi_survey/dysh_script.py",
        "gbtidl_script": "scripts/hi_survey/gbtidl.pro",
        "data_path": "/home/astro-util/HIsurvey/Session02",
        "data_alias": {"example": "survey"},
        "has_output": False,
        "verify_script": "scripts/hi_survey/verify.py",
        "verify_needs_stdout": True,
    },
    "nod_kfpa": {
        "dysh_script": "scripts/nod_kfpa/dysh_script.py",
        "gbtidl_script": None,
        "data_path": "/home/dysh/example_data/nod-KFPA/data/TGBT22A_503_02.raw.vegas",
        "data_alias": {"example": "nod"},
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


def _resolve_data_path(name: str, cfg: dict) -> str | None:
    """Resolve a benchmark's data path.

    Resolution order is the ``$DYSH_BENCH_DATA_PATH`` environment override, the canonical
    GBO path if it exists on this host, then the `dysh.util.files.dysh_data` alias.

    Parameters
    ----------
    name : str
        Benchmark name, used only in messages.
    cfg : dict
        The benchmark's entry in `BENCHMARKS`. Uses the ``data_path`` (canonical path or `None`)
        and ``data_alias`` (keyword arguments for `dysh_data`, optional) keys.

    Returns
    -------
    str or None
        The resolved path, or `None` if no data could be found. `None` is also returned for
        benchmarks that use no data (``data_path`` is `None`) and have no environment override.
    """
    env_path = os.environ.get("DYSH_BENCH_DATA_PATH")
    if env_path:
        return env_path
    gbo_path = cfg.get("data_path")
    if gbo_path is None:
        return None
    if Path(gbo_path).exists():
        return gbo_path
    alias = cfg.get("data_alias")
    if alias:
        from dysh.util.files import dysh_data

        resolved = dysh_data(**alias)
        if resolved is not None:
            console.print(f"[dim]{name}: canonical path missing; using dysh_data alias {alias} -> {resolved}[/]")
            return str(resolved)
    console.print(f"[yellow]warn:[/] no data found for [bold]{name}[/] (tried {gbo_path} and alias {alias})")
    return None


# ---------------------------------------------------------------------------
# Page-cache eviction
# ---------------------------------------------------------------------------


def _evict_from_pagecache(path: Path, verbose: bool = False) -> None:
    """Ask the OS to drop `path`, or every file under it, from the page cache.

    Calls ``posix_fadvise(DONTNEED)`` on each file. This is a hint, not a guarantee; use
    ``drop_caches`` as root for a guaranteed cold cache. Failures are reported as warnings,
    not raised.

    Parameters
    ----------
    path : `~pathlib.Path`
        A single file, or a directory to walk recursively.
    verbose : bool, optional
        Print each evicted file.
    """
    if not hasattr(os, "posix_fadvise"):
        console.print("[yellow]warn:[/] posix_fadvise not available on this platform — skipping eviction")
        return
    for fpath in [path] if path.is_file() else sorted(path.rglob("*")):
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
    """Check whether GBTIDL can be run on this host.

    Returns
    -------
    bool
        `True` if ``gbtidl`` is on ``PATH``.
    """
    return shutil.which("gbtidl") is not None


def _resolve_path(path: str | Path) -> Path:
    """Make a registry path absolute.

    Parameters
    ----------
    path : str or `~pathlib.Path`
        An absolute path, or one relative to `HERE`.

    Returns
    -------
    `~pathlib.Path`
        The absolute, normalized path (``..`` components removed). It need not exist.
    """
    path = Path(path)
    return path if path.is_absolute() else (HERE / path).resolve()


def _script_argv(script: Script) -> list[str]:
    """Normalize a registry script entry to an argv list with an absolute script path.

    Parameters
    ----------
    script : str or list of str
        A path, or an argv list whose first item is the path. Only the first item is resolved
        (see `_resolve_path`); the remaining items are arguments and are passed unchanged.

    Returns
    -------
    list of str
        The script's absolute path followed by its arguments.
    """
    argv = [script] if isinstance(script, str) else list(script)
    return [str(_resolve_path(argv[0])), *argv[1:]]


def _make_cmd(tool: str, script: Script) -> list[str]:
    """Build the command line that runs a benchmark script.

    dysh scripts run with ``--frozen`` so that ``uv run`` uses ``uv.lock`` as is and never
    re-locks or rewrites it, keeping the environment identical between the "before" and "after"
    runs of a performance comparison.

    Parameters
    ----------
    tool : {"dysh", "gbtidl"}
        Which tool runs the script. Any value other than ``"dysh"`` is treated as GBTIDL.
    script : str or list of str
        The script path (``.py`` for dysh, ``.pro`` for GBTIDL), or an argv list with the script
        path first and its arguments after. A relative path is resolved against `HERE`.

    Returns
    -------
    list of str
        ``["uv", "run", "--frozen", "python", *argv]`` for dysh, ``["gbtidl", *argv]`` otherwise.
    """
    argv = _script_argv(script)
    if tool == "dysh":
        return ["uv", "run", "--frozen", "python", *argv]
    else:
        return ["gbtidl", *argv]


# ---------------------------------------------------------------------------
# Single-iteration helpers
# ---------------------------------------------------------------------------


def _parse_script_seconds(output: str) -> float | None:
    """Extract the whole-script time from a script's output.

    Parameters
    ----------
    output : str
        Captured stdout containing a ``DYSH_BENCH_SCRIPT_MS=`` or ``GBTIDL_BENCH_SCRIPT_MS=``
        marker. IDL's ``D`` exponent is accepted.

    Returns
    -------
    float or None
        Script time in seconds (first marker found), or `None` if there is no marker.
    """
    match = SCRIPT_MS_RE.search(output)
    if match is None:
        return None
    return float(match.group("ms").replace("D", "E").replace("d", "e")) / 1000.0


def _parse_stage_seconds(output: str) -> dict[str, float]:
    """Extract per-stage times from a script's output.

    Parameters
    ----------
    output : str
        Captured stdout containing ``*_BENCH_STAGE_MS[<stage>]=<ms>`` markers.

    Returns
    -------
    dict
        Stage name -> time in seconds. If a stage name repeats, the last value wins.
    """
    stages: dict[str, float] = {}
    for match in STAGE_MS_RE.finditer(output):
        stages[match.group("stage")] = float(match.group("ms").replace("D", "E").replace("d", "e")) / 1000.0
    return stages


def _read_peak_rss_kb(pid: int) -> int:
    """Read a process's peak resident set size from procfs.

    Parameters
    ----------
    pid : int
        Process id.

    Returns
    -------
    int
        The larger of ``VmHWM`` and ``VmRSS`` in kB, or 0 if the process is gone or
        ``/proc/<pid>/status`` cannot be read (e.g. non-Linux).
    """
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
    """Copy lines from a stream into a list until the stream closes.

    Intended as a thread target so a child process cannot block on a full pipe.

    Parameters
    ----------
    stream : file-like
        Text stream to read, e.g. ``Popen.stdout``.
    output_lines : list of str
        List that each line is appended to.
    verbose : bool
        Also echo each line to the console.
    """
    for line in stream:
        output_lines.append(line)
        if verbose:
            console.print(line, end="")


def _run_script(cmd: list[str], env: dict, verbose: bool = False, require_script_marker: bool = True) -> dict:
    """Run one benchmark script, measuring wall time, peak memory, and its own timing markers.

    Parameters
    ----------
    cmd : list of str
        Command to run, from `_make_cmd`. It runs with `HERE` as its working directory.
    env : dict
        Environment for the child process.
    verbose : bool, optional
        Stream the child's output live.
    require_script_marker : bool, optional
        Treat a missing ``BENCH_SCRIPT_MS`` marker as a failure.

    Returns
    -------
    dict
        ``elapsed_s`` (wall time of the whole process), ``stdout`` (stdout and stderr combined),
        ``stderr`` (always empty; stderr is merged into stdout), ``script_body_s`` (from the
        marker, or `None`), ``stage_s`` (stage name -> seconds), and ``peak_rss_mb``.

    Raises
    ------
    SystemExit
        If the process exits nonzero, or if `require_script_marker` is set and no script
        marker was emitted.
    """
    if verbose:
        console.print(f"[dim]$ {' '.join(cmd)}[/]")
    t0 = time.perf_counter()
    proc = subprocess.Popen(
        cmd,
        env=env,
        cwd=HERE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    output_lines: list[str] = []
    if proc.stdout is None:
        raise RuntimeError("proc.stdout is None")
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
    script: Script,
    data_path: str | None,
    has_output: bool,
    tmpdir: str | None,
    verbose: bool,
    require_script_marker: bool,
) -> dict:
    """Run a script once against a fresh, page-cache-evicted copy of its data.

    The data (a directory tree or a single file) is copied to a unique temporary directory,
    evicted from the page cache, and the script is pointed at the copy through
    ``$DYSH_BENCH_DATA_PATH``. The copy is removed afterwards.

    Parameters
    ----------
    tool : {"dysh", "gbtidl"}
        Which tool runs the script.
    script : str or list of str
        The script path, or an argv list (see `_make_cmd`).
    data_path : str or None
        Data file or directory to copy. `None` for benchmarks that use no data.
    has_output : bool
        Set ``$DYSH_BENCH_OUT_PATH`` to a file inside the temporary directory.
    tmpdir : str or None
        Parent directory for the temporary copy; `None` uses the system default.
    verbose : bool
        Print progress and stream the child's output.
    require_script_marker : bool
        Treat a missing ``BENCH_SCRIPT_MS`` marker as a failure.

    Returns
    -------
    dict
        The result of `_run_script`.
    """
    tmp_root = Path(tempfile.mkdtemp(prefix=f"dysh_cold_{uuid.uuid4().hex[:8]}_", dir=tmpdir))
    try:
        env = dict(os.environ)
        if data_path:
            tmp_data = tmp_root / Path(data_path).name
            if verbose:
                console.print(f"  copying {data_path} -> {tmp_data}")
            if Path(data_path).is_dir():
                shutil.copytree(data_path, tmp_data)
            else:
                shutil.copy2(data_path, tmp_data)
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
    script: Script,
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
    """Run a script repeatedly and collect the timing of each run.

    In ``"warm"`` mode one untimed warm-up run precedes the timed runs, which then share the
    same data path and page cache. In ``"cold"`` mode each timed run uses a fresh evicted copy
    of the data (see `_run_one_cold`).

    Parameters
    ----------
    tool : {"dysh", "gbtidl"}
        Which tool runs the script.
    label : str
        Name shown in progress output.
    script : str or list of str
        The script path, or an argv list (see `_make_cmd`).
    data_path : str or None
        Data to run against; `None` for benchmarks that use no data.
    has_output : bool
        Set ``$DYSH_BENCH_OUT_PATH`` so the script can write an output file.
    cache_mode : {"warm", "cold"}
        Cache behavior.
    n_iterations : int
        Number of timed runs.
    tmpdir : str or None
        Parent directory for cold-mode copies.
    verbose : bool
        Stream script output instead of showing a progress bar.
    overall_task : int or None
        Task id in `progress` to advance, or `None` when there is no progress bar.
    progress : `rich.progress.Progress` or None
        Progress bar; `None` in verbose mode.
    require_script_marker : bool
        Treat a missing ``BENCH_SCRIPT_MS`` marker as a failure.

    Returns
    -------
    list of dict
        One `_run_script` result per timed run (the warm-up is not included).
    """
    env = dict(os.environ)
    if data_path:
        env["DYSH_BENCH_DATA_PATH"] = data_path
    if has_output:
        env["DYSH_BENCH_OUT_PATH"] = "bench_output.fits"

    cmd = _make_cmd(tool, script)

    def _run():
        """Run the script once in the warm-cache environment.

        Returns
        -------
        dict
            The result of `_run_script`.
        """
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


def _golden_path(cfg: dict) -> Path:
    """Locate the golden capture for a benchmark.

    Parameters
    ----------
    cfg : dict
        The benchmark's entry in `BENCHMARKS`; uses ``dysh_script``.

    Returns
    -------
    `~pathlib.Path`
        Absolute path of ``golden.txt`` next to the benchmark's dysh script. It may not exist.
    """
    return _resolve_path(_script_argv(cfg["dysh_script"])[0]).parent / "golden.txt"


def _capture_stdout(tool: str, name: str, script: Script, env: dict, verbose: bool) -> str:
    """Run one script and return its stdout.

    Parameters
    ----------
    tool : {"dysh", "gbtidl"}
        Which tool runs the script.
    name : str
        Benchmark name, used in messages.
    script : str or list of str
        The script path, or an argv list (see `_make_cmd`).
    env : dict
        Environment for the child process.
    verbose : bool
        Print a note that the capture is starting.

    Returns
    -------
    str
        The captured stdout.

    Raises
    ------
    SystemExit
        If the script exits nonzero; its output is printed first.
    """
    if verbose:
        console.print(f"[dim]verify: capturing {tool} stdout for {name}[/]")
    result = subprocess.run(_make_cmd(tool, script), env=env, cwd=HERE, capture_output=True, text=True)
    if result.returncode != 0:
        console.print(f"[red bold]VERIFY FAIL[/] {name} ({tool} exited with code {result.returncode})")
        if result.stdout:
            console.print(result.stdout, end="")
        if result.stderr:
            console.print(result.stderr, end="", style="red")
        sys.exit(1)
    return result.stdout


def _capture_golden(name: str, cfg: dict, data_path: str | None, verbose: bool) -> None:
    """Capture dysh's stdout on the current build as the golden reference.

    Writes ``golden.txt`` next to the dysh script. Benchmarks without a stdout verifier are
    skipped.

    Parameters
    ----------
    name : str
        Benchmark name.
    cfg : dict
        The benchmark's entry in `BENCHMARKS`.
    data_path : str or None
        Data path passed to the script through ``$DYSH_BENCH_DATA_PATH``.
    verbose : bool
        Print a note that the capture is starting.

    Raises
    ------
    SystemExit
        If the script exits nonzero.
    """
    if not cfg.get("verify_needs_stdout", False):
        console.print(f"[yellow]golden:[/] [bold]{name}[/] has no stdout verifier, skipping")
        return
    env = dict(os.environ)
    if data_path:
        env["DYSH_BENCH_DATA_PATH"] = data_path
    golden = _golden_path(cfg)
    golden.write_text(_capture_stdout("dysh", name, cfg["dysh_script"], env, verbose))
    console.print(f"[green]golden:[/] wrote {golden} (data: {data_path})")


def _run_verify(
    name: str,
    cfg: dict,
    data_path: str | None,
    has_gbtidl: bool,
    verbose: bool,
    dysh_stdout: str | None = None,
    gbtidl_stdout: str | None = None,
) -> None:
    """Run a benchmark's ``verify.py`` once, capturing script stdout as needed.

    The reference is GBTIDL when available, otherwise the committed golden capture. The check
    is skipped, with a message, if there is no ``verify.py``, the verifier does not use stdout,
    or neither GBTIDL nor a golden file is available. Output already captured by a timed run
    can be passed in, which avoids re-running the script.

    Parameters
    ----------
    name : str
        Benchmark name.
    cfg : dict
        The benchmark's entry in `BENCHMARKS`.
    data_path : str or None
        Data path passed to the scripts through ``$DYSH_BENCH_DATA_PATH``.
    has_gbtidl : bool
        Whether to run GBTIDL to produce the reference instead of using the golden file.
    verbose : bool
        Print progress notes.
    dysh_stdout : str, optional
        Output of an earlier dysh run of the script. If `None`, the script is run again.
    gbtidl_stdout : str, optional
        Output of an earlier GBTIDL run of the script. If `None` and `has_gbtidl`, the script is
        run again.

    Raises
    ------
    SystemExit
        If a script fails or the verifier reports a mismatch.
    """
    verify_script = cfg.get("verify_script")
    if not verify_script or not _resolve_path(verify_script).exists():
        console.print(f"[yellow]verify:[/] no verify.py for [bold]{name}[/], skipping")
        return

    if not cfg.get("verify_needs_stdout", False):
        console.print(f"[yellow]verify:[/] [bold]{name}[/] verifier does not use stdout, skipping")
        return

    env = dict(os.environ)
    if data_path:
        env["DYSH_BENCH_DATA_PATH"] = data_path

    golden = _golden_path(cfg)
    if not has_gbtidl and not golden.exists():
        console.print(f"[yellow]verify:[/] neither gbtidl nor {golden} available — cannot verify [bold]{name}[/]")
        return

    dysh_out = Path(tempfile.mktemp(suffix=f".{name}.dysh.out"))
    gbtidl_out = Path(tempfile.mktemp(suffix=f".{name}.gbtidl.out"))
    try:
        if dysh_stdout is None:
            dysh_stdout = _capture_stdout("dysh", name, cfg["dysh_script"], env, verbose)
        dysh_out.write_text(dysh_stdout)
        if has_gbtidl:
            if gbtidl_stdout is None:
                gbtidl_stdout = _capture_stdout("gbtidl", name, cfg["gbtidl_script"], env, verbose)
            gbtidl_out.write_text(gbtidl_stdout)
            reference = gbtidl_out
        else:
            console.print(f"[dim]verify: gbtidl not available, comparing against {golden}[/]")
            reference = golden

        result = subprocess.run(
            _make_cmd("dysh", [verify_script, str(dysh_out), str(reference)]),
            cwd=HERE,
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


# ---------------------------------------------------------------------------
# Stats + table
# ---------------------------------------------------------------------------


def _stats(times: list[float]) -> dict:
    """Summarize a list of timings.

    Parameters
    ----------
    times : list of float
        Timings in seconds (or any single unit).

    Returns
    -------
    dict
        ``n``, ``mean``, ``std`` (sample standard deviation, 0 for a single value), ``min``,
        ``max``, and the raw ``times``.
    """
    n = len(times)
    m = mean(times)
    s = stdev(times) if n > 1 else 0.0
    return {"n": n, "mean": m, "std": s, "min": min(times), "max": max(times), "times": times}


def _stats_from_runs(runs: list[dict]) -> dict:
    """Summarize the runs of one benchmark.

    Parameters
    ----------
    runs : list of dict
        Results from `_run_iterations`.

    Returns
    -------
    dict
        The `_stats` summary of the total wall time, plus these keys when available:
        ``peak_rss_mb`` (memory), ``script_body`` and ``startup_overhead`` (only if every run
        reported a script time), and ``stages`` (per-stage summaries, only for stages present
        in every run).
    """
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
    """Fill in timing for scripts that only measure process startup.

    For the ``exit`` benchmark the whole run is startup overhead, so the script body time is
    zero and the total time is the startup overhead.

    Parameters
    ----------
    stats : dict
        Result of `_stats_from_runs`. Modified in place.

    Returns
    -------
    dict
        The same `stats`, with ``script_body`` and ``startup_overhead`` set.
    """
    zeroes = [0.0] * stats["n"]
    stats["script_body"] = _stats(zeroes)
    stats["startup_overhead"] = _stats(list(stats["times"]))
    return stats


def _print_results(results: dict, modes: list[str], all_columns: bool = False) -> None:
    """Print the summary table and per-stage breakdown tables.

    Parameters
    ----------
    results : dict
        Nested ``results[benchmark][mode][tool]`` dictionaries of `_stats_from_runs` output.
    modes : list of str
        Cache modes to print, in order.
    all_columns : bool, optional
        Also show startup time, peak RSS, standard deviation, and startup speedup.
    """
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

            def _print_single_tool_stage_table(
                tool: str, stages: dict, script_mean: float | None, _name: str = name, _mode: str = mode
            ) -> None:
                """Print the stage breakdown of one tool, when there is no other tool to compare against.

                Parameters
                ----------
                tool : str
                    Tool name, used in the table title and column headers.
                stages : dict
                    Stage name -> `_stats` summary. Nothing is printed if empty.
                script_mean : float or None
                    Mean script time, used for the percentage column; the column is blank if `None` or 0.
                _name : str, optional
                    Benchmark name for the title. The default binds the loop variable of the enclosing
                    function; do not pass it.
                _mode : str, optional
                    Cache mode for the title. Bound like `_name`.
                """
                nonlocal stage_tables_printed
                if not stages:
                    return
                stage_tables_printed += 1
                stage_table = Table(
                    title=f"{_name} {_mode} {tool} stage breakdown",
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
    """Command-line entry point.

    Reads its options from ``sys.argv`` (see ``--help``): times the selected benchmarks with
    dysh and, if available, GBTIDL; optionally verifies results or captures golden files;
    prints the summary tables; and optionally writes JSON.

    Raises
    ------
    SystemExit
        If a benchmark script fails or verification fails.
    """
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
        help="after timing, run each benchmark's verify.py against gbtidl (or golden.txt if gbtidl is absent)",
    )
    parser.add_argument(
        "--capture-golden",
        action="store_true",
        help="instead of timing, capture dysh stdout as scripts/<name>/golden.txt for later --verify runs",
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

    # Resolve data up front; drop benchmarks whose data cannot be found.
    resolved_data: dict[str, str | None] = {}
    selected: list[str] = []
    for name in args.benchmarks:
        cfg = BENCHMARKS[name]
        resolved_data[name] = args.data_path or _resolve_data_path(name, cfg)
        if resolved_data[name] is None and cfg["data_path"] is not None:
            console.print(f"[yellow]skip:[/] [bold]{name}[/] — no data available on this host")
            continue
        selected.append(name)

    if args.capture_golden:
        for name in selected:
            _capture_golden(name, BENCHMARKS[name], resolved_data[name], args.verbose)
        return

    results: dict = {}

    # Compute total steps: each tool run = n_iterations timed + 1 warmup (warm mode only)
    steps_per_run = args.iterations + (1 if args.mode != "cold" else 0)
    n_total_steps = sum(
        steps_per_run * (1 + (1 if has_gbtidl and BENCHMARKS[n]["gbtidl_script"] else 0)) * len(modes) for n in selected
    )

    def _run_all(overall, progress):
        """Time every selected benchmark in every requested mode and store the results.

        Fills the enclosing function's ``results`` dictionary, and verifies each benchmark
        afterwards if ``--verify`` was given.

        Parameters
        ----------
        overall : int or None
            Task id of the overall progress bar, or `None` in verbose mode.
        progress : `rich.progress.Progress` or None
            Progress bar, or `None` in verbose mode.
        """
        for name in selected:
            cfg = BENCHMARKS[name]
            data_path = resolved_data[name]
            has_output = cfg.get("has_output", False)
            results[name] = {}
            dysh_stdout = gbtidl_stdout = None  # first timed run's output, reused by --verify

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
                if dysh_stdout is None and runs:
                    dysh_stdout = runs[0]["stdout"]
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
                    if gbtidl_stdout is None and runs:
                        gbtidl_stdout = runs[0]["stdout"]
                    gbtidl_stats = _stats_from_runs(runs)
                    if cfg.get("script_body_zero"):
                        gbtidl_stats = _apply_zero_script_body(gbtidl_stats)
                    results[name][mode]["gbtidl"] = gbtidl_stats

            if args.verify:
                _run_verify(name, cfg, data_path, has_gbtidl, args.verbose, dysh_stdout, gbtidl_stdout)

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
