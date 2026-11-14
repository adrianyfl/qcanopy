#!/usr/bin/env python3
"""Profile qcanopy::Simulator (src/sim/sim.cpp) on labeled circuit pairs.

Mirrors profiling/qcec: generate datasets with qcanopy_generate across a size
grid, run every pair through the simulator with one or more strategies, check
each result against the generator's proven label (and, for small sizes, an
independent reference simulation), and write a terminal summary plus a
self-contained HTML report.

    sim_profile.py run [--sweep FILE] [--set KEY=VALUE ...] [--jobs N] [--timeout S]
                       [--repeats N] [--out DIR | --run-dir DIR] [--shard K/N]
                       [--generator PATH] [--bench PATH] [--resume RUN_DIR]
    sim_profile.py summary [RUN_DIR]
    sim_profile.py report  [RUN_DIR]
    sim_profile.py merge   OUT_DIR SHARD_DIR [SHARD_DIR ...]

Standard library only (Python 3.11+ for tomllib).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import hashlib
import itertools
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import threading
import time
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DEFAULT_SWEEP = HERE / "sweeps" / "smoke.toml"
DEFAULT_OUT = REPO / "generated" / "sim-profile"
TEMPLATE = HERE / "report_template.html"
EXE = ".exe" if os.name == "nt" else ""
STRATEGY_KINDS = ("flat", "dag", "normalized")

# Generator options the sweep may set; keys map to --kebab-case flags.
GENERATOR_KEYS = (
    "qubits", "gates", "pairs", "prefix_fraction", "equivalent_fraction", "seed",
    "gate_set", "angle_mode", "prefix_group_size", "verify",
)
ABBREV = {
    "qubits": "q", "gates": "g", "pairs": "p", "prefix_fraction": "pf",
    "equivalent_fraction": "ef", "seed": "s", "gate_set": "", "angle_mode": "",
    "prefix_group_size": "pg", "verify": "v",
}
# Exit codes of a worker stopped by Ctrl+C (POSIX SIGINT, Windows STATUS_CONTROL_C_EXIT).
INTERRUPTED = (-2, 130, 0xC000013A, -1073741510)


class Interrupted(Exception):
    pass


# Workers currently running, so an interrupt (Ctrl+C, or Slurm's warning signal
# before the time limit) can stop them instead of waiting for them to finish.
_LIVE: set = set()
_LIVE_LOCK = threading.Lock()
_STOPPING = threading.Event()


def run_worker(cmd: list[str], timeout: float):
    """subprocess.run(capture_output=True, text=True, timeout=...) that can be stopped."""
    if _STOPPING.is_set():
        raise Interrupted()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with _LIVE_LOCK:
        _LIVE.add(proc)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        raise
    finally:
        with _LIVE_LOCK:
            _LIVE.discard(proc)
    if _STOPPING.is_set():
        raise Interrupted()  # killed by stop_workers(); leave the check for --resume
    return subprocess.CompletedProcess(cmd, proc.returncode, out, err)


def stop_workers() -> None:
    _STOPPING.set()
    with _LIVE_LOCK:
        for proc in list(_LIVE):
            try:
                proc.kill()
            except OSError:
                pass


TERMINAL = ("correct", "undetected", "incorrect", "timeout", "error", "skipped", "skipped_memory")


# ------------------------------------------------------------------ helpers

def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def find_executable(name: str, explicit: str | None, subdir: str) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.is_file():
            sys.exit(f"error: {p} not found")
        return p.resolve()
    candidates = []
    for build in sorted(REPO.glob("build*")):
        for cfg in ("", "Release", "RelWithDebInfo"):
            candidates.append(build / subdir / cfg / (name + EXE))
    found = [c for c in candidates if c.is_file()]
    if not found:
        sys.exit(f"error: could not find {name}{EXE} under {REPO}/build*/{subdir}; "
                 f"build it (see profiling/sim/README.md) or pass its path explicitly")
    return max(found, key=lambda p: p.stat().st_mtime).resolve()


def physical_memory_mb() -> float | None:
    try:
        if hasattr(os, "sysconf") and "SC_PHYS_PAGES" in os.sysconf_names:
            return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2**20
        if os.name == "nt":
            import ctypes

            class MS(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
            ms = MS()
            ms.dwLength = ctypes.sizeof(MS)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
            return ms.ullTotalPhys / 2**20
    except Exception:
        pass
    return None


def allocated_memory_mb() -> tuple[float | None, str]:
    """Memory this process may use: the Slurm allocation when inside a job, else physical RAM."""
    env = os.environ
    try:
        if env.get("SLURM_MEM_PER_NODE"):
            return float(env["SLURM_MEM_PER_NODE"]), "SLURM_MEM_PER_NODE"
        if env.get("SLURM_MEM_PER_CPU"):
            cpus = int(env.get("SLURM_CPUS_ON_NODE") or env.get("SLURM_CPUS_PER_TASK") or 1)
            return float(env["SLURM_MEM_PER_CPU"]) * cpus, "SLURM_MEM_PER_CPU x CPUs"
    except ValueError:
        pass
    return physical_memory_mb(), "physical RAM"


def slurm_info() -> dict | None:
    keys = ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_PARTITION",
            "SLURMD_NODENAME", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU")
    info = {k: os.environ[k] for k in keys if k in os.environ}
    return info or None


def cpu_name() -> str:
    try:
        if sys.platform.startswith("linux"):
            for line in Path("/proc/cpuinfo").read_text().splitlines():
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
        if os.name == "nt":
            return os.environ.get("PROCESSOR_IDENTIFIER", platform.processor())
        if sys.platform == "darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                  capture_output=True, text=True).stdout.strip()
    except Exception:
        pass
    return platform.processor() or "unknown"


def git_commit() -> str:
    try:
        r = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return r.stdout.strip() + ("+dirty" if dirty else "") if r.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def file_sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16] if path.is_file() else None


def median(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


# ------------------------------------------------------------------ sweep

def load_sweep(path: Path, args) -> dict:
    with open(path, "rb") as f:
        sw = tomllib.load(f)
    for item in args.set or []:
        apply_override(sw, item)
    sw.setdefault("name", path.stem)
    sw.setdefault("repeats", 3)
    sw.setdefault("timeout", 120.0)
    sw.setdefault("warmup", True)
    sw.setdefault("probes", 3)
    sw.setdefault("tol", 1e-8)
    sw.setdefault("reference_max_qubits", 12)
    sw.setdefault("reference_tol", 1e-8)
    sw.setdefault("kernel_scan", True)
    sw.setdefault("generator", {})
    sw.setdefault("strategies", {"flat": {}})
    if "memory_budget_mb" not in sw:
        mem, source = allocated_memory_mb()
        # Inside a Slurm allocation the limit is enforced, so use most of it;
        # on a shared workstation leave half for everything else.
        frac = 0.85 if source.startswith("SLURM") else 0.5
        sw["memory_budget_mb"] = round(mem * frac) if mem else 4096
        sw["memory_budget_source"] = f"{frac:g} x {source}" if mem else "fallback"
    if args.timeout is not None:
        sw["timeout"] = args.timeout
    if args.repeats is not None:
        sw["repeats"] = args.repeats
    unknown = set(sw["generator"]) - set(GENERATOR_KEYS)
    if unknown:
        sys.exit(f"error: unknown [generator] keys in {path}: {sorted(unknown)}")
    for name, st in sw["strategies"].items():
        kind = st.get("strategy", name if name in STRATEGY_KINDS else None)
        if kind not in STRATEGY_KINDS:
            sys.exit(f"error: strategy '{name}' needs strategy = one of {STRATEGY_KINDS}")
        st["strategy"] = kind
    return sw


def apply_override(sw: dict, item: str) -> None:
    """--set a.b.c=VALUE, VALUE in TOML syntax: 3600, [20, 24, 28], "mixed", true, {}."""
    if "=" not in item:
        sys.exit(f"error: --set needs KEY=VALUE, got {item!r}")
    key, raw = item.split("=", 1)
    try:
        value = tomllib.loads(f"v = {raw.strip()}")["v"]
    except tomllib.TOMLDecodeError:
        value = raw.strip()  # bare words such as clifford-t
    parts = key.strip().split(".")
    node = sw
    for part in parts[:-1]:
        node = node.setdefault(part, {})
        if not isinstance(node, dict):
            sys.exit(f"error: --set {key}: {part} is not a table")
    node[parts[-1]] = value


def series_of(ds: dict) -> str:
    return json.dumps({k: v for k, v in ds.items() if k not in ("gates", "pairs", "id")}, sort_keys=True)


def shard_filter(datasets: list[dict], shard: tuple[int, int] | None) -> list[dict]:
    """Keep the datasets of shard K of N. A series (every gate count of one
    qubit count / gate set / ...) stays in one shard so auto-skip still works."""
    if not shard:
        return datasets
    k, n = shard
    series = sorted({series_of(d) for d in datasets},
                    key=lambda s: (json.loads(s).get("qubits", 0), s))
    mine = {s for i, s in enumerate(series) if i % n == k}
    return [d for d in datasets if series_of(d) in mine]


def parse_shard(text: str | None) -> tuple[int, int] | None:
    if not text:
        return None
    try:
        k, n = (int(x) for x in text.split("/"))
    except ValueError:
        sys.exit(f"error: --shard must look like K/N, got {text!r}")
    if not (n >= 1 and 0 <= k < n):
        sys.exit(f"error: --shard {text}: need 0 <= K < N")
    return k, n


def expand_grid(sw: dict) -> list[dict]:
    gen = sw["generator"]
    keys = [k for k in GENERATOR_KEYS if k in gen]
    values = [gen[k] if isinstance(gen[k], list) else [gen[k]] for k in keys]
    varying = {k for k, v in zip(keys, values) if len(v) > 1}
    datasets = []
    for combo in itertools.product(*values):
        params = dict(zip(keys, combo))
        params.setdefault("qubits", 4)
        params.setdefault("gates", 100)
        params.setdefault("pairs", 32)
        budget = sw.get("pair_gate_budget")
        if budget:
            params["pairs"] = max(2, min(params["pairs"], budget // max(1, params["gates"])))
        parts = [f"q{params['qubits']}", f"g{params['gates']}"]
        parts += [f"{ABBREV[k]}{params[k]}" for k in keys if k in varying and k not in ("qubits", "gates")]
        params["id"] = "-".join(str(p).replace(".", "p") for p in parts)
        datasets.append(params)
    datasets.sort(key=lambda d: (d["gates"], d["qubits"], d["id"]))
    return datasets


def series_key(ds: dict, strategy: str) -> str:
    rest = {k: v for k, v in ds.items() if k not in ("gates", "pairs", "id")}
    return strategy + "|" + json.dumps(rest, sort_keys=True)


def state_mb(qubits: int) -> float:
    return 16 * 2**qubits / 2**20


def peak_estimate_mb(qubits: int) -> float:
    # Two live state vectors (A and B) plus a reference or a returned temporary.
    return 3 * state_mb(qubits) + 64


# ------------------------------------------------------------------ manifest

def _walk(obj, path=()):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _walk(v, path + (str(k),))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk(v, path + (str(i),))
    else:
        yield path, obj


def _label(rec: dict) -> str | None:
    for path, v in _walk(rec):
        key = path[-1].lower() if path else ""
        if isinstance(v, str) and any(t in key for t in ("relation", "expected", "label")):
            s = v.lower().replace("-", "_")
            if "not" in s or "non" in s or "inequiv" in s:
                return "not_equivalent"
            if "equiv" in s:
                return "equivalent"
    for path, v in _walk(rec):
        if isinstance(v, bool) and "equivalent" in (path[-1].lower() if path else ""):
            return "equivalent" if v else "not_equivalent"
    return None


def _side(path: tuple, value: str) -> str | None:
    import re
    tokens = [t.lower() for p in path for t in re.split(r"[_\-./]", p) if t]
    if tokens and tokens[-1] in ("a", "b"):
        return tokens[-1]
    for t in reversed(tokens):
        if t in ("a", "b"):
            return t
    stem = Path(value).stem.lower()
    if stem.endswith(("_a", "-a", ".a")) or stem == "a":
        return "a"
    if stem.endswith(("_b", "-b", ".b")) or stem == "b":
        return "b"
    return None


def read_manifest(dataset_dir: Path) -> list[dict]:
    """Return [{case_id, label, a, b}] from qcanopy_generate's manifest.jsonl.

    The parser looks for the documented content (expected relation, case id,
    relative circuit paths) by key names rather than one exact layout, so it
    keeps working if field names are refined. If it cannot find something it
    stops with the offending record, and this function is the place to adapt.
    """
    manifest = dataset_dir / "manifest.jsonl"
    if not manifest.is_file():
        raise RuntimeError(f"{manifest} missing (generation failed or incomplete)")
    cases = []
    for lineno, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        rec = json.loads(line)
        files = [(p, v) for p, v in _walk(rec)
                 if isinstance(v, str) and v.lower().endswith(".json") and "manifest" not in v.lower()]
        if not files:
            continue  # header record
        label = _label(rec)
        sides = {}
        for p, v in files:
            s = _side(p, v)
            if s and s not in sides:
                sides[s] = v
        if len(sides) < 2 and len(files) >= 2:
            sides = {"a": files[0][1], "b": files[1][1]}
        case_id = next((v for p, v in _walk(rec)
                        if isinstance(v, str) and v.startswith("case-") and not v.endswith(".json")), None)
        if label is None or len(sides) < 2:
            raise RuntimeError(f"{manifest}:{lineno}: could not find label and A/B paths in "
                               f"{line[:300]}; adapt read_manifest() in sim_profile.py")
        cases.append({
            "case_id": case_id or f"case-{len(cases):06d}",
            "label": label,
            "a": str((dataset_dir / sides["a"]).resolve()),
            "b": str((dataset_dir / sides["b"]).resolve()),
        })
    if not cases:
        raise RuntimeError(f"{manifest}: no case records")
    return cases


def generate(generator: Path, ds: dict, out: Path, timeout: float) -> None:
    if (out / "manifest.jsonl").is_file():
        return
    if out.exists():
        shutil.rmtree(out)  # an interrupted earlier attempt; the generator needs an empty dir
    cmd = [str(generator), "--format", "ir", "--output", str(out)]
    for k in GENERATOR_KEYS:
        if k in ds:
            cmd += ["--" + k.replace("_", "-"), str(ds[k])]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"qcanopy_generate failed ({r.returncode}): {(r.stderr or r.stdout).strip()[:2000]}")


# ------------------------------------------------------------------ checks

def classify(res: dict, label: str, sw: dict) -> tuple[str, str | None]:
    if res.get("error"):
        return "error", None
    reasons = []
    ref_tol = sw["reference_tol"]
    for side in ("a", "b"):
        e = res.get(f"ref_err_{side}")
        if e is not None and e > ref_tol:
            reasons.append(f"state of {side.upper()} differs from reference by {e:.3g}")
    if (res.get("max_norm_err") or 0) > 1e-6:
        reasons.append(f"state norm drifted by {res['max_norm_err']:.3g}")
    verdict = res.get("verdict")
    if label == "equivalent" and verdict != "equivalent":
        reasons.append(f"labeled equivalent but min fidelity {res.get('min_fidelity'):.6g}, "
                       f"phase spread {res.get('phase_spread'):.3g}")
    if reasons:
        return "incorrect", "; ".join(reasons)
    if label == "not_equivalent" and verdict == "equivalent":
        return "undetected", "labeled not equivalent but every probe state matched"
    return "correct", None


def run_check(bench: Path, case: dict, ds: dict, strat_name: str, st: dict, sw: dict) -> dict:
    timeout = float(st.get("timeout", sw["timeout"]))
    cmd = [str(bench), "--a", case["a"], "--b", case["b"], "--strategy", st["strategy"],
           "--repeats", str(st.get("repeats", sw["repeats"])),
           "--probes", str(st.get("probes", sw["probes"])),
           "--probe-seed", str(ds.get("seed", 759)),
           "--tol", str(sw["tol"]),
           "--reference-max-qubits", str(st.get("reference_max_qubits", sw["reference_max_qubits"]))]
    if not sw["warmup"]:
        cmd.append("--no-warmup")
    rec = base_record(ds, case, strat_name, st)
    t0 = time.perf_counter()
    try:
        p = run_worker(cmd, timeout)
    except subprocess.TimeoutExpired:
        rec.update(outcome="timeout", process_s=time.perf_counter() - t0,
                   error=f"killed after {timeout:g} s")
        return rec
    rec["process_s"] = time.perf_counter() - t0
    if p.returncode in INTERRUPTED:
        raise Interrupted()  # Ctrl+C reached the worker too; leave the check for --resume
    try:
        res = json.loads(p.stdout.strip().splitlines()[-1])
    except Exception:
        res = {"error": f"exit {p.returncode}; no JSON on stdout; stderr: {p.stderr.strip()[:1000]}"}
    if p.returncode != 0 and not res.get("error"):
        res["error"] = f"exit {p.returncode}: {p.stderr.strip()[:1000]}"
    outcome, reason = classify(res, case["label"], sw)
    rec.update(outcome=outcome, reason=reason)
    rec.update({k: v for k, v in res.items() if k not in ("mode", "n", "strategy")})
    if outcome != "error":
        n = res["n"]
        sa, sb = res["sim_a_median_s"], res["sim_b_median_s"]
        ga, gb = res["simulated_gates_a"], res["simulated_gates_b"]
        rec["sim_s"] = (sa + sb) / 2  # per circuit
        work = (ga + gb) * 2**n
        rec["ns_per_gate"] = (sa + sb) / max(1, ga + gb) * 1e9
        rec["ns_per_gate_amp"] = (sa + sb) / work * 1e9 if work else None
        rec["gb_s"] = 32 * work / (sa + sb) / 1e9 if (sa + sb) > 0 else None
        rec["prep_s"] = (res["prep_a_s"] + res["prep_b_s"]) / 2 if st["strategy"] != "flat" else None
    return rec


def base_record(ds, case, strat_name, st) -> dict:
    rec = {"dataset": ds["id"], "case_id": case["case_id"], "strategy": strat_name,
           "kind": st["strategy"], "label": case["label"]}
    rec.update({k: ds[k] for k in GENERATOR_KEYS if k in ds})
    return rec


# ------------------------------------------------------------------ run

class Recorder:
    def __init__(self, path: Path):
        self.path = path
        self.lock = threading.Lock()

    def write(self, rec: dict) -> None:
        with self.lock, open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")


def load_results(path: Path) -> list[dict]:
    """results.jsonl with retried checks collapsed to their latest record."""
    latest = {}
    for r in load_jsonl(path):
        latest[(r.get("dataset"), r.get("case_id"), r.get("strategy"))] = r
    return list(latest.values())


def load_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass  # a line cut off by an interrupted run
    return out


def cmd_run(args) -> None:
    import signal
    # Slurm sends SIGTERM at the time limit (or a chosen signal earlier with
    # --signal); treat it like Ctrl+C so the run can be resumed cleanly.
    def on_term(signum, frame):
        raise KeyboardInterrupt
    for name in ("SIGTERM", "SIGUSR1"):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), on_term)
    try:
        _cmd_run(args)
    except (KeyboardInterrupt, Interrupted):
        stop_workers()
        log("\ninterrupted; finish with: sim_profile.py run --resume <run dir printed above>")
        sys.exit(130)


def _cmd_run(args) -> None:
    generator = find_executable("qcanopy_generate", args.generator, "src/generator")
    bench = find_executable("sim_bench", args.bench, "profiling/sim")

    if args.resume:
        run_dir = Path(args.resume).resolve()
        info = json.loads((run_dir / "run.json").read_text())
        sweep_path = Path(info["sweep_file"])
        sw = info["sweep"]
        shard = tuple(info["shard"]) if info.get("shard") else None
        if args.set or args.shard:
            sys.exit("error: --set and --shard cannot change a run being resumed")
        if args.timeout is not None:
            sw["timeout"] = args.timeout
    else:
        sweep_path = Path(args.sweep or DEFAULT_SWEEP).resolve()
        sw = load_sweep(sweep_path, args)
        shard = parse_shard(args.shard)
        if args.run_dir:
            run_dir = Path(args.run_dir).resolve()
            if run_dir.exists() and any(run_dir.iterdir()):
                sys.exit(f"error: {run_dir} is not empty; use --resume {run_dir} to continue it")
        else:
            stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
            suffix = f"-shard{shard[0]}of{shard[1]}" if shard else ""
            run_dir = Path(args.out or DEFAULT_OUT).resolve() / f"{stamp}-{sw['name']}{suffix}"
        run_dir.mkdir(parents=True, exist_ok=True)
        sim_src = REPO / "src" / "sim" / "sim.cpp"
        info = {
            "sweep_file": str(sweep_path), "sweep": sw,
            "started": dt.datetime.now().isoformat(timespec="seconds"),
            "git_commit": git_commit(), "sim_cpp_sha256": file_sha256(sim_src),
            "bench": str(bench), "generator": str(generator),
            "host": {"platform": platform.platform(), "machine": platform.machine(),
                     "cpu": cpu_name(), "logical_cpus": os.cpu_count(),
                     "memory_mb": physical_memory_mb(), "python": platform.python_version(),
                     "node": platform.node(), "slurm": slurm_info()},
            "jobs": args.jobs, "shard": list(shard) if shard else None,
        }
        (run_dir / "run.json").write_text(json.dumps(info, indent=2))

    results = Recorder(run_dir / "results.jsonl")
    kernels = Recorder(run_dir / "kernels.jsonl")
    # Completed checks are skipped on --resume; errors are retried.
    done = {(r["dataset"], r["case_id"], r["strategy"]) for r in load_jsonl(results.path)
            if r.get("outcome") != "error"}
    kernels_done = {r["n"] for r in load_jsonl(kernels.path) if "n" in r}
    all_datasets = expand_grid(sw)
    datasets = shard_filter(all_datasets, shard)
    budget = sw["memory_budget_mb"]
    threshold = sw.get("skip_after_timeouts")

    # Auto-skip state: timeouts[(series, gates)] and blocked[series] = from gates.
    timeouts: dict = {}
    blocked: dict = {}
    for r in load_jsonl(results.path):
        if r.get("outcome") == "timeout" and threshold:
            for ds in datasets:
                if ds["id"] == r["dataset"]:
                    key = (series_key(ds, r["strategy"]), ds["gates"])
                    timeouts[key] = timeouts.get(key, 0) + 1
                    if timeouts[key] >= threshold:
                        s = key[0]
                        blocked[s] = min(blocked.get(s, math.inf), ds["gates"])
    state_lock = threading.Lock()

    t_start = time.perf_counter()
    total_checks = 0
    log(f"run: {run_dir}")
    log(f"bench: {bench}\ngenerator: {generator}\nmemory budget: {budget} MB "
        f"({sw.get('memory_budget_source', 'set in sweep')}), jobs: {args.jobs}")
    if shard:
        log(f"shard {shard[0]}/{shard[1]}: {len(datasets)} of {len(all_datasets)} datasets")

    # Kernel scans: one per distinct qubit count that fits.
    if sw["kernel_scan"]:
        scan_qubits = sw["kernel_scan"] if isinstance(sw["kernel_scan"], list) else \
            sorted({d["qubits"] for d in all_datasets})
        if shard:
            scan_qubits = [q for i, q in enumerate(sorted(scan_qubits)) if i % shard[1] == shard[0]]
        for n in scan_qubits:
            if n in kernels_done or peak_estimate_mb(n) > budget:
                continue
            log(f"kernel scan  n={n}")
            try:
                p = run_worker([str(bench), "--kernel-scan", "--qubits", str(n), "--repeats", "3"],
                               sw["timeout"] * 4)
                rec = json.loads(p.stdout.strip().splitlines()[-1])
                rec.setdefault("n", n)
            except subprocess.TimeoutExpired:
                rec = {"n": n, "error": "timeout"}
            except Interrupted:
                raise
            except Exception as e:  # noqa: BLE001
                rec = {"n": n, "error": str(e)}
            kernels.write(rec)

    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for ds in datasets:
            fits = peak_estimate_mb(ds["qubits"]) <= budget
            ds_dir = run_dir / "datasets" / ds["id"]
            pending = [(name, st) for name, st in sw["strategies"].items()]
            if fits:
                log(f"dataset {ds['id']}: generating ({ds['pairs']} pairs)")
                try:
                    generate(generator, ds, ds_dir, timeout=max(600.0, sw["timeout"]))
                    cases = read_manifest(ds_dir)
                except Exception as e:  # noqa: BLE001
                    log(f"  ERROR {e}")
                    for name, st in pending:
                        results.write({**base_record(ds, {"case_id": "dataset", "label": None}, name, st),
                                       "outcome": "error", "error": str(e)})
                    continue
            else:
                cases = [{"case_id": f"case-{i:06d}", "label": None} for i in range(ds["pairs"])]

            futures = []
            for name, st in pending:
                skey = series_key(ds, name)
                for case in cases:
                    if (ds["id"], case["case_id"], name) in done:
                        continue
                    total_checks += 1
                    if not fits:
                        results.write({**base_record(ds, case, name, st), "outcome": "skipped_memory",
                                       "reason": f"needs ~{peak_estimate_mb(ds['qubits']):.0f} MB, "
                                                 f"budget {budget} MB"})
                        continue

                    def task(case=case, name=name, st=st, skey=skey):
                        with state_lock:
                            if skey in blocked and ds["gates"] >= blocked[skey]:
                                rec = {**base_record(ds, case, name, st), "outcome": "skipped",
                                       "reason": f"auto-skip after {threshold} timeouts at "
                                                 f"{blocked[skey]} gates"}
                                results.write(rec)
                                return rec
                        rec = run_check(bench, case, ds, name, st, sw)
                        results.write(rec)
                        if rec["outcome"] == "timeout" and threshold:
                            with state_lock:
                                k = (skey, ds["gates"])
                                timeouts[k] = timeouts.get(k, 0) + 1
                                if timeouts[k] >= threshold:
                                    blocked[skey] = min(blocked.get(skey, math.inf), ds["gates"])
                        return rec

                    futures.append(pool.submit(task))
            counts: dict = {}
            try:
                for f in cf.as_completed(futures):
                    o = f.result()["outcome"]
                    counts[o] = counts.get(o, 0) + 1
            except BaseException:
                pool.shutdown(wait=False, cancel_futures=True)
                raise
            if futures:
                log("  " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
            if not fits:
                log(f"dataset {ds['id']}: skipped (state needs ~{peak_estimate_mb(ds['qubits']):.0f} MB "
                    f"> budget {budget} MB)")

    info = json.loads((run_dir / "run.json").read_text())
    info["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    info["elapsed_s"] = info.get("elapsed_s", 0) + (time.perf_counter() - t_start)
    (run_dir / "run.json").write_text(json.dumps(info, indent=2))
    print_summary(run_dir)
    write_report(run_dir)
    log(f"\nreport: {run_dir / 'report.html'}")


# ------------------------------------------------------------------ summary

def latest_run(arg: str | None) -> Path:
    if arg:
        return Path(arg).resolve()
    runs = sorted(p for p in DEFAULT_OUT.glob("*") if (p / "run.json").is_file())
    if not runs:
        sys.exit(f"error: no runs under {DEFAULT_OUT}")
    return runs[-1]


def table(headers: list[str], rows: list[list]) -> str:
    cells = [[str(h) for h in headers]] + [["" if c is None else str(c) for c in r] for r in rows]
    widths = [max(len(r[i]) for r in cells) for i in range(len(headers))]
    right = [i > 0 and all(_isnum(r[i]) for r in cells[1:] if r[i]) for i in range(len(headers))]
    fmt = lambda r: "  ".join(c.rjust(w) if right[i] else c.ljust(w) for i, (c, w) in enumerate(zip(r, widths)))
    return "\n".join([fmt(cells[0]), "  ".join("-" * w for w in widths)] + [fmt(r) for r in cells[1:]])


def _isnum(s: str) -> bool:
    try:
        float(s.replace(",", ""))
        return True
    except ValueError:
        return False


def fmt_time(s):
    if s is None:
        return None
    if s < 1e-3:
        return f"{s * 1e6:.1f}us"
    if s < 1:
        return f"{s * 1e3:.2f}ms"
    return f"{s:.2f}s"


def fmt(v, spec):
    return None if v is None else format(v, spec)


def print_summary(run_dir: Path) -> None:
    info = json.loads((run_dir / "run.json").read_text())
    recs = load_results(run_dir / "results.jsonl")
    print(f"\n{info['sweep']['name']}  ({run_dir.name})")
    print(f"commit {info.get('git_commit')}  sim.cpp {info.get('sim_cpp_sha256')}  "
          f"{info['host']['cpu']}, {info['host']['logical_cpus']} threads")

    strategies = list(dict.fromkeys(r["strategy"] for r in recs))
    print("\nOutcomes")
    rows = []
    for s in strategies:
        rs = [r for r in recs if r["strategy"] == s]
        rows.append([s, len(rs)] + [sum(r["outcome"] == o for r in rs) or "" for o in TERMINAL])
    print(table(["strategy", "checks", *TERMINAL], rows))

    print("\nTiming (medians over pairs; time is per circuit)")
    groups: dict = {}
    for r in recs:
        if r.get("sim_s") is None:
            continue
        k = (r["strategy"], r["qubits"], r["gates"], r.get("gate_set", ""), r.get("angle_mode", ""))
        groups.setdefault(k, []).append(r)
    rows = []
    for k in sorted(groups, key=lambda k: (strategies.index(k[0]), k[1], k[2], k[3], k[4])):
        rs = groups[k]
        rows.append([k[0], k[1], k[2], k[3] or None, len(rs),
                     fmt_time(median(r["sim_s"] for r in rs)),
                     fmt(median(r["ns_per_gate_amp"] for r in rs), ".3f"),
                     fmt(median(r["gb_s"] for r in rs), ".2f"),
                     fmt_time(median(r["prep_s"] for r in rs)),
                     fmt_time(median(r["load_s"] for r in rs)),
                     fmt(median(r.get("peak_rss_mb") for r in rs), ".1f")])
    print(table(["strategy", "qubits", "gates", "gate set", "pairs", "sim", "ns/gate·amp",
                 "GB/s", "prep", "load", "peak MB"], rows))

    ks = [k for k in load_jsonl(run_dir / "kernels.jsonl") if k.get("rows")]
    if ks:
        print("\nKernel scan (ns per amplitude; median over targets, [min..max])")
        kinds = ["h", "z", "cx", "cz", "swap"]
        rows = []
        for k in sorted(ks, key=lambda k: k["n"]):
            row = [k["n"]]
            for kind in kinds:
                v = [x["ns_per_amp"] for x in k["rows"] if x["kind"] == kind and x["ns_per_amp"] is not None]
                row.append(f"{statistics.median(v):.2f} [{min(v):.2f}..{max(v):.2f}]" if v else None)
            rows.append(row)
        print(table(["qubits", *kinds], rows))

    bad = [r for r in recs if r["outcome"] in ("incorrect", "undetected", "error")]
    if bad:
        print(f"\nNeeds attention ({len(bad)})")
        rows = [[r["outcome"], r["strategy"], r["dataset"], r["case_id"], r.get("label"),
                 (r.get("reason") or r.get("error") or "")[:90]] for r in bad[:25]]
        print(table(["outcome", "strategy", "dataset", "case", "label", "detail"], rows))
        if len(bad) > 25:
            print(f"... {len(bad) - 25} more in results.jsonl")


# ------------------------------------------------------------------ report

def write_report(run_dir: Path) -> None:
    info = json.loads((run_dir / "run.json").read_text())
    data = {"run": info, "results": load_results(run_dir / "results.jsonl"),
            "kernels": load_jsonl(run_dir / "kernels.jsonl"), "dir": run_dir.name}
    html = TEMPLATE.read_text(encoding="utf-8")
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    (run_dir / "report.html").write_text(html.replace("/*__DATA__*/null", payload), encoding="utf-8")


# ------------------------------------------------------------------ merge

def cmd_merge(out: Path, shards: list[Path]) -> None:
    """Combine shard runs (e.g. one per Slurm array task) into one run directory."""
    infos = []
    for d in shards:
        if not (d / "run.json").is_file():
            sys.exit(f"error: {d} has no run.json")
        infos.append(json.loads((d / "run.json").read_text()))
    if out.exists() and any(out.iterdir()):
        sys.exit(f"error: {out} is not empty")
    out.mkdir(parents=True, exist_ok=True)
    seen_shards = sorted(tuple(i["shard"]) for i in infos if i.get("shard"))
    expected = infos[0].get("shard", [0, 1])[1] if infos[0].get("shard") else 1
    missing = sorted(set(range(expected)) - {s[0] for s in seen_shards})
    if missing:
        log(f"warning: shards {missing} of {expected} are missing; the report will be partial")
    unfinished = [str(d) for d, i in zip(shards, infos) if not i.get("finished")]
    if unfinished:
        log(f"warning: unfinished shard runs (resume them first?): {unfinished}")
    merged = dict(infos[0])
    merged.update(shard=None, merged_from=[str(d) for d in shards],
                  started=min(i["started"] for i in infos),
                  finished=max((i.get("finished") or "") for i in infos) or None,
                  elapsed_s=max(i.get("elapsed_s") or 0 for i in infos),
                  hosts=sorted({(i["host"].get("node") or "?") for i in infos}),
                  shard_runs=[{"dir": str(d), "shard": i.get("shard"), "host": i["host"]}
                              for d, i in zip(shards, infos)])
    if len({i.get("sim_cpp_sha256") for i in infos}) > 1 or len({i["host"]["cpu"] for i in infos}) > 1:
        log("warning: shards differ in sim.cpp hash or CPU model; timings are not directly comparable")
    (out / "run.json").write_text(json.dumps(merged, indent=2))
    for name in ("results.jsonl", "kernels.jsonl"):
        with open(out / name, "w", encoding="utf-8") as f:
            for d in shards:
                for rec in (load_results if name == "results.jsonl" else load_jsonl)(d / name):
                    f.write(json.dumps(rec) + "\n")
    print_summary(out)
    write_report(out)
    log(f"\nreport: {out / 'report.html'}")


# ------------------------------------------------------------------ main

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="generate datasets and profile the simulator")
    r.add_argument("--sweep", help=f"sweep TOML (default {DEFAULT_SWEEP.relative_to(REPO)})")
    r.add_argument("--jobs", type=int, default=1, help="concurrent checks (default 1; keep 1 for timing)")
    r.add_argument("--timeout", type=float, help="hard per-check timeout in seconds")
    r.add_argument("--repeats", type=int, help="timed executions per circuit")
    r.add_argument("--out", help=f"parent directory for runs (default {DEFAULT_OUT.relative_to(REPO)})")
    r.add_argument("--generator", help="path to qcanopy_generate")
    r.add_argument("--bench", help="path to sim_bench")
    r.add_argument("--resume", metavar="RUN_DIR", help="finish an interrupted run")
    r.add_argument("--set", action="append", metavar="KEY=VALUE",
                   help="override a sweep value, TOML syntax; repeatable. "
                        "e.g. --set generator.qubits=[24,26,28] --set timeout=3600")
    r.add_argument("--shard", metavar="K/N", help="run only shard K (0-based) of N, e.g. a Slurm array task")
    r.add_argument("--run-dir", help="exact run directory (instead of a timestamped one under --out)")
    s = sub.add_parser("summary", help="print the terminal tables for a run")
    s.add_argument("run_dir", nargs="?")
    p = sub.add_parser("report", help="re-render report.html for a run")
    p.add_argument("run_dir", nargs="?")
    mg = sub.add_parser("merge", help="combine shard runs into one report")
    mg.add_argument("out_dir")
    mg.add_argument("shard_dirs", nargs="+")
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "summary":
        print_summary(latest_run(args.run_dir))
    elif args.cmd == "merge":
        cmd_merge(Path(args.out_dir).resolve(), [Path(d).resolve() for d in args.shard_dirs])
    else:
        d = latest_run(args.run_dir)
        write_report(d)
        print(d / "report.html")


if __name__ == "__main__":
    main()
