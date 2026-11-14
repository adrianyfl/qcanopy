#!/usr/bin/env python3
"""Profile MQT QCEC on circuit pairs produced by ``qcanopy_generate``.

Subcommands
-----------
run      Generate datasets for every point of a sweep grid, run ``mqt.qcec.verify``
         on every pair with every configured checker strategy, write
         ``results.jsonl`` + ``run.json``, then render the report.
report   (Re-)render ``report.html`` for an existing run directory.
summary  Print the terminal summary for an existing run directory.

Every check runs in its own forked worker process. This gives a hard timeout
that also stops a stuck native checker, isolation from crashes, and per-check
peak-RSS / CPU-time measurements. See profiling/qcec/README.md.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import itertools
import json
import multiprocessing as mp
import multiprocessing.connection as mpc
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time
import tomllib
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_GENERATOR = REPO / "build" / "src" / "generator" / "qcanopy_generate"
DEFAULT_OUT = REPO / "generated" / "qcec-profile"
TEMPLATE = HERE / "report_template.html"

# Sweep [generator] key -> qcanopy_generate flag.
GEN_FLAGS = {
    "qubits": "--qubits", "gates": "--gates", "pairs": "--pairs",
    "prefix_fraction": "--prefix-fraction", "prefix_group_size": "--prefix-group-size",
    "equivalent_fraction": "--equivalent-fraction", "gate_set": "--gate-set",
    "angle_mode": "--angle-mode", "seed": "--seed", "verify": "--verify",
}
# Sweep strategy options that are pybind enums in mqt.qcec, given as strings in TOML.
ENUM_OPTIONS = {
    "alternating_scheme": "ApplicationScheme", "construction_scheme": "ApplicationScheme",
    "simulation_scheme": "ApplicationScheme", "state_type": "StateType",
}
OUTCOMES = ("correct", "probable", "incorrect", "inconclusive", "timeout", "error", "skipped")


# --------------------------------------------------------------------------- sweep
def load_sweep(path: Path) -> dict[str, Any]:
    with path.open("rb") as f:
        sweep = tomllib.load(f)
    gen = sweep.get("generator", {})
    unknown = set(gen) - set(GEN_FLAGS)
    if unknown:
        sys.exit(f"sweep: unknown [generator] keys: {sorted(unknown)}")
    if not sweep.get("strategies"):
        sweep["strategies"] = {"default": {}}
    sweep.setdefault("name", path.stem)
    sweep.setdefault("repeats", 1)
    sweep.setdefault("timeout", 60.0)
    sweep.setdefault("warmup", True)
    sweep.setdefault("skip_after_timeouts", 0)   # 0 disables auto-skip
    for name, opts in sweep["strategies"].items():
        if not isinstance(opts, dict):
            sys.exit(f"sweep: [strategies.{name}] must be a table")
    budget = sweep.get("pair_gate_budget")
    if budget is not None and (not isinstance(budget, int) or budget <= 0):
        sys.exit("sweep: pair_gate_budget must be a positive integer")
    if "hard_timeout" not in sweep:
        sweep["hard_timeout"] = sweep["timeout"] * 1.5 + 10.0
    return sweep


def grid(gen: dict[str, Any], pair_gate_budget: int | None = None) -> list[dict[str, Any]]:
    """Cartesian product over list-valued generator keys, ordered by gate count.

    With ``pair_gate_budget``, each dataset's pair count is reduced so that
    pairs * gates stays within the budget (minimum 2 pairs, so both labels can
    appear). qcanopy_generate keeps the whole batch in memory (~0.8 KB/gate).
    """
    keys = list(gen)
    axes = [v if isinstance(v, list) else [v] for v in gen.values()]
    configs = [dict(zip(keys, combo)) for combo in itertools.product(*axes)]
    if pair_gate_budget:
        for cfg in configs:
            gates = max(1, int(cfg.get("gates", 100)))
            cfg["pairs"] = max(2, min(int(cfg.get("pairs", 32)), pair_gate_budget // gates))
    return sorted(configs, key=lambda c: c.get("gates", 100))  # stable: keeps TOML order otherwise


def skip_group(cfg: dict[str, Any], strategy: str) -> tuple:
    """Everything that identifies a scaling series except the gate count."""
    return (strategy, *sorted((k, str(v)) for k, v in cfg.items() if k not in ("gates", "pairs")))


# Strategy keys consumed by the profiler instead of being passed to mqt.qcec.verify.
PROFILER_KEYS = ("hard_timeout",)


def strategy_timeouts(sweep: dict[str, Any], options: dict[str, Any]) -> tuple[float, float]:
    """(cooperative, hard) timeout for a strategy.

    A strategy may set its own ``timeout`` (passed to QCEC) and ``hard_timeout``.
    Some QCEC checkers (notably construction) overshoot the cooperative timeout
    by a lot, so the hard kill follows the strategy's own timeout.
    """
    if "timeout" not in options:
        return float(sweep["timeout"]), float(options.get("hard_timeout", sweep["hard_timeout"]))
    coop = float(options["timeout"])
    default_hard = coop * 1.5 + 10.0 if coop > 0 else float(sweep["hard_timeout"])
    return coop, float(options.get("hard_timeout", default_hard))


def dataset_id(cfg: dict[str, Any]) -> str:
    digest = hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:8]
    parts = [f"q{cfg.get('qubits', 4)}", f"g{cfg.get('gates', 100)}", str(cfg.get("gate_set", "mixed"))]
    return "_".join(parts) + "_" + digest


def generate(generator: Path, cfg: dict[str, Any], dest: Path) -> None:
    if (dest / "manifest.jsonl").exists():
        return  # already generated (resume)
    if dest.exists():
        shutil.rmtree(dest)  # incomplete earlier attempt; generator refuses nonempty dirs
    cmd = [str(generator), "--format", "qasm", "--output", str(dest)]
    for k, v in cfg.items():
        cmd += [GEN_FLAGS[k], str(v)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"qcanopy_generate failed for {cfg}:\n{proc.stderr}")


def read_manifest(ds_dir: Path) -> tuple[dict, list[dict]]:
    lines = [json.loads(l) for l in (ds_dir / "manifest.jsonl").read_text().splitlines() if l.strip()]
    return lines[0], lines[1:]


# --------------------------------------------------------------------------- worker
def _proc_status_kb(field: str) -> int | None:
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith(field + ":"):
                return int(line.split()[1])
    except OSError:
        pass
    return None


def _reset_peak_rss() -> bool:
    try:  # Linux >= 4.0: writing 5 resets VmHWM to the current RSS.
        Path("/proc/self/clear_refs").write_text("5")
        return True
    except OSError:
        return False


def _coerce_options(options: dict[str, Any]) -> dict[str, Any]:
    from mqt.qcec import pyqcec
    out = dict(options)
    for key, enum_name in ENUM_OPTIONS.items():
        if isinstance(out.get(key), str):
            out[key] = getattr(pyqcec, enum_name).__members__[out[key]]
    return out


def _decided_by(result_json: dict, equivalence: str) -> str:
    checkers = result_json.get("checkers", [])
    decisive = [c for c in checkers if c.get("equivalence") != "no_information"]
    if not decisive:
        return "preprocessing" if equivalence != "no_information" else "none"
    match = [c for c in decisive if c.get("equivalence") == equivalence]
    return (match or decisive)[0].get("checker", "unknown")


_WARM: tuple | None = None


def _worker(conn, a_path: str, b_path: str, options: dict, timeout: float, warmup: bool) -> None:
    import resource
    from time import perf_counter
    try:
        from mqt import qcec
        from mqt.core import load

        opts = _coerce_options(options)
        opts.setdefault("timeout", timeout)
        if warmup and _WARM is not None:
            # Pay one-time lazy initialisation (allocators, code pages) outside the timing.
            qcec.verify(_WARM[0], _WARM[1], **{k: v for k, v in opts.items() if k != "timeout"})
        precise_rss = _reset_peak_rss()
        rss0 = _proc_status_kb("VmRSS")

        t0 = perf_counter()
        a, b = load(a_path), load(b_path)
        t_load = perf_counter() - t0

        ru0 = resource.getrusage(resource.RUSAGE_SELF)
        t0 = perf_counter()
        res = qcec.verify(a, b, **opts)
        wall = perf_counter() - t0
        ru1 = resource.getrusage(resource.RUSAGE_SELF)

        rj = res.json()
        equivalence = str(rj.get("equivalence", res.equivalence.name))
        peak = _proc_status_kb("VmHWM") if precise_rss else None
        if peak is None:
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        conn.send({
            "status": "ok",
            "equivalence": equivalence,
            "wall_s": wall,
            "load_s": t_load,
            "preprocessing_s": rj.get("preprocessing_time", res.preprocessing_time),
            "check_s": rj.get("check_time", res.check_time),
            "cpu_s": (ru1.ru_utime - ru0.ru_utime) + (ru1.ru_stime - ru0.ru_stime),
            "peak_rss_mb": peak / 1024 if peak else None,
            "rss_delta_mb": (peak - rss0) / 1024 if (peak and rss0 and precise_rss) else None,
            "decided_by": _decided_by(rj, equivalence),
            "simulations_performed": rj.get("simulations", {}).get("performed"),
            "checkers": [{"checker": c.get("checker"), "equivalence": c.get("equivalence"),
                          "runtime": c.get("runtime")} for c in rj.get("checkers", [])
                         if c.get("equivalence") != "no_information" or c.get("runtime")],
        })
    except BaseException as e:  # noqa: BLE001 - report everything to the parent
        try:
            conn.send({"status": "error", "error": f"{type(e).__name__}: {e}"})
        except Exception:
            pass
    finally:
        conn.close()


def _prepare_warmup() -> None:
    """Build the warm-up pair in the parent so forked children inherit it."""
    global _WARM
    from mqt.core import load
    qasm = 'OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\n'
    _WARM = (load(qasm), load(qasm))


# --------------------------------------------------------------------------- classification
def verdict(equivalence: str) -> tuple[str | None, bool]:
    """Map a QCEC EquivalenceCriterion name to (relation, decisive)."""
    if equivalence in ("equivalent", "equivalent_up_to_global_phase", "equivalent_up_to_phase"):
        return "equivalent", True
    if equivalence == "not_equivalent":
        return "non_equivalent", True
    if equivalence == "probably_equivalent":
        return "equivalent", False
    if equivalence == "probably_not_equivalent":
        return "non_equivalent", False
    return None, False


def classify(expected: str, rec: dict, timeout: float) -> str:
    if rec["status"] == "hard_timeout":
        return "timeout"
    if rec["status"] == "skipped":
        return "skipped"
    if rec["status"] != "ok":
        return "error"
    rel, decisive = verdict(rec["equivalence"])
    if rel is None:
        return "timeout" if timeout > 0 and rec["wall_s"] >= 0.98 * timeout else "inconclusive"
    if rel != expected:
        return "incorrect"
    return "correct" if decisive else "probable"


# --------------------------------------------------------------------------- run
def host_info() -> dict[str, Any]:
    cpu = platform.processor() or ""
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    mem_gb = None
    try:
        mem_gb = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30, 1)
    except (ValueError, OSError):
        pass
    return {"hostname": platform.node(), "platform": platform.platform(), "cpu": cpu,
            "logical_cpus": os.cpu_count(), "memory_gb": mem_gb, "python": platform.python_version()}


def git_info() -> dict[str, Any]:
    def git(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True).stdout.strip()
    try:
        return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
                "dirty": bool(git("status", "--porcelain"))}
    except OSError:
        return {}


def cmd_run(args: argparse.Namespace) -> None:
    import mqt.core
    import mqt.qcec

    sweep = load_sweep(args.sweep)
    if args.timeout is not None:
        sweep["timeout"] = args.timeout
        sweep["hard_timeout"] = args.timeout * 1.5 + 10.0
    if args.repeats is not None:
        sweep["repeats"] = args.repeats
    generator = Path(args.generator)
    if not generator.exists():
        sys.exit(f"generator not found: {generator}\nBuild it first: cmake --build build")

    if args.resume:
        run_dir = Path(args.resume)
        if not (run_dir / "run.json").exists():
            sys.exit(f"--resume: {run_dir} has no run.json")
        sweep = json.loads((run_dir / "run.json").read_text())["sweep"]
    else:
        stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        run_dir = Path(args.out) / f"{stamp}-{sweep['name']}"
        run_dir.mkdir(parents=True, exist_ok=False)

    versions = {"mqt.qcec": mqt.qcec.__version__, "mqt.core": mqt.core.__version__}
    if args.resume:
        meta = json.loads((run_dir / "run.json").read_text())
        if meta["versions"] != versions:
            sys.exit(f"--resume: run used {meta['versions']} but this environment has {versions}")
        meta.setdefault("resumed", []).append(dt.datetime.now().isoformat(timespec="seconds"))
    else:
        meta = {
            "schema": 1,
            "sweep": sweep,
            "sweep_file": str(args.sweep),
            "started": dt.datetime.now().isoformat(timespec="seconds"),
            "versions": versions,
            "host": host_info(),
            "git": git_info(),
            "jobs": args.jobs,
            "generator": str(generator),
        }
    (run_dir / "run.json").write_text(json.dumps(meta, indent=2))

    # 1. Generate datasets.
    datasets = []
    configs = grid(sweep.get("generator", {}), sweep.get("pair_gate_budget"))
    print(f"[qcec-profile] run dir: {run_dir}")
    for i, cfg in enumerate(configs, 1):
        ds = dataset_id(cfg)
        print(f"[qcec-profile] generating {i}/{len(configs)}: {ds} ({cfg.get('pairs', 32)} pairs)", flush=True)
        generate(generator, cfg, run_dir / "datasets" / ds)
        _, cases = read_manifest(run_dir / "datasets" / ds)
        datasets.append((ds, cfg, cases))

    # 2. Build task list: repeat-major so drift hits every strategy equally, and
    #    ascending gate count (datasets are sorted) so auto-skip can prune.
    tasks = []
    for rep in range(sweep["repeats"]):
        for ds, cfg, cases in datasets:
            for case in cases:
                for strat, opts in sweep["strategies"].items():
                    tasks.append({"key": f"{ds}/{case['case_id']}/{strat}/{rep}", "dataset": ds,
                                  "cfg": cfg, "case": case, "strategy": strat, "options": opts, "repeat": rep})
    results_path = run_dir / "results.jsonl"
    done = set()
    if results_path.exists():
        done = {json.loads(l)["key"] for l in results_path.read_text().splitlines() if l.strip()}
    pending = [t for t in tasks if t["key"] not in done]
    print(f"[qcec-profile] {len(tasks)} checks ({len(done)} already done), "
          f"{len(sweep['strategies'])} strategies, jobs={args.jobs}", flush=True)

    # 3. Execute.
    if sweep.get("warmup", True):
        _prepare_warmup()
    method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    ctx = mp.get_context(method)
    skip_after = int(sweep.get("skip_after_timeouts", 0))
    timeouts_at: dict[tuple, dict[int, int]] = {}   # skip group -> gates -> timeout count
    skip_above: dict[tuple, int] = {}               # skip group -> skip gate counts above this
    if results_path.exists() and skip_after:
        for l in results_path.read_text().splitlines():
            r = json.loads(l) if l.strip() else None
            if r and r["outcome"] == "timeout":
                task = next(t for t in tasks if t["key"] == r["key"])
                g = skip_group(task["cfg"], task["strategy"])
                timeouts_at.setdefault(g, {}).setdefault(r["gates"], 0)
                timeouts_at[g][r["gates"]] += 1
                if timeouts_at[g][r["gates"]] >= skip_after:
                    skip_above[g] = min(skip_above.get(g, r["gates"]), r["gates"])
    running: dict[Any, tuple] = {}
    n_done, n_total, t_start = 0, len(pending), time.perf_counter()
    bad = 0
    with results_path.open("a") as out:
        def finish(task: dict, payload: dict) -> None:
            nonlocal n_done, bad
            rec = make_record(task, payload, run_dir, strategy_timeouts(sweep, task["options"])[0])
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n_done += 1
            if skip_after and rec["outcome"] == "timeout":
                g = skip_group(task["cfg"], task["strategy"])
                cnt = timeouts_at.setdefault(g, {})
                cnt[rec["gates"]] = cnt.get(rec["gates"], 0) + 1
                if cnt[rec["gates"]] >= skip_after and rec["gates"] < skip_above.get(g, float("inf")):
                    skip_above[g] = rec["gates"]
                    print(f"\n[qcec-profile] {task['strategy']}: {cnt[rec['gates']]} timeouts at "
                          f"{rec['gates']} gates; skipping the rest of this series from {rec['gates']} gates up",
                          flush=True)
            if rec["outcome"] in ("incorrect", "error"):
                bad += 1
                print(f"\n[qcec-profile] {rec['outcome'].upper()}: {rec['key']} "
                      f"expected={rec['expected']} got={rec.get('equivalence')} {rec.get('error') or ''}")
            if sys.stdout.isatty() or n_done % 25 == 0 or n_done == n_total:
                el = time.perf_counter() - t_start
                eta = el / n_done * (n_total - n_done)
                print(f"\r[qcec-profile] {n_done}/{n_total}  elapsed {el:6.1f}s  eta {eta:6.1f}s  "
                      f"flagged {bad}   ", end="" if sys.stdout.isatty() else "\n", flush=True)

        queue = list(pending)
        while queue or running:
            while queue and len(running) < args.jobs:
                task = queue.pop(0)
                limit = skip_above.get(skip_group(task["cfg"], task["strategy"]))
                if limit is not None and task["cfg"].get("gates", 100) >= limit:
                    finish(task, {"status": "skipped", "error": f"timed out at {limit} gates"})
                    continue
                ds_dir = run_dir / "datasets" / task["dataset"]
                timeout, hard = strategy_timeouts(sweep, task["options"])
                parent, child = ctx.Pipe(duplex=False)
                p = ctx.Process(target=_worker, args=(
                    child, str(ds_dir / task["case"]["files"]["a"]), str(ds_dir / task["case"]["files"]["b"]),
                    {k: v for k, v in task["options"].items() if k not in PROFILER_KEYS},
                    timeout, sweep.get("warmup", True)), daemon=True)
                p.start()
                child.close()
                running[parent] = (p, task, time.perf_counter(), hard)
            ready = mpc.wait(list(running), timeout=0.25)
            for conn in ready:
                p, task, _, _ = running.pop(conn)
                try:
                    payload = conn.recv()
                except EOFError:
                    p.join(5)
                    payload = {"status": "error", "error": f"worker died (exit code {p.exitcode})"}
                conn.close()
                p.join(5)
                finish(task, payload)
            now = time.perf_counter()
            for conn, (p, task, t0, hard) in list(running.items()):
                if now - t0 > hard:
                    p.kill()
                    p.join(5)
                    conn.close()
                    running.pop(conn)
                    finish(task, {"status": "hard_timeout", "error": f"killed after {hard:.0f}s", "wall_s": now - t0})
    print()
    meta["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    (run_dir / "run.json").write_text(json.dumps(meta, indent=2))
    print_summary(run_dir)
    path = render_report(run_dir)
    print(f"\n[qcec-profile] report: {path}")


def make_record(task: dict, payload: dict, run_dir: Path, timeout: float) -> dict:
    cfg, case = task["cfg"], task["case"]
    sharing = case.get("sharing", {})
    rel = lambda p: str(Path("datasets") / task["dataset"] / p)  # noqa: E731
    rec = {
        "key": task["key"], "dataset": task["dataset"], "case_id": case["case_id"],
        "strategy": task["strategy"], "repeat": task["repeat"],
        "qubits": cfg.get("qubits", 4), "gates": cfg.get("gates", 100),
        "gate_set": cfg.get("gate_set", "mixed"), "angle_mode": cfg.get("angle_mode", "mixed"),
        "prefix_fraction": cfg.get("prefix_fraction", 0.5),
        "expected": case["expected_relation"],
        "transform": case.get("transform", {}).get("kind"),
        "a_gates": sharing.get("raw", {}).get("a_gates"), "b_gates": sharing.get("raw", {}).get("b_gates"),
        "raw_prefix": sharing.get("raw", {}).get("common_prefix_fraction"),
        "norm_prefix": sharing.get("normalized", {}).get("common_prefix_fraction"),
        "norm_a_gates": sharing.get("normalized", {}).get("a_gates"),
        "norm_b_gates": sharing.get("normalized", {}).get("b_gates"),
        "files": {"a": rel(case["files"]["a"]), "b": rel(case["files"]["b"])},
        **payload,
    }
    rec["outcome"] = classify(rec["expected"], rec, timeout)
    return rec


# --------------------------------------------------------------------------- summary
def load_run(run_dir: Path) -> tuple[dict, list[dict]]:
    meta = json.loads((run_dir / "run.json").read_text())
    path = run_dir / "results.jsonl"
    recs = [json.loads(l) for l in path.read_text().splitlines() if l.strip()] if path.exists() else []
    return meta, recs


def _q(xs: list[float], p: float) -> float:
    xs = sorted(xs)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    lo, hi = int(k), min(int(k) + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def fmt_t(s: float | None) -> str:
    if s is None or s != s:
        return "-"
    if s < 1e-3:
        return f"{s * 1e6:.0f}µs"
    if s < 1:
        return f"{s * 1e3:.1f}ms"
    return f"{s:.2f}s"


def print_summary(run_dir: Path) -> None:
    meta, recs = load_run(run_dir)
    if not recs:
        print("no results yet")
        return
    strategies = list(meta["sweep"]["strategies"])
    v = meta["versions"]
    print(f"\n== MQT QCEC profile: {meta['sweep']['name']}  (qcec {v['mqt.qcec']}, mqt.core {v['mqt.core']}) ==")
    sw = max(10, *(len(s) + 2 for s in strategies))
    hdr = f"{'strategy':<{sw}}{'checks':>7}{'correct':>8}{'prob':>6}{'WRONG':>6}{'incon':>6}{'t/o':>5}{'err':>5}{'skip':>5}" \
          f"{'median':>10}{'p90':>10}{'max':>10}{'total':>10}"
    print(hdr)
    print("-" * len(hdr))
    for s in strategies:
        rs = [r for r in recs if r["strategy"] == s]
        if not rs:
            continue
        c = {o: sum(r["outcome"] == o for r in rs) for o in OUTCOMES}
        ts = [r["wall_s"] for r in rs if r.get("wall_s") is not None]
        print(f"{s:<{sw}}{len(rs):>7}{c['correct']:>8}{c['probable']:>6}{c['incorrect']:>6}{c['inconclusive']:>6}"
              f"{c['timeout']:>5}{c['error']:>5}{c['skipped']:>5}{fmt_t(_q(ts, .5)):>10}{fmt_t(_q(ts, .9)):>10}"
              f"{fmt_t(max(ts) if ts else None):>10}{fmt_t(sum(ts)):>10}")

    # Median wall time per (qubits, gates) x strategy.
    points = sorted({(r["qubits"], r["gates"]) for r in recs})
    if len(points) > 1:
        print("\nmedian verify time by size")
        print(f"{'qubits':>7}{'gates':>7}" + "".join(f"{s[:13]:>14}" for s in strategies))
        for q, g in points:
            row = f"{q:>7}{g:>7}"
            for s in strategies:
                ts = [r["wall_s"] for r in recs if r["strategy"] == s and r["qubits"] == q
                      and r["gates"] == g and r.get("wall_s") is not None]
                row += f"{fmt_t(statistics.median(ts)) if ts else '-':>14}"
            print(row)

    flagged = [r for r in recs if r["outcome"] in ("incorrect", "error", "timeout")]
    if flagged:
        print(f"\nflagged checks ({len(flagged)}):")
        for r in flagged[:20]:
            print(f"  {r['outcome']:<10} {r['key']}  expected={r['expected']} got={r.get('equivalence')} "
                  f"{r.get('error') or ''}")
        if len(flagged) > 20:
            print(f"  ... {len(flagged) - 20} more (see report.html)")


# --------------------------------------------------------------------------- report
REPORT_FIELDS = ("dataset", "case_id", "strategy", "repeat", "qubits", "gates", "gate_set", "angle_mode",
                 "prefix_fraction", "expected", "transform", "a_gates", "b_gates", "norm_a_gates", "norm_b_gates",
                 "raw_prefix", "norm_prefix",
                 "files", "status", "equivalence", "outcome", "wall_s", "load_s", "preprocessing_s", "check_s",
                 "cpu_s", "peak_rss_mb", "rss_delta_mb", "decided_by", "simulations_performed", "error")


def render_report(run_dir: Path) -> Path:
    meta, recs = load_run(run_dir)
    data = {"meta": meta, "records": [{k: r.get(k) for k in REPORT_FIELDS} for r in recs]}
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html = TEMPLATE.read_text().replace("/*__DATA__*/null", blob)
    out = run_dir / "report.html"
    out.write_text(html)
    return out


def latest_run(out: Path) -> Path:
    runs = sorted(p for p in out.glob("*") if (p / "run.json").exists())
    if not runs:
        sys.exit(f"no runs under {out}")
    return runs[-1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="generate datasets and profile QCEC")
    r.add_argument("--sweep", type=Path, default=HERE / "sweeps" / "smoke.toml")
    r.add_argument("--generator", default=str(DEFAULT_GENERATOR))
    r.add_argument("--out", default=str(DEFAULT_OUT), help="parent directory for run directories")
    r.add_argument("--jobs", type=int, default=1, help="concurrent checks (keep 1 for clean timings)")
    r.add_argument("--timeout", type=float, help="override sweep timeout (seconds)")
    r.add_argument("--repeats", type=int, help="override sweep repeats")
    r.add_argument("--resume", help="continue an interrupted run directory")
    for name, hlp in (("report", "render report.html"), ("summary", "print terminal summary")):
        s = sub.add_parser(name, help=hlp)
        s.add_argument("run_dir", nargs="?", help="run directory (default: latest)")
        s.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
        return
    run_dir = Path(args.run_dir) if args.run_dir else latest_run(Path(args.out))
    if args.cmd == "report":
        print(render_report(run_dir))
    else:
        print_summary(run_dir)


if __name__ == "__main__":
    main()
