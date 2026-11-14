# MQT QCEC profiler

Profiles [MQT QCEC](https://github.com/munich-quantum-toolkit/qcec) on labeled
circuit pairs from `qcanopy_generate`. One command generates datasets across a
size grid, runs `mqt.qcec.verify` on every pair with one or more checker
configurations, checks each verdict against the generator's proven label, and
writes a terminal summary plus a self-contained HTML report.

## Setup

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
python3 -m venv .venv
.venv/bin/pip install -r profiling/qcec/requirements.txt
```

The profiler needs Python 3.11+ (it uses `tomllib`) and a Linux host for exact
per-check peak-RSS numbers. Elsewhere it falls back to `ru_maxrss`.

## Run

```sh
.venv/bin/python profiling/qcec/qcec_profile.py run                       # sweeps/smoke.toml, ~1-2 min
.venv/bin/python profiling/qcec/qcec_profile.py run --sweep profiling/qcec/sweeps/scaling.toml
.venv/bin/python profiling/qcec/qcec_profile.py run --sweep profiling/qcec/sweeps/strategies.toml
```

Each run creates `generated/qcec-profile/<timestamp>-<sweep name>/`:

```text
run.json         sweep, qcec/mqt.core versions, host, git commit, timings
results.jsonl    one record per check (written as it goes)
datasets/<id>/   the qcanopy_generate QASM output for each grid point
report.html      the report; open it in a browser
```

Other commands:

```sh
qcec_profile.py summary [RUN_DIR]     # print the terminal tables again (default: latest run)
qcec_profile.py report  [RUN_DIR]     # re-render report.html (e.g. after editing the template)
qcec_profile.py run --resume RUN_DIR  # finish an interrupted run, skipping completed checks
```

`run` options: `--jobs N` (concurrent checks, default 1), `--timeout S`, `--repeats N`,
`--out DIR`, `--generator PATH`.

## Sweep files

A sweep is a TOML file. Scalar `[generator]` values apply to every dataset,
and list values form a cartesian product. Keys map one-to-one to
`qcanopy_generate` options. Each `[strategies.<name>]` table is passed to
`mqt.qcec.verify` as keyword arguments (any `ConfigurationOptions` key works;
enum options like `state_type = "stabilizer"` or `alternating_scheme = "lookahead"` take the member name as a
string). An empty table means QCEC defaults. A strategy may set its own
`timeout` (passed to QCEC) and `hard_timeout` (used by the profiler only).

```toml
name = "my-sweep"
repeats = 3            # every check runs this many times; the report uses the median
timeout = 60.0         # QCEC's cooperative timeout, per check
hard_timeout = 100.0   # the worker is killed after this (default 1.5 * timeout + 10)
warmup = true          # untimed tiny check in each worker before the timed one
pair_gate_budget = 4_000_000  # optional: pairs = min(pairs, budget // gates), at least 2
skip_after_timeouts = 2       # optional: after 2 timeouts at N gates, skip that series from N up

[generator]
qubits = [4, 8, 16]
gates = [100, 400]
gate_set = ["clifford-t", "mixed"]
pairs = 16
seed = 759

[strategies.default]
[strategies.alternating]
run_alternating_checker = true
run_construction_checker = false
run_simulation_checker = false
run_zx_checker = false
```

Datasets run in ascending gate count. `skip_after_timeouts` uses that order:
once a strategy reaches the threshold at N gates for a given qubit count, gate
set, and other settings, its remaining checks at N gates and above are recorded
as `skipped` instead of each waiting for the timeout. Skipping lasts across repeats.

Presets:

| Sweep | Purpose |
| --- | --- |
| `smoke` | Quick sanity check of every checker at small and medium sizes |
| `scaling` | Qubits × gates growth for the default, alternating, and construction checkers |
| `strategies` | Every checker at 10 qubits across gate pools and angle modes |
| `large_gates` | 1k to 100k base gates (half-decade steps) at 8, 16, and 32 qubits |

## Large gate counts

`large_gates.toml` is set up for very long circuits. Resource limits to know:

- **Generator memory.** `qcanopy_generate` keeps the whole batch in memory, about
  0.8 KB per gate, so 4 pairs at 1M gates peak around 3.3 GB. `pair_gate_budget`
  reduces the pair count automatically for the largest datasets.
- **Disk.** A circuit takes about 38 MB of QASM per 1M gates; the `large_gates` preset (up to 100k) writes about 250 MB.
- **Parse time.** `mqt.core.load` takes about 3 s per 1M-gate file. Every check
  re-parses its two files; this is recorded as `load_s` and is not in `wall_s`.
- **Checker cliffs.** On random non-Clifford circuits, the simulation-only and
  construction checkers time out from about 1k gates at 16 qubits, and ZX from
  about 16k–100k. These strategies get a 60 s timeout, and auto-skip stops them
  from consuming the run.
- **Timeout overshoot.** QCEC's construction checker can run well past its
  cooperative timeout, so the profiler's hard kill follows each strategy's own timeout.

Choose **Verify time per gate** in the report's metric selector. A flat line
means linear scaling in gate count; a rising line means superlinear growth.

## What is measured

Every check runs in a fresh process forked from the profiler, so one slow or
crashing check cannot affect the next one, and the hard timeout can stop native
code. In each worker:

| Field | Meaning |
| --- | --- |
| `wall_s` | `perf_counter` around `verify()` only, without circuit parsing |
| per gate (report only) | `wall_s / (gates in A + gates in B)` |
| `load_s` | `mqt.core.load` of both QASM files |
| `preprocessing_s`, `check_s` | QCEC's own timers from `Results.json()` |
| `cpu_s` | user + system CPU of the worker during `verify()`, all threads. `cpu_s / wall_s` above 1 means checkers ran in parallel |
| `peak_rss_mb`, `rss_delta_mb` | the worker's peak RSS (`VmHWM`, reset just before loading), and its increase over the pre-load RSS |
| `decided_by` | the checker whose result became the final verdict, or `preprocessing` if QCEC decided before any checker ran |

Outcomes compare QCEC's verdict with the label the generator proves:

| Outcome | Meaning |
| --- | --- |
| `correct` | Decisive verdict that matches the label. `equivalent_up_to_global_phase` counts as equivalent, because labels allow one global phase |
| `probable` | `probably_equivalent` / `probably_not_equivalent` that matches the label. Expected for simulation-only and ZX-only runs |
| `incorrect` | The verdict contradicts the label. Investigate: this points to a checker bug or a QASM import problem |
| `inconclusive` | `no_information` before the timeout |
| `timeout` | `no_information` at the cooperative timeout, or the worker was killed |
| `error` | An exception or a worker crash; the message is in `error` |
| `skipped` | Not run, because auto-skip had already stopped this series |

## Interpreting results

- Keep `--jobs 1` for timing. Concurrent checks compete for cores and memory
  bandwidth, and QCEC's default configuration already runs checkers in
  parallel threads.
- Small checks (well under 1 ms) are dominated by fixed overhead. Compare them
  through `check_s` rather than `wall_s`.
- The generator's equivalent pairs differ by only one local rewrite, so the
  alternating checker is expected to do well. These datasets measure QCEC on
  near-identical pairs with shared prefixes, not on arbitrary compiled circuits.
- `results.jsonl` is plain JSON lines. Load it into pandas for custom analysis:
  `pd.read_json(path, lines=True)`.
