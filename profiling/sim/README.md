# Simulator profiler

Profiles the state-vector simulator (`src/sim/sim.cpp`, `qcanopy::Simulator`) on
labeled circuit pairs from `qcanopy_generate`. One command generates datasets
across a size grid, runs every pair through the simulator, checks each result
against the generator's proven label
and an independent reference simulation, and writes a terminal summary plus a
self-contained HTML report.

## Setup

1. **Build in Release** (timings from a Debug build are meaningless):

   ```sh
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ```

   With the workspace toolchain from `docs/generator.md` (PowerShell):

   ```powershell
   .\.tools\cmake\data\bin\cmake.exe -S . -B build-llvm -G Ninja -DCMAKE_BUILD_TYPE=Release
   .\.tools\cmake\data\bin\cmake.exe --build build-llvm
   ```

2. **Check it.** Python 3.11+ is the only other requirement (standard library
   only, no `pip install`):

   ```sh
   python profiling/sim/sim_profile.py run      # smoke sweep, under a minute
   ```

## Sweep files

A sweep is a TOML file. Scalar `[generator]` values apply to every dataset, and
list values form a cartesian product. Keys map one-to-one to `qcanopy_generate`
options (the profiler always uses `--format ir`). Each circuit is timed with
`Simulator::execute(const Circuit&)` as generated.

```toml
name = "my-sweep"
repeats = 3                  # timed executions of each circuit; the median is reported
timeout = 120.0              # the worker is killed after this many seconds
probes = 3                   # |0...0> plus 2 random product-state inputs per pair
tol = 1e-8                   # fidelity tolerance for the equivalence verdict
reference_max_qubits = 12    # independent reference simulation up to this size
reference_tol = 1e-8         # max |amplitude difference| allowed against the reference
memory_budget_mb = 8192      # default: half of RAM
kernel_scan = true           # or a list of qubit counts, or false
pair_gate_budget = 400_000   # optional: pairs = min(pairs, budget // gates), at least 2
skip_after_timeouts = 2      # optional: after 2 timeouts at N gates, skip that series from N up
warmup = true                # untimed tiny circuit in each worker before timing

[generator]
qubits = [8, 16, 20]
gates = [100, 400]
gate_set = ["clifford-t", "mixed"]
pairs = 8
seed = 759
```

Datasets run in ascending gate count, then qubit count. `skip_after_timeouts`
works as in the QCEC profiler: once the threshold is reached at N gates for a
given qubit count, gate set and other settings, the remaining checks at N gates
and above are recorded as `skipped`.

Presets:

| Sweep | Purpose |
| --- | --- |
| `smoke` | 3 to 14 qubits, 50 and 200 gates |
| `scaling` | 4 to 24 qubits × 100 and 400 gates, with kernel scans |
| `large_gates` | 1k to 100k base gates (half-decade steps) at 8, 12 and 16 qubits |

## Changing parameters

Every value comes from the sweep file. To change one, either copy a preset and
edit it, or override it on the command line with `--set KEY=VALUE` (repeatable;
the value is TOML, so lists use brackets and strings may be bare):

```sh
python profiling/sim/sim_profile.py run --sweep profiling/sim/sweeps/scaling.toml \
    --set 'generator.qubits=[24,26,28]' --set 'generator.gates=[50,200]' \
    --set generator.pairs=2 --set repeats=1 --set probes=2 \
    --set timeout=7200 --set memory_budget_mb=120000 --set name=big
```

Dotted keys reach into tables (`generator.gate_set=clifford-t`). `--timeout`
and `--repeats` are shortcuts for the top-level keys. Overrides are recorded in `run.json` and a
resumed run keeps them.

| To increase | Key | Limit and cost |
| --- | --- | --- |
| Qubits | `generator.qubits` | Generator maximum 32. Memory is 16·2^n bytes per state and time per gate doubles with each qubit. Raise `memory_budget_mb` too, or larger sizes are recorded as `skipped_memory` |
| Gates | `generator.gates` | Linear in time. The generator keeps the batch in memory (~0.8 KB per gate); use `pair_gate_budget` to cap pairs × gates |
| Pairs per grid point | `generator.pairs` | Linear in time; more pairs tighten the medians |
| Grid points | any `generator` key as a list | Cartesian product; every combination is a dataset |
| Timing precision | `repeats` | Each circuit is executed `repeats` times; at large n use 1 |
| Equivalence coverage | `probes` | One extra execution of A and B per probe beyond the first |
| Reference check size | `reference_max_qubits` | Cost ~4·gates·2^n per circuit; also capped at 4e8 operations by the bench |
| Time per check | `timeout` | Hard kill; combine with `skip_after_timeouts` |
| Memory allowed | `memory_budget_mb` | Default: half of RAM |
| Kernel scans | `kernel_scan` | `true` (every qubit count), a list such as `[16, 20, 24]`, or `false` |

Memory and time per check at larger sizes (one core, ~3 ns per amplitude per gate):

| Qubits | One state | Budget needed | One 100-gate execution |
| --- | --- | --- | --- |
| 24 | 256 MB | 0.8 GB | ~5 s |
| 26 | 1 GB | 3.1 GB | ~20 s |
| 28 | 4 GB | 12 GB | ~1.5 min |
| 30 | 16 GB | 48 GB | ~5 min |
| 32 | 64 GB | 192 GB | ~20 min |

A check executes each circuit `repeats` times plus once per extra probe, so at
30 qubits with `repeats = 1`, `probes = 2` and 100 gates, one pair is about
4 × 5 min. Keep `pairs` small and `gates` modest there.

## What is measured

Every check runs `sim_bench` in a fresh process, so one slow or crashing check
cannot affect the next, the hard timeout can stop native code, and peak memory
is per check.

| Field | Meaning |
| --- | --- |
| `sim_s` | per circuit: mean of the A and B medians over `repeats`, input `\|0...0>` |
| `sim_a_s`, `sim_b_s` | every timed repetition |
| `ns_per_gate`, `ns_per_gate_amp` | `sim_s` divided by gates, and by gates × 2^n |
| `gb_s` | effective bandwidth: 32 bytes (one complex read and write) per amplitude per gate |
| `load_s` | `read_ir` of both files |
| `gates_a/b` | gates in each circuit |
| `peak_rss_mb` | peak RSS of the worker: `VmHWM` on Linux (reset after warm-up), `PeakWorkingSetSize` on Windows, `ru_maxrss` on macOS. On Windows and macOS it includes the warm-up and cannot be reset |
| `min_fidelity`, `phase_spread` | minimum `\|<A psi\|B psi>\|` over probes, and the spread of its phase |
| `ref_err_a/b` | max amplitude difference against the reference simulation (null above the size cap) |
| `process_s` | whole worker process, as seen by the runner |

Outcomes:

| Outcome | Meaning |
| --- | --- |
| `correct` | The probe verdict matches the label, and both states match the reference when it ran |
| `incorrect` | An equivalent pair gave different states, a state's norm drifted by more than 1e-6, or a state differs from the reference. This is a simulator bug |
| `undetected` | A non-equivalent pair whose probe states all matched. Expected only with `probes = 1` (the generator's appended X can leave `\|0...0>` unchanged up to phase) |
| `timeout` | The worker was killed at `timeout` |
| `error` | An exception, allocation failure or crash; the message is in `error`. Retried by `--resume` |
| `skipped` | Not run, because auto-skip had already stopped this series |
| `skipped_memory` | Not run, because the state would exceed `memory_budget_mb` |

### Kernel scan

For each qubit count, the bench times circuits made of one repeated gate and
subtracts an empty circuit (state allocation and initialization): `h` and `z`
on every target, `cx`, `cz` and `swap` with control 0 and every target. sim.cpp
has no gate-specific paths, so diagonal (`z`, `cz`) and permutation (`swap`)
gates currently cost the same as dense ones. The scan shows what a diagonal or
permutation fast path would save, and a rise at high targets shows the cost of
the 2^q stride.

## Interpreting results

- Keep `--jobs 1` for timing. Concurrent workers compete for memory bandwidth,
  which is the simulator's main cost above a few qubits.
- Below about 10 qubits a circuit runs in microseconds and timer, allocation
  and loop overhead dominate; compare those sizes through the kernel scan.
- `ns_per_gate_amp` should be roughly flat in gates at fixed n. Steps as n
  grows mark the state leaving L2, L3 and then fitting only in DRAM (about
  16, 19–22 qubits on typical desktop parts).
- The generator's pairs are near-identical circuits with shared prefixes. They
  are good for checking correctness on every gate kind and for scaling curves,
  not a model of compiled application circuits.
- `results.jsonl` is plain JSON lines: `pd.read_json(path, lines=True)`.