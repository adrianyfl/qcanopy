# Simulator profiler

Profiles the state-vector simulator (`src/sim/sim.cpp`, `qcanopy::Simulator`) on
labeled circuit pairs from `qcanopy_generate`. One command generates datasets
across a size grid, runs every pair through the simulator with one or more
execution strategies, checks each result against the generator's proven label
and an independent reference simulation, and writes a terminal summary plus a
self-contained HTML report.

It follows the layout of `profiling/qcec` (same sweep format, run directory,
`run`/`summary`/`report`/`--resume` commands and outcome table), with the
checker replaced by the simulator and the metrics chosen for a state-vector
code: time per gate per amplitude, effective memory bandwidth, and a per-kernel
scan.

## Setup

1. **Put the files in the repository.** Unzip so that this folder is
   `profiling/sim/` next to `profiling/qcec/`. The simulator (`src/sim/`,
   `include/qcanopy/sim.hpp`) and the generator (`src/generator/`) must be on the
   same branch, so merge the `sim` branch into the branch that has the generator
   (or the other way round) first.
2. **Register the bench with CMake.** In the top-level `CMakeLists.txt`, after
   the existing `add_subdirectory` lines for the IR, DAG and generator:

   ```cmake
   add_subdirectory(profiling/sim)
   ```

   `profiling/sim/CMakeLists.txt` links `qcanopy::sim` if that target exists and
   otherwise compiles `src/sim/sim.cpp` into the bench. If sim.cpp lives
   elsewhere, configure with `-DQCANOPY_SIM_SOURCE=path/to/sim.cpp`.
3. **Build in Release** (timings from a Debug build are meaningless):

   ```sh
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ```

   With the workspace toolchain from `docs/generator.md` (PowerShell):

   ```powershell
   .\.tools\cmake\data\bin\cmake.exe -S . -B build-llvm -G Ninja -DCMAKE_BUILD_TYPE=Release
   .\.tools\cmake\data\bin\cmake.exe --build build-llvm
   ```

4. **Check it.** Python 3.11+ is the only other requirement (standard library
   only, no `pip install`):

   ```sh
   python profiling/sim/sim_profile.py run      # smoke sweep, under a minute
   ```

   The runner finds `qcanopy_generate` and `sim_bench` under any `build*/`
   directory (newest wins); pass `--generator` / `--bench` to choose explicitly.
   If the build fails in `sim_bench_adapter.hpp`, or the runner stops in
   `read_manifest()`, see [Adapting to the repository](#adapting-to-the-repository).
5. Add `generated/` to `.gitignore` if it is not there already.

## Run

```sh
python profiling/sim/sim_profile.py run                                   # sweeps/smoke.toml, under a minute
python profiling/sim/sim_profile.py run --sweep profiling/sim/sweeps/scaling.toml
python profiling/sim/sim_profile.py run --sweep profiling/sim/sweeps/strategies.toml
```

The runner finds `qcanopy_generate` and `sim_bench` under any `build*/`
directory (newest wins); pass `--generator` / `--bench` to choose explicitly.

Each run creates `generated/sim-profile/<timestamp>-<sweep name>/`:

```text
run.json         sweep, host, git commit, sim.cpp hash, timings
results.jsonl    one record per check (written as it goes)
kernels.jsonl    one kernel scan per qubit count
datasets/<id>/   the qcanopy_generate IR output for each grid point
report.html      the report; open it in a browser
```

Other commands:

```sh
sim_profile.py summary [RUN_DIR]     # print the terminal tables again (default: latest run)
sim_profile.py report  [RUN_DIR]     # re-render report.html (e.g. after editing the template)
sim_profile.py run --resume RUN_DIR  # finish an interrupted run; completed checks are kept, errors retried
sim_profile.py merge OUT_DIR SHARD_DIR...   # combine shard runs into one report
```

`run` options: `--set KEY=VALUE` (override any sweep value, see below),
`--jobs N` (concurrent checks, default 1), `--timeout S`, `--repeats N`,
`--out DIR` or `--run-dir DIR`, `--shard K/N`, `--generator PATH`, `--bench PATH`.
Ctrl+C or SIGTERM stops the running worker and leaves the run resumable.

`run.json` records a short SHA-256 of `src/sim/sim.cpp` next to the commit, so
two reports from before and after an optimization can be told apart even with
uncommitted edits.

## Sweep files

A sweep is a TOML file. Scalar `[generator]` values apply to every dataset, and
list values form a cartesian product. Keys map one-to-one to `qcanopy_generate`
options (the profiler always uses `--format ir`). Each `[strategies.<name>]`
table chooses how circuits reach the simulator:

| `strategy` | What is timed |
| --- | --- |
| `flat` | `Simulator::execute(const Circuit&)` on the raw circuit |
| `dag` | `Simulator::execute(const CircuitDAG&)`, including its `to_circuit()` flattening. Building the DAG is `prep_s` |
| `normalized` | `execute` on `normalized_copy(circuit)`. Normalization is `prep_s`; fewer gates may be simulated |

A table named `flat`, `dag` or `normalized` needs no `strategy` key. A strategy
may override `timeout`, `repeats`, `probes` and `reference_max_qubits`.

```toml
name = "my-sweep"
repeats = 3                  # timed executions of each circuit; the median is reported
timeout = 120.0              # the worker is killed after this many seconds
probes = 3                   # |0...0> plus 2 random product-state inputs per pair
tol = 1e-8                   # fidelity tolerance for the equivalence verdict
reference_max_qubits = 12    # independent reference simulation up to this size
reference_tol = 1e-8         # max |amplitude difference| allowed against the reference
memory_budget_mb = 8192      # default: half of RAM, or 85% of the Slurm allocation
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

[strategies.flat]
[strategies.normalized]
```

Datasets run in ascending gate count, then qubit count. `skip_after_timeouts`
works as in the QCEC profiler: once a strategy reaches the threshold at N gates
for a given qubit count, gate set and other settings, its remaining checks at
N gates and above are recorded as `skipped`.

Presets:

| Sweep | Purpose |
| --- | --- |
| `smoke` | Every strategy at 3 to 14 qubits, 50 and 200 gates |
| `scaling` | 4 to 24 qubits × 100 and 400 gates, flat only, with kernel scans |
| `strategies` | flat vs dag vs normalized at 12 and 18 qubits across gate pools and angle modes |
| `large_gates` | 1k to 100k base gates (half-decade steps) at 8, 12 and 16 qubits |

## Memory

A state vector takes 16 · 2^n bytes: 16 MB at 20 qubits, 256 MB at 24, 4 GB at
28, 64 GB at 32. A check holds two states (A and B) plus a temporary, so the
runner estimates 3 · 16 · 2^n bytes and records datasets above
`memory_budget_mb` as `skipped_memory` without generating them. Raise the
budget deliberately; the bench has no out-of-core fallback and an allocation
failure is reported as an `error`.

Time scales the same way: one gate costs roughly 2^n × (2–5 ns) on one core,
so 400 gates at 24 qubits is about 30 s per execution, and a check runs each
circuit `repeats` times plus once per extra probe.

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

Dotted keys reach into tables (`strategies.flat.timeout=3600`); a new
`--set strategies.dag={}` adds a strategy. `--timeout` and `--repeats` are
shortcuts for the top-level keys. Overrides are recorded in `run.json` and a
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
| Time per check | `timeout`, or `strategies.<name>.timeout` | Hard kill; combine with `skip_after_timeouts` |
| Memory allowed | `memory_budget_mb` | Default: half of RAM, or 85% of the Slurm allocation inside a job |
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

## Running on a Slurm cluster

The scripts in `profiling/sim/slurm/` run from the repository root.

1. **Environment.** Edit `slurm/env.sh` for your site: load a C++17 compiler,
   CMake ≥ 3.18 and Python ≥ 3.11 (modules or conda). Every job sources it.
2. **Build on a compute node**, so the binaries match the nodes that run them:

   ```sh
   sbatch profiling/sim/slurm/build.sbatch          # builds into build-slurm/
   ```

   For `-march=native`, add `CMAKE_ARGS="-DCMAKE_CXX_FLAGS=-march=native"`, and
   only if build and profiling jobs land on the same CPU type (pin both with
   `-p`/`--constraint`).
3. **One sweep in one job:**

   ```sh
   sbatch profiling/sim/slurm/profile.sbatch                              # scaling sweep
   SWEEP=profiling/sim/sweeps/strategies.toml sbatch profiling/sim/slurm/profile.sbatch
   sbatch --mem=200G --time=24:00:00 -p bigmem profiling/sim/slurm/profile.sbatch \
       --set 'generator.qubits=[28,30]' --set generator.pairs=2 --set repeats=1
   ```

   The run goes to `generated/sim-profile/<sweep>-slurm/` (set `RUN_DIR` to
   choose). Five minutes before the time limit Slurm signals the job; the runner
   stops its worker and exits cleanly. **Resubmit the same command to resume**:
   completed checks are kept. Set `RUN_DIR` to start a separate run.
4. **Split a sweep over several nodes** with an array job and an automatic merge:

   ```sh
   profiling/sim/slurm/submit_array.sh profiling/sim/sweeps/scaling.toml 4
   SBATCH_ARGS="--mem=200G --time=12:00:00 -p bigmem" \
     profiling/sim/slurm/submit_array.sh profiling/sim/sweeps/scaling.toml 4 \
     --set 'generator.qubits=[24,26,28,30]'
   ```

   Each task runs `sim_profile.py run --shard K/N` into `<base>/shard-K/`, and a
   dependent job merges them into `<base>/merged/report.html`. Shards are whole
   series (all gate counts of one qubit count / gate set), dealt round-robin by
   qubit count, so auto-skip keeps working and the largest sizes spread across
   tasks. To finish shards that hit the time limit, resubmit only those IDs
   (the script prints the command), then rerun `merge.sbatch`. `merge` also
   works by hand: `sim_profile.py merge OUT_DIR SHARD_DIR...`.

Resources:

- `--mem` sets the largest qubit count: the runner reads the allocation from
  `SLURM_MEM_PER_NODE` (or `SLURM_MEM_PER_CPU` × CPUs) and uses 85% of it as
  `memory_budget_mb`; use the table above. Sizes that don't fit are recorded as
  `skipped_memory`, not killed by the OOM handler.
- `--cpus-per-task=2` is enough: the simulator is single-threaded. Keep the
  default `--jobs 1`.
- `--exclusive` (set in the scripts) keeps other jobs off the node's memory
  bandwidth, which is what the simulator is limited by. Drop it if the queue
  wait is long and accept noisier timings.
- Compare timings only within one CPU model. `run.json` and the report record
  the node, CPU and Slurm job, and `merge` warns if shards ran on different CPU
  models or different sim.cpp versions.

## What is measured

Every check runs `sim_bench` in a fresh process, so one slow or crashing check
cannot affect the next, the hard timeout can stop native code, and peak memory
is per check.

| Field | Meaning |
| --- | --- |
| `sim_s` | per circuit: mean of the A and B medians over `repeats`, input `\|0...0>` |
| `sim_a_s`, `sim_b_s` | every timed repetition |
| `ns_per_gate`, `ns_per_gate_amp` | `sim_s` divided by simulated gates, and by gates × 2^n |
| `gb_s` | effective bandwidth: 32 bytes (one complex read and write) per amplitude per gate |
| `prep_s` | DAG build or normalization time (null for `flat`) |
| `load_s` | `read_ir` of both files |
| `simulated_gates_a/b` | gates actually applied (differs from `gates_a/b` after normalization) |
| `peak_rss_mb` | peak RSS of the worker: `VmHWM` on Linux (reset after warm-up), `PeakWorkingSetSize` on Windows, `ru_maxrss` on macOS. On Windows and macOS it includes the warm-up and cannot be reset |
| `min_fidelity`, `phase_spread` | minimum `\|<A psi\|B psi>\|` over probes, and the spread of its phase |
| `ref_err_a/b` | max amplitude difference against the reference simulation (null above the size cap) |
| `process_s` | whole worker process, as seen by the runner |

Outcomes:

| Outcome | Meaning |
| --- | --- |
| `correct` | The probe verdict matches the label, and both states match the reference when it ran |
| `incorrect` | An equivalent pair gave different states, a state's norm drifted by more than 1e-6, or a state differs from the reference. This is a simulator, DAG or normalization bug |
| `undetected` | A non-equivalent pair whose probe states all matched. Expected only with `probes = 1` (the generator's appended X can leave `\|0...0>` unchanged up to phase) |
| `timeout` | The worker was killed at `timeout` |
| `error` | An exception, allocation failure or crash; the message is in `error`. Retried by `--resume` |
| `skipped` | Not run, because auto-skip had already stopped this series |
| `skipped_memory` | Not run, because the state would exceed `memory_budget_mb` |

### How correctness is checked

A simulator only computes `U|psi>` for one input, so the pair check uses probes.
For each probe the same random layer of single-qubit U gates is prepended to A
and B, which preserves equivalence. The verdict is *equivalent* when every
probe gives fidelity ≥ 1 − `tol` **and** the same overlap phase, which matches
the generator's definition (equal up to one global phase). Probe 0 is the plain
`|0...0>` input and is the one that is timed.

The label check alone only tests that the simulator is self-consistent. A bug
that relabels every two-qubit gate the same way (for example swapping the
`idx01`/`idx10` roles in `apply_2q_gate`) keeps equivalent pairs equivalent and
passes it. The reference simulation catches that: it applies each gate by a
plain loop over all basis indices, using the IR's documented convention (local
index `(bit(q[1]) << 1) | bit(q[0])`, row-major matrices, qubit 0 least
significant), without sim.cpp's index-insertion arithmetic. Keep
`reference_max_qubits` above zero in at least one sweep you run regularly.

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
- `normalized` pays `prep_s` once per circuit to save simulated gates; for
  repeated execution of one circuit, compare its `sim_s` with `flat`.
- The generator's pairs are near-identical circuits with shared prefixes. They
  are good for checking correctness on every gate kind and for scaling curves,
  not a model of compiled application circuits.
- `results.jsonl` is plain JSON lines: `pd.read_json(path, lines=True)`.

## Adapting to the repository

`sim_bench_adapter.hpp` is the only file that calls QCanopy APIs beyond
`sim.hpp`: `testgen::read_ir(std::istream&)`, `testgen::normalized_copy(const
Circuit&)` and `CircuitDAG(const Circuit&)`. If a signature differs, change it
there. On the Python side, `read_manifest()` in `sim_profile.py` finds each
case's label and A/B paths by key name in `manifest.jsonl`; if the generator's
field names change and it cannot find them, it stops and prints the record.
