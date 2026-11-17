# QCEC profiling results

Two runs of the MQT QCEC 3.12.0 **simulation checker** on random circuit pairs from
`qcanopy_generate` (mixed gate set, mixed angles, 50% shared prefix, half equivalent and
half non-equivalent). Each run directory has:

| File | Contents |
| --- | --- |
| `run.json` | sweep settings, QCEC/mqt.core versions, host, git commit |
| `results.jsonl` | one JSON record per check (see the profiler README for the fields) |
| `report.html` | the profiler's report |
| `progress_chart.html` | heatmap of median check time per (qubits, gates) point, with a table view |

The generated QASM circuits are not committed. Running the same sweep file regenerates them
exactly, because the generator is deterministic for a given seed and commit (`bfe229b`).

## Setup

- AMD Ryzen 7 9700X (8 cores, 16 threads), Docker on WSL2, container capped at 14 GB of memory.
- One check at a time (`--jobs 1`). The simulation checker runs 16 random-input simulations
  in parallel threads, so each check uses all 16 threads.
- `preprocessing_s` is QCEC's circuit optimization time; `check_s` is the simulation itself.

| Run | Sweep | Qubits | Pairs per point | Limit |
| --- | --- | --- | --- | --- |
| `20261114-032142-gate-scaling` | `sweeps/gate_scaling.toml` | 4–20 | 8 | 60 s timeout; a series stops after 2 timeouts |
| `20261114-033901-gate-scaling-nolimit` | `sweeps/gate_scaling_nolimit.toml` | 4–32 | 4 | none; checks run until done or out of memory |

Base gate counts are 100, 300, 1k, 3k and 10k per circuit (a pair has about twice that).

## Results (no-timeout run)

Median simulation time, equivalent / non-equivalent pairs:

| Qubits | 100 gates | 300 | 1k | 3k | 10k |
| --- | --- | --- | --- | --- | --- |
| 4 | 5 / 2 ms | 5 / 4 ms | 8 / 7 ms | 22 / 18 ms | 63 / 51 ms |
| 8 | 5 / 3 ms | 12 / 12 ms | 59 / 51 ms | 173 / 177 ms | 0.63 / 0.50 s |
| 12 | 5 / 4 ms | 0.12 / 0.21 s | 1.0 / 0.87 s | 4.7 / 3.6 s | 13 / 11 s |
| 16 | 5 / 11 ms | 0.92 / 1.2 s | 42 / 39 s | 2.9 / 2.4 min | 8.9 / 8.5 min |
| 20 | 5 / 4 ms | 1.3 / 12 s | 21 / 19 min | 1.6 / 1.5 h | 5.6 / 5.4 h |
| 24 | 4 / 3 ms | 4.8 / 13 s | out of memory | out of memory | out of memory |
| 28 | 6 / 7 ms | 1.4 / 5.4 min | out of memory | skipped | skipped |
| 32 | 6 / 5 ms | 3.2 / 5.8 min | out of memory (1), skipped (3) | skipped | skipped |

"Skipped" means the profiler did not run the check because a lower qubit count had already
run out of memory at that gate count.

## Findings

- **Optimization is negligible.** QCEC's preprocessing takes 3–19 ms at every size; simulation
  ranges from milliseconds to hours.
- **Simulation cost is set by entanglement, not gate count alone.** Each qubit count grows much
  faster than linearly in gates until the decision diagram reaches its maximum size (about 2^n
  nodes), then roughly linearly. At 20 qubits that happens by 1k gates: from 1k to 10k gates
  the time grows about linearly and peak memory stays at 8.3–8.9 GB.
- **24 qubits and up run out of 14 GB** from 1k gates on. At 300 gates they still finish,
  because the circuits are too short to fill the diagram.
- **Three wrong verdicts at 100 gates** (20q in the first run; 20q and 24q in this one):
  non-equivalent pairs reported as `probably_equivalent`. The simulation checker uses
  basis-state inputs, which can miss a difference such as an X applied to a qubit already in
  |±⟩. The first-run case was confirmed by hand; the other two have not been checked.
- These are worst-case inputs for this checker. QCEC's default configuration (with the
  alternating checker) has not been profiled on these pairs yet.

## Regenerating the charts

```sh
python profiling/qcec/progress_chart.py profiling/qcec/results/<run dir>
```
