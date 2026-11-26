# QCanopy: Batched, Reuse-Aware GPU Simulation for Large-Scale Quantum Circuit Equivalence Checking

QCanopy checks whether two quantum circuits are equal up to a global phase. The
`qcanopy` command-line tool runs the whole pipeline in one step: it parses both
circuits, uses rewriting to remove the parts it can prove equal, and simulates only
what is left. Each stage (parser, rewriter, simulator, test-case generator) can also
be run on its own, or used as a C++ library.

```sh
$ qcanopy check a.qasm b.qasm
equivalent up to global phase e^{i*0} (proved by the rewriter in 3 steps)
  gates: 100 vs 101, left after rewriting 0 vs 0
```

## Build

Requirements: CMake 3.18+ and a C++17 compiler with working `<filesystem>` and
integer `from_chars` support. Python 3.9+ (standard library only) is optional: `check`
uses it to run the [rewriter](rewriting/README.md), and falls back to plain simulation
without it. Everything currently runs on the CPU; CUDA is not required.

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The tool is `build/src/cli/qcanopy`. The examples below assume it is on your `PATH`,
for example after `export PATH="$PWD/build/src/cli:$PATH"`.

## Quick start

Generate a labeled dataset of 32 circuit pairs (16 equivalent, 16 not), then check
some of them:

```sh
qcanopy generate --format qasm --output generated/quickstart-qasm --seed 759
cd generated/quickstart-qasm/cases
qcanopy check case-000000.a.qasm case-000000.b.qasm
qcanopy check case-000003.a.qasm case-000003.b.qasm
```

```text
NOT equivalent: stimulus #0 (haar) on the rewriter's residual: output is not a phase multiple of the input (deviation 1.3)
  gates: 100 vs 101, left after rewriting 0 vs 1
equivalent up to global phase e^{i*0} (proved by the rewriter in 3 steps)
  gates: 100 vs 101, left after rewriting 0 vs 0
```

Each pair's expected answer is in `generated/quickstart-qasm/manifest.jsonl`
(`expected_relation`). `qcanopy A B` is short for `qcanopy check A B`.

## Commands

| Command | What it does |
|---|---|
| `qcanopy check A B` | **Main command.** Equivalence check up to global phase: parse, rewrite, simulate |
| `qcanopy parse FILE...` | Summarize circuits, or convert between OpenQASM 3.0 and JSON IR |
| `qcanopy rewrite normalize FILE` | C++ normalization of one circuit |
| `qcanopy rewrite peel A B` | C++ normalization, then common prefix/suffix removal |
| `qcanopy rewrite prove A B` | Python pair prover on its own |
| `qcanopy rewrite classify FILE` | Python boundary classifier on its own |
| `qcanopy sim FILE` | Simulate a circuit and print its output amplitudes |
| `qcanopy generate ...` | Generate labeled equivalent / non-equivalent test-case pairs |

Every command takes `--help`. Circuit files are OpenQASM 3.0, or JSON IR v1 if their
names end in `.json`. Exit status is 0 on success and 2 on a usage, parse or I/O error;
commands with a yes/no answer use 1 for "no" (see below).

### `qcanopy check`

```text
qcanopy check A B [options]
```

`check` runs four steps:

1. **Parse** both files. They must have the same number of qubits.
2. **Rewrite.** The Python pair prover (`rewriting/run.py pair`) removes common
   boundaries, gate/inverse pairs and verified library identities from both circuits.
   If nothing is left, the pair is proved equivalent.
3. **Simulate** whatever the rewriter left on a batch of input states (by default 8
   Haar-random and 8 random product states). The leftovers are equivalent exactly when
   the original circuits are, so only they need simulating.
4. **Measure the phase.** The rewriter does not track global phase, so after an
   "equivalent" answer the full circuits are run on one random state to measure it.
   This run also catches a wrong rewriter proof.

| Option | Meaning |
|---|---|
| `--stimulus LIST` | Input states, comma-separated (default `haar,product`): `haar`, `product`, `basis` (random basis states), `all-basis` (every basis state; at most 12 qubits) |
| `--count N` | States per random kind (default 8) |
| `--seed UINT64` | Seed for the random states (default 759) |
| `--tol X` | Tolerance on the residual norm and phase (default 1e-9) |
| `--no-rewrite` | Skip the rewriter and simulate the full circuits |
| `--overlap` | Let the rewriter try overlapping rule substitutions (slower, proves more) |
| `--json` | Print the result as one JSON object |
| `--quiet` | Print nothing; report the answer by exit status only |

**Exit status:** 0 equivalent, 1 not equivalent, 2 error.

**How sure is the answer?**

- **Not equivalent** is always definitive. The output names the input state that
  told the circuits apart.
- **Equivalent** is a proof when the rewriter proved it, or when `--stimulus all-basis`
  passed. With random stimuli it means no state told the circuits apart; for
  non-equivalent circuits that happens with probability zero (up to `--tol`). The
  output says "not a proof" in this case.

**Without the rewriter.** If Python is missing or the rewriter fails, `check` prints a
warning on stderr and simulates the full circuits, so the answer is still correct.
`QCANOPY_PYTHON` sets the interpreter (default: the Python 3 found when CMake ran)
and `QCANOPY_REWRITER` sets the directory holding `run.py` (default: `rewriting/` in
this repository). The rewriter only runs on Linux and macOS.

### `qcanopy parse`

```sh
qcanopy parse cases/*.qasm                     # one summary line per file
qcanopy parse case-000000.a.qasm --emit ir     # OpenQASM -> JSON IR
qcanopy parse circuit.json --emit qasm         # JSON IR -> OpenQASM
```

The summary prints each file's qubit count, gate count and global phase, then
`parsed N of M files`. Exit status: 0 if every file parsed, 1 if any failed, 2 on a
usage error. Parse errors read `file:line:column: message`.

The parser accepts OpenQASM 3.0 as the generator writes it (`OPENQASM 3.0;`,
`include "stdgates.inc";`, one `qubit[n]` register, the gates `id h x y z s sdg t tdg
rx ry rz p U cx cz swap`, and `gphase`), plus comments and arithmetic such as `pi/2`.
It does not read OpenQASM 2.0, other gates, or circuits with more than 32 qubits.

### `qcanopy rewrite`

| Subcommand | Output |
|---|---|
| `normalize FILE [--emit qasm\|ir]` | The circuit with canonical gate spellings, fused rotations and neighbouring inverses cancelled. Exact, including global phase. |
| `peel A B [--emit qasm\|ir] [--no-normalize]` | What is left of each circuit after normalizing and removing their common prefix and suffix. |
| `prove A B [--overlap] [--json]` | The Python prover's verdict and proof steps. Exit 0 proved, 1 inconclusive. |
| `classify FILE [--overlap] [--json]` | A summary of the Python classifier's matches and regions. Exit 0 complete, 1 stopped early. |

`normalize` and `peel` are C++ passes ([peeling algorithm](src/rewrite/algorithm.md)).
`prove` and `classify` run the [Python rewriting module](rewriting/README.md);
`--json` prints its full output. An inconclusive `prove` does not mean the circuits
differ: use `check` for an answer.

### `qcanopy sim`

```sh
qcanopy sim circuit.qasm                       # from |0...0>, 16 most likely states
qcanopy sim circuit.qasm --input basis:5 --all
qcanopy sim circuit.qasm --input haar:42 --json
```

`--input` is `zero` (default), `basis:INDEX`, `haar:SEED` or `product:SEED`. Basis
states are written `|q[n-1] ... q[0]>`. `--top K` shows the K most likely states
(default 16) and `--all` shows every state with nonzero probability. The state vector
takes 2^n × 16 bytes of memory (16 GiB at 30 qubits).

### `qcanopy generate`

```sh
qcanopy generate --format qasm --output generated/run --qubits 4 --gates 100 --pairs 32 --seed 759
```

Writes labeled circuit pairs as OpenQASM 3.0 or JSON IR v1, plus a `manifest.jsonl`
with labels, seeds, construction proofs, verification status and sharing measurements.
The output directory must be absent or empty.

- Control qubit count (1–32), base gate count, gate set, angle mode, label balance and seed.
- Plant shared prefixes within pairs and across groups of pairs.
- Dense verification checks every label through 5 qubits; larger cases keep their
  construction proofs and record that numerical verification was skipped.

Gate count is not dependency depth, and a shared prefix does not make two circuits
equivalent. See the [generator quick start](docs/generator_quickstart.md) and the
[generator reference](docs/test_case_generator.md) for every option.

`qcanopy_generate` (in `build/src/generator/`) and `qcanopy_parse` (in
`build/src/parser/`) are still built for existing scripts. They behave exactly like
`qcanopy generate` and `qcanopy parse`.

## Using the C++ libraries

Each stage is a CMake target you can link:

| Target | Header | Provides |
|---|---|---|
| `qcanopy::ir` | `qcanopy/ir.hpp`, `qcanopy/dag.hpp` | `Circuit`, gates, and the dependency DAG `CircuitDAG` |
| `qcanopy::parser` | `qcanopy/parser.hpp` | `parse_qasm_file()`, `parse_qasm_file_dag()` |
| `qcanopy::rewrite` | `qcanopy/normalize.hpp`, `qcanopy/peel.hpp` | `normalize()`, `peel()` |
| `qcanopy::sim` | `qcanopy/sim.hpp` | `Simulator`, `StateVector`, `Stimulus` |
| `qcanopy::equiv` | `qcanopy/equiv.hpp`, `qcanopy/stimulus.hpp` | `EquivalenceChecker`, stimulus generators |
| `qcanopy::generator` | `qcanopy/generator.hpp`, `qcanopy/dataset_io.hpp` | Pair generation and dataset I/O (`qcanopy::testgen::read_ir()`, `write_ir()`, `write_qasm()`) |

```cpp
#include "qcanopy/equiv.hpp"
#include "qcanopy/parser.hpp"
#include "qcanopy/stimulus.hpp"

const qcanopy::Circuit a = qcanopy::parse_qasm_file("a.qasm");
const qcanopy::Circuit b = qcanopy::parse_qasm_file("b.qasm");
const qcanopy::EquivResult r =
    qcanopy::EquivalenceChecker{}.check(a, b, qcanopy::haar_random_states(a.n, 16, /*seed=*/759));
// r.equivalent, r.global_phase, and the failing stimulus if not equivalent
```

`EquivalenceChecker` simulates only; the rewriter is used by the `qcanopy` tool.

## Repository layout

| Path | Contents |
|---|---|
| `include/qcanopy/` | Public C++ headers |
| `src/ir`, `src/parser`, `src/rewrite`, `src/sim`, `src/equiv`, `src/generator` | The C++ libraries above |
| `src/cli/` | The `qcanopy` tool: one file per command (`commands.hpp` lists them) |
| `rewriting/` | Python rewriting module: verified rule library, classifier and pair prover |
| `tests/` | C++ unit tests and CLI tests (run with `ctest`) |
| `docs/` | Generator and parser guides |

## Documentation

- [Generator quick start](docs/generator_quickstart.md): build, commands, key options,
  output files, verification and common errors.
- [Generator reference](docs/test_case_generator.md): the full option list,
  constructions, reproducibility, file formats and C++ API.
- [Parser quick start](docs/parser_quickstart.md): parsing datasets, the C++ calls,
  accepted input, and how parsed circuits relate to IR files.
- [Peeling algorithm](src/rewrite/algorithm.md): prefix/suffix removal.
- [Python rewriting module](rewriting/README.md): SWAP removal, boundary
  classification, overlap search and pair proofs, plus how its rule library was verified.

## Current limits

- The simulator is a CPU state-vector simulator. GPU execution and reuse of shared
  computation across circuits are planned but not implemented yet.
- The rewriter runs as a separate Python process (about 0.2–0.5 s per pair). For very
  small circuits `check --no-rewrite` can be faster.
- The rewriter runs on Linux and macOS only; on Windows `check` always simulates the
  full circuits.

Datasets written under `generated/` are Git-ignored. Generate them locally rather
than expecting sample datasets in a fresh clone.
