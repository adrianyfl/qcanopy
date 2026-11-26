# QCanopy: Batched, Reuse-Aware GPU Simulation for Large-Scale Quantum Circuit Equivalence Checking

QCanopy provides a quantum circuit IR and dependency DAG, normalization and
prefix/suffix peeling passes, a deterministic circuit-pair generator, and a
parser for the generator's OpenQASM 3.0 files.
`qcanopy_generate` creates labeled equivalent and non-equivalent pairs for
equivalence-checker tests and experiments on shared computation.
`qcanopy_parse` and the `qcanopy::parser` library read those files back
into the IR.

**New users:** follow the [test-case generator quick start](docs/generator_quickstart.md)
to build the tool, generate a dataset, and understand its files, then the
[parser quick start](docs/parser_quickstart.md) to read the dataset back.

## Build

Requirements: CMake 3.18+ and a C++17 compiler with working `<filesystem>` and
integer `from_chars` support. The generator runs on the CPU and does not require CUDA.

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Check two circuits

The `qcanopy` command-line tool (`build/src/cli/qcanopy`) checks whether two
circuits are equal up to a global phase:

```sh
./build/src/cli/qcanopy check a.qasm b.qasm      # or just: qcanopy a.qasm b.qasm
```

`check` runs the whole pipeline: it parses both files, lets the
[Python rewriter](rewriting/README.md) remove what it can prove equal, and simulates
whatever is left on random (or, with `--stimulus all-basis`, all basis) input states.
Exit status: 0 equivalent, 1 not equivalent, 2 error. Run `qcanopy check --help` for
the options. If Python 3.9+ is missing, `check` warns and simulates the full circuits.

Each stage can also be used on its own:

| Command | What it does |
|---|---|
| `qcanopy parse FILE... [--emit qasm\|ir]` | Summarize circuits, or convert between OpenQASM 3.0 and JSON IR |
| `qcanopy rewrite normalize FILE` | C++ normalization: canonical gates, fused rotations, cancelled inverses |
| `qcanopy rewrite peel A B` | C++ normalization + common prefix/suffix peeling; prints both residuals |
| `qcanopy rewrite prove A B` | Python pair prover: proved (exit 0) or inconclusive (exit 1) |
| `qcanopy rewrite classify FILE` | Python boundary classifier |
| `qcanopy sim FILE [--input SPEC]` | Simulate a circuit and print its most likely output amplitudes |
| `qcanopy generate ...` | Generate labeled test-case pairs (same options as `qcanopy_generate`) |

Circuit files are OpenQASM 3.0, or JSON IR v1 if their names end in `.json`.
`qcanopy_parse` and `qcanopy_generate` remain available and behave like
`qcanopy parse` and `qcanopy generate`.

## Generate your first dataset

On Linux/macOS with a single-configuration build:

```sh
./build/src/generator/qcanopy_generate --format ir --output generated/quickstart-ir --seed 759
./build/src/generator/qcanopy_generate --format qasm --output generated/quickstart-qasm --seed 759
```

For Windows executable paths and PowerShell examples, see the
[quick start](docs/generator_quickstart.md#2-generate-test-cases).

Each command produces **32 pairs: 16 equivalent and 16 non-equivalent**, using
4 qubits and 100 base gates by default. Each dataset contains 64 circuit files
and a `manifest.jsonl` with labels, seeds, construction proofs, verification
status, and raw/normalized sharing measurements. Both formats represent the
same cases when generation arguments match. Each output directory must be absent
or empty.

## Parse a generated dataset

```sh
./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/*.qasm
```

This prints the qubit count, gate count, and global phase of each circuit, then
`parsed 64 of 64 files`. From C++, include `qcanopy/parser.hpp` and link
`qcanopy::parser`:

```cpp
qcanopy::Circuit    c   = qcanopy::parse_qasm_file("generated/quickstart-qasm/cases/case-000000.a.qasm");
qcanopy::CircuitDAG dag = qcanopy::parse_qasm_file_dag("generated/quickstart-qasm/cases/case-000000.a.qasm");
```

The parser accepts OpenQASM 3.0 as the generator writes it, plus comments and
arithmetic such as `pi/2`. It does not read OpenQASM 2.0 files, gates outside the
generator's set, or circuits with more than 32 qubits. Errors are reported as
`file:line:column: message`.

## Generator capabilities

- Select portable JSON IR v1 or OpenQASM 3.0 output with `--format ir|qasm`.
- Control qubit count, base gate count, gate pool, angles, label balance, and seed.
- Plant shared prefixes within pairs and across groups of pairs.
- Generate circuits with **1-32 qubits**. Default dense verification checks the
  complete operator through **5 qubits**; larger cases retain construction proofs
  and record that numerical verification was skipped.
- Reuse generation and export functions through the `qcanopy::generator` C++ library.

Gate count is not dependency depth. The prefix option controls planted sharing;
actual sharing can change after transformations and normalization. A full shared
prefix does not establish circuit equivalence.

## Documentation

- [Quick user guide](docs/generator_quickstart.md): build, commands, key options,
  output files, verification, and common errors.
- [Generator reference](docs/test_case_generator.md): complete CLI contract,
  constructions, reproducibility, serialization schemas, and C++ API.
- [Parser quick start](docs/parser_quickstart.md): parsing generated datasets,
  the C++ calls, accepted input, how parsed circuits relate to IR files, and errors.
- [Peeling algorithm](src/rewrite/algorithm.md): matching prefix/suffix removal.
- [Python DAG rewriting interface](rewriting/README.md): SWAP removal with output
  mappings, fixed-library boundary classification, optional overlap, and bounded
  circuit-pair proofs. This standalone Python module exchanges JSON with the
  generator; it does not require a C++ build or provide simulator integration.

Datasets written under `generated/` are Git-ignored. Generate them locally rather
than expecting sample datasets or workspace-specific build tools in a fresh clone.
