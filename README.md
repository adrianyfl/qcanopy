# QCanopy: Batched, Reuse-Aware GPU Simulation for Large-Scale Quantum Circuit Equivalence Checking

QCanopy provides a quantum circuit IR and dependency DAG, normalization and
prefix/suffix peeling passes, a deterministic circuit-pair generator, and a
parser for OpenQASM 3.0 and 2.0 circuit files.
`qcanopy_generate` creates labeled equivalent and non-equivalent pairs for
equivalence-checker tests and experiments on shared computation.
`qcanopy_parse` and the `qcanopy::parser` library read the generator's
OpenQASM 3.0 files, and OpenQASM 2.0 benchmark files, into the IR.

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

## Parse circuit files

```sh
./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/*.qasm
```

This prints the qubit count, gate count, and global phase of each circuit, then
`parsed 64 of 64 files`. A file that fails is reported as
`file:line:column: message` and the rest are still parsed. From C++, include
`qcanopy/parser.hpp` and link `qcanopy::parser`:

```cpp
qcanopy::Circuit    c   = qcanopy::parse_qasm_file("generated/quickstart-qasm/cases/case-000000.a.qasm");
qcanopy::CircuitDAG dag = qcanopy::parse_qasm_file_dag("generated/quickstart-qasm/cases/case-000000.a.qasm");
```

The version header on the first line of a file selects the dialect:

| | OpenQASM 3.0 | OpenQASM 2.0 |
| --- | --- | --- |
| Source | `qcanopy_generate --format qasm` | Benchmark circuits, such as the Quokka# set |
| Library | `stdgates.inc` | `qelib1.inc` |
| Registers | One `qubit[n]` | Any number of `qreg`; `creg` is ignored |
| Gates | `id h x y z s sdg t tdg rx ry rz p U cx cz swap` | `U`, `CX`, all of `qelib1.inc` (including `ccx`), and gates the file defines with `gate` |

Both accept comments and arithmetic such as `pi/2`. Gates on three or more
qubits are expanded into one- and two-qubit gates from their definitions; `ccx`
becomes 15 gates.

Not supported: circuits with more than 32 qubits, `measure`, `reset`, `barrier`,
`if`, and gates that are in neither the library nor the file. Of the 7,327
Quokka# benchmark files, 6,365 parse; 947 of the rest exceed 32 qubits.

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
- [Parser quick start](docs/parser_quickstart.md): parsing OpenQASM 3.0 and 2.0
  files, the C++ calls, accepted input, how parsed circuits relate to their
  source, and errors.
- [Peeling algorithm](src/rewrite/algorithm.md): matching prefix/suffix removal.

Datasets written under `generated/` are Git-ignored. Generate them locally rather
than expecting sample datasets or workspace-specific build tools in a fresh clone.
