# Circuit-pair test-case generator

For your first dataset, start with the [quick user guide](generator_quickstart.md).
This page is the detailed format, construction, and API reference.

`qcanopy_generate` produces deterministic labeled circuit pairs using the current
`qcanopy::Circuit` IR. The reusable `qcanopy::generator` library separates generation,
verification, and export. It depends on the existing IR and DAG normalization
libraries, with no third-party C++ or runtime dependencies.

Integration baseline: upstream `main` at
`049e25b5c72dd769e62cfa6a0b506b5241043f0a` ("Added dag to circuit, and basic dag
normalization"). The build filename `src/ir/CMakeLIsts.txt` was corrected to
`CMakeLists.txt` for case-sensitive systems. The unfinished `sim` branch is unused.

## Build and run

Use CMake 3.18+ and a compiler with a working C++17 standard library, including
`<filesystem>` and integer `from_chars` (for example current GCC, Clang, or MSVC).

```sh
cmake -S . -B build -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The executable is `build/src/generator/qcanopy_generate` on single-configuration
Unix builds, `build/src/generator/qcanopy_generate.exe` on Ninja/Windows, and
`build/src/generator/Release/qcanopy_generate.exe` with Visual Studio Release builds.
Use the appropriate path below:

```sh
build/src/generator/qcanopy_generate --format qasm --qubits 4 --gates 100 --pairs 32 --prefix-fraction 0.75 --equivalent-fraction 0.5 --prefix-group-size 4 --seed 759 --output generated/qasm-run
build/src/generator/qcanopy_generate --format ir --qubits 4 --gates 100 --pairs 32 --prefix-fraction 0.75 --equivalent-fraction 0.5 --prefix-group-size 4 --seed 759 --output generated/ir-run
```

These runs each produce 32 pairs (16 equivalent and 16 non-equivalent), 64 circuit
files, and a `manifest.jsonl` with 33 records. Identical generation arguments produce
the same circuits and case metadata in both formats. The selected format only
changes serialization and file paths. Generated data and workspace-local build tools
are ignored by Git.

The destination must be absent or an empty directory. Existing data is refused.
`manifest.jsonl` is published last; incomplete writes retain a `.partial` manifest
and must not be treated as a completed dataset. A failure can leave partial circuit
files; choose another destination or inspect those files before removing them.
Do not have multiple writers target the same directory concurrently.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--format qasm\|ir` | required | OpenQASM 3.0 or portable circuit JSON v1 |
| `--output DIR` | required | Absent or empty destination |
| `--qubits N` | 4 | 1 through `kMaxQubits` (currently 32) |
| `--gates L` | 100 | Base gate count; zero allowed; **not dependency depth** |
| `--pairs B` | 32 | Positive batch size |
| `--prefix-fraction R` | 0.5 | Plant `floor(R*L)` initial gates, `0 <= R <= 1` |
| `--equivalent-fraction E` | 0.5 | Exactly `floor(E*B+0.5)` equivalent labels |
| `--seed UINT64` | 759 | Decimal unsigned 64-bit master seed |
| `--gate-set SET` | mixed | `clifford`, `clifford-t`, `parameterized`, or `mixed` |
| `--angle-mode MODE` | mixed | `discrete`, `continuous`, or `mixed` |
| `--prefix-group-size G` | 1 | Positive size of consecutive groups sharing a prefix |
| `--verify small\|none` | small | Dense verification up to five qubits, or skip |
| `--help` | | Show usage |

Integer parsing accepts decimal digits without signs or trailing text and rejects
overflow. Counts must also fit the platform's `size_t`, vector capacity, and DAG's
32-bit node indices, with room for inserted gates. Fractions must be finite and in
range. Unknown, missing, or duplicate options fail. Exit status is 0 for success
(including help) and nonzero for invalid arguments, generation/verification errors,
or I/O failures. The entire batch is retained in memory; select batch size and base
gate count appropriate to available RAM. No exponential state is allocated for
large-qubit generation or sharing measurement.

## Gate mix and transformations

Sampling is uniform over eligible gate kinds, then uniform over operand qubits.
Two-qubit operands are distinct; one-qubit circuits exclude two-qubit gates.

* `clifford`: H, X, Y, Z, S, Sdg, SX, SXdg, and CX/CZ/SWAP.
* `clifford-t`: Clifford pool plus T and Tdg.
* `parameterized`: RX, RY, RZ, P, U, and CX/CZ/SWAP.
* `mixed`: union of these pools.

The gate-set restriction applies to base and prefix generation. Equivalent rewrites
and the appended negative X can introduce other named gates (for example S -> T T
in a Clifford base). Every generated operation is exportable to both formats.
Required parameter slots are populated through `make_gate`; unused slots are zero.
Discrete angles are the 17 multiples `k*pi/4`, `k in [-8,8]`. Continuous angles use
`(2*u-1)*2*pi`, with values within `1e-6` of a multiple of `pi/4` shifted by `1e-5`
to avoid normalization tolerance boundaries. Mixed mode chooses each mode with
probability 1/2 independently for each required angle.

Transformation requests cycle with `case_index % 7`: identical, inverse-pair
insertion, X -> H Z H, S -> T T, SWAP -> three CX gates, adjacent disjoint swap,
and a circuit-wide phase change of `pi/7`. Replacements use the first eligible
suffix location. If a replacement or commutation has no eligible suffix location,
it falls back to inserting G followed by its inverse. Insertion uses a sampled
existing gate (H for an empty circuit) at a seeded position at or after the planted
boundary. Metadata records the **actual** operation, zero-based location, removed
and inserted gates (including parameters), phase adjustment, and proof. An inverse
circuit reverses the sequence, inverts each gate, and negates the global phase.

Labels are shuffled deterministically after assigning the rounded equivalent count.
Every negative case first receives one of the equivalent constructions, then gets
an X on a seeded qubit appended to B. If `V = exp(i*phi) A`, then
`B = X(q) V` cannot equal a scalar multiple of A: multiplication by the inverse of
A would require X to be scalar. This proof holds for every supported qubit count.
No unchecked random gate perturbation receives a negative label.

Some equivalent pairs normalize to the same gate sequence (inverse cancellation,
rotation fusion); others can retain different structure (H Z H versus X, or three
CX versus SWAP). Structural difference after normalization does not imply semantic
inequivalence.

## Prefixes and normalized sharing

Consecutive groups of G pairs explicitly reuse the same generated prefix vector;
the final group may contain fewer than G cases. `prefix_identity` is a dataset-local
group identity, reproducible with its prefix seed and config; it is not a global
content hash. Suffix and transformation streams remain independent per case.

`--prefix-fraction` is a **planted-prefix control**, not a target final percentage.
Transformations preserve the planted gates. An accidental suffix match can extend
the common prefix, and length-changing transformations change the denominator.
Raw sharing uses `same_gate` with absolute matrix tolerance `1e-12`, including its
upstream symmetric-operand semantics for CZ/SWAP. Gate parameters are validated
against factory matrices. Global phase is tracked separately and does not shorten
the gate prefix.

For each raw/normalized pair the manifest records both gate counts, the longest
shared initial gate count, and that count divided by the shorter gate count. Two
empty gate lists have fraction 1; exactly one empty list has fraction 0.
At R=1, a negative pair can have fraction 1 because the appended X follows the
entire shorter circuit. Fraction 1 therefore does not establish equivalence.

Normalized statistics use copies converted through `CircuitDAG`, `normalize`, and
`to_circuit`. Canonical DAG traversal can reorder disjoint gates, and normalization
can cancel or fuse across a prefix boundary. Planted raw sharing does not guarantee
normalized sharing. Output files retain the original raw circuits. Dependency depth
and targeted normalized-prefix bins are not computed.

## IR JSON v1

Each file is one JSON object with exactly these members:

```json
{"schema_version":1,"n":2,"global_phase":0.25,"gates":[{"name":"h","qubits":[0],"params":[]},{"name":"cx","qubits":[0,1],"params":[]},{"name":"rz","qubits":[1],"params":[0.5]}]}
```

The sequence is in application order. Names are `id`, `h`, `x`, `y`, `z`, `s`, `sdg`,
`t`, `tdg`, `sx`, `sxdg`, `rx`, `ry`, `rz`, `p`, `u`, `cx`, `cz`, `swap`,
`unitary_1q`, and `unitary_2q`. Named gates carry exactly their required number of
parameters (0, 1, or 3) and operands (1 or 2); matrices are rebuilt by `make_gate`.
Generic gates carry `params:[]` and an additional `matrix` array of exactly 4 or
16 row-major complex entries, each `[real,imaginary]`. Generic matrices must be
finite and unitary within `1e-9`. Generic gates are supported by the library's
IR reader/writer, but are not sampled by the default generator.

Qubit 0 is the least significant state bit. The local two-qubit index is
`(bit(q[1]) << 1) | bit(q[0])`; CX's first operand is control. `global_phase` and
angles are radians. Doubles use 17 significant decimal digits with the classic
locale; finite doubles round-trip. Enum ordinals, memory layout, and pointers are
not part of the schema. The strict reader rejects unknown versions/names, missing,
duplicate or extra members, malformed JSON, non-finite/out-of-range numbers,
wrong parameter counts, inconsistent dimensions, and nonunitary generic gates.
Nesting is capped at 64 levels. This is a portable new format, not an upstream
QCanopy serialization convention.

## OpenQASM 3.0 export

Files declare `OPENQASM 3.0;`, include `stdgates.inc`, and allocate `qubit[n] q`.
Gate order and operand order are preserved. The exporter emits one accumulated
`gphase` instruction and the original operation sequence, lowering SX/SXdg to
`rx(±pi/2)` while adding phase `±pi/4`.

The IR's phase-free U is emitted as OpenQASM's uppercase built-in `U`, with phase
`-theta/2` added to the circuit's exported `gphase`. OpenQASM 3 defines its built-in
U with an additional factor `exp(i*theta/2)`. This correction is necessary to
preserve the **exact** IR operator. Named standard gates use their standard names,
including `id`. See the [OpenQASM 3.0 gate definitions](https://openqasm.com/versions/3.0/language/gates.html)
and [tagged standard library](https://github.com/openqasm/openqasm/blob/3.0.0/examples/stdgates.inc).

Generic U1Q/U2Q export fails clearly before publishing a dataset; use IR for those
matrices. There is no production QASM parser, QASM 2 export, or generic decomposition.
The optional validation tool below parses only for tests and interprets the
generator's straight-line subset; it is not a general import API.

## Manifest and reproducibility

`manifest.jsonl` starts with a header containing schema/generator versions, build
base source commit (or `unknown`), output format/dialect, complete config, PRNG name,
verification cap, and prefix tolerance. The source commit names the checkout HEAD
at CMake configuration time, not a commit containing uncommitted generator edits.
Each following case record includes expected relation, IDs, all component seeds,
planted/actual sharing, raw/normalized sizes, phases, transform evidence, negative
mutation proof, verification status, and relative circuit paths. Seeds are decimal
JSON strings so consumers preserve all 64 bits. Other counts are JSON numbers.
Fields have schema version 1; case IDs start at `case-000000`.

The PRNG is SplitMix64. In unsigned 64-bit arithmetic, mixing is:

```text
mix(z): z=(z xor (z>>30))*0xbf58476d1ce4e5b9
        z=(z xor (z>>27))*0x94d049bb133111eb
        return z xor (z>>31)
derive(master,domain,index) = mix(master xor mix(domain) xor mix(index+0x9e3779b97f4a7c15))
next(): state += 0x9e3779b97f4a7c15; return mix(state)
```

Domains 1, 2, and 3 derive prefix/group, suffix/case, and transform/case seeds;
domain 4/index 0 derives the label-shuffle seed. The negative-qubit seed is
`derive(transform_seed,5,0)`. Bounded sampling rejects `next()` values below
`(-bound mod 2^64) % bound`, then takes `% bound`. Floating samples use the upper
53 bits times `2^-53`. Labels use descending Fisher-Yates. `generator.cpp` specifies
the exact draw order; this version does not use `std::uniform_*` or `std::hash`.
Random choices and angle sampling are specified across platforms. Transcendental
gate-matrix calculations and decimal formatting can differ in the last bit across
standard libraries; byte-for-byte output reproducibility is tested on a fixed
toolchain. Generation semantics do not depend on output format.

With `--verify small`, all cases at N<=5 compare complete dense operators with
one global phase for the entire operator, at absolute tolerance `1e-8`. Both
normalized copies are also checked against the original operators **including**
phase. The reference applies little-endian local gates to every basis column.
Status is `passed`, `skipped_by_request`, or `skipped_size_limit`; skipped cases
have null numerical method/tolerance and retain construction proofs. Dense memory
is O(4^N); it is never allocated above the cap. Structural validation still runs
on all sizes, including finite parameters, named factory consistency, and generic
unitarity. Verification tolerance is distinct from upstream normalization's
`kAngleTol=1e-10` and the prefix tolerance.

## Library use and independent export checks

Link `qcanopy::generator`; include `qcanopy/generator.hpp` and
`qcanopy/dataset_io.hpp`. Example:

```cpp
qcanopy::testgen::GeneratorConfig config;
config.qubits = 3;
config.pairs = 12;
auto pairs = qcanopy::testgen::generate_batch(config);
qcanopy::testgen::write_dataset({"ir", "generated/library-run"}, config, pairs);
// write_ir(stream, pairs[0].a); read_ir(stream); write_qasm(stream, pairs[0].b);
```

`generate_pair(config,i)` agrees with batch element i, including group prefix and
shuffled label. `apply_equivalent_transform`, `inverse_circuit`, `common_prefix`,
`normalized_copy`, `measure_sharing`, and the capped dense reference are also
available. Generation and serialization errors throw standard exceptions. Pass
the generated, unmodified cases and matching config to `write_dataset`; modifying
circuits afterward invalidates their construction/measurement metadata.

The normal CTest suite has no Python dependency. For an optional independent check
of matching small-case runs, install the pinned parser in a Python environment:

```sh
python -m pip install "openqasm3[parser]==1.0.1"
python tools/validate_exports.py --ir generated/ir-run --qasm generated/qasm-run
```

This tool uses the OpenQASM parser's AST, independently computes complete
operators in Python, checks exact export phase, verifies all labels up to one
global phase, and compares format-independent manifest records. It requires
1..5 qubits and handles only emitted operations. It does not execute GPU
simulation or numerically verify large-qubit cases.

## Validation in this workspace

The implementation was built in Release mode with LLVM-MinGW/Clang 23.1.5, CMake
3.31.6, and Ninja 1.11.1. The existing system MinGW GCC 8.1 had incomplete Windows
`<filesystem>` support, so a checksum-verified portable compiler was placed under
ignored `.tools/`. No global toolchain settings were changed. The working build
is `build-llvm`; its executables use static C++ runtime linking.

From this checkout in PowerShell, the prepared build can be reused directly:

```powershell
.\.tools\cmake\data\bin\cmake.exe --build build-llvm --config Release
.\.tools\cmake\data\bin\ctest.exe --test-dir build-llvm -C Release --output-on-failure
.\build-llvm\src\generator\qcanopy_generate.exe --format ir --seed 759 --output generated\new-ir
.\build-llvm\src\generator\qcanopy_generate.exe --format qasm --seed 759 --output generated\new-qasm
```

All seven CTest tests passed. Matching `generated/ir-run` and `generated/qasm-run`
samples were generated using the 4-qubit, 100-base-gate, 32-pair, 0.75-prefix,
group-size-4 examples above: each has 16 equivalent and 16 non-equivalent pairs,
with all 32 pairs numerically verified. `openqasm3[parser]==1.0.1` parsed all 64
QASM files; the independent Python evaluator confirmed exact export semantics,
all 32 labels, and matching cross-format metadata. That run covered every sampled
named gate kind (identity is covered separately in the exporter unit test).
