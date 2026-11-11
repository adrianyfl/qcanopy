# Test-case generator quick start

Use `qcanopy_generate` to create circuit pairs with known expected labels for an
equivalence checker. Choose JSON IR or OpenQASM 3.0; generation happens in the
same in-memory IR for both formats. Negative pairs append an X gate to an
equivalent construction, giving their non-equivalence labels a proof.

## 1. Build the executable

Start in a clone of this repository:

```sh
git clone https://github.com/adrianyfl/qcanopy.git
cd qcanopy
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Install CMake 3.18+ and a compiler with complete C++17 library support before
building. No CUDA, Python, or third-party C++ library is needed for the generator
or normal CTest suite. The optional independent QASM validator has separate
[Python setup instructions](test_case_generator.md#library-use-and-independent-export-checks).

## 2. Generate test cases

Pick the executable for your build:

| Build | Executable from the repository root |
| --- | --- |
| Linux/macOS, single configuration | `./build/src/generator/qcanopy_generate` |
| Windows, Visual Studio Release | `.\build\src\generator\Release\qcanopy_generate.exe` |
| Windows, Ninja single configuration | `.\build\src\generator\qcanopy_generate.exe` |

For Linux/macOS:

```sh
./build/src/generator/qcanopy_generate --format ir --output generated/quickstart-ir --seed 759
./build/src/generator/qcanopy_generate --format qasm --output generated/quickstart-qasm --seed 759
```

For Windows PowerShell with a Visual Studio Release build:

```powershell
.\build\src\generator\Release\qcanopy_generate.exe --format ir --output generated\quickstart-ir --seed 759
.\build\src\generator\Release\qcanopy_generate.exe --format qasm --output generated\quickstart-qasm --seed 759
```

For Ninja on Windows, use the executable path without `Release`. A build directory
with another name also needs that name in the path.

**`--format` and `--output` are required.** Each destination must be absent or
empty. To rerun, choose a new destination. The commands above generate the same
32 pairs in two formats, with 16 equivalent and 16 non-equivalent labels in each.
They use the defaults of 4 qubits, 100 base gates, a 0.5 planted-prefix fraction,
and prefix group size 1. The two formats do not represent 64 distinct pairs.

## 3. Adjust the experiment

For a batch with 75 planted prefix gates and groups of four pairs sharing that
prefix, run this example using the executable for your platform:

```sh
./build/src/generator/qcanopy_generate --format ir --output generated/shared-ir --qubits 4 --gates 100 --pairs 32 --prefix-fraction 0.75 --prefix-group-size 4 --equivalent-fraction 0.5 --gate-set mixed --angle-mode mixed --seed 759
```

| Option | Default | What it controls |
| --- | --- | --- |
| `--qubits N` | `4` | Circuit size, from 1 through 32 qubits |
| `--gates L` | `100` | Base gate count; transformed B may have a different length |
| `--pairs B` | `32` | Positive number of pairs |
| `--equivalent-fraction E` | `0.5` | Equivalent count is `floor(E * B + 0.5)`, followed by a seeded label shuffle |
| `--prefix-fraction R` | `0.5` | Plant `floor(R * L)` initial gates in common |
| `--prefix-group-size G` | `1` | Consecutive groups of G pairs share the planted prefix |
| `--gate-set SET` | `mixed` | Base gate pool: `clifford`, `clifford-t`, `parameterized`, or `mixed` |
| `--angle-mode MODE` | `mixed` | `discrete` multiples of pi/4, `continuous`, or `mixed` |
| `--seed UINT64` | `759` | Decimal unsigned master seed |
| `--verify MODE` | `small` | Dense verification for small circuits, or `none` |

Fractions must be in [0,1]. Keep generation options and the seed identical when
comparing formats. The same generator version reproduces the same constructions;
standard-library floating-point details can differ in the last bit across platforms.
Use `--help` for the command's option summary.

## 4. Read the generated files

For the default 32-pair IR run:

```text
generated/quickstart-ir/
  manifest.jsonl             1 configuration header + 32 case records
  cases/
    case-000000.a.ir.json    circuit A
    case-000000.b.ir.json    circuit B
    ...                     64 circuit files in total
```

QASM mode has the same layout with `.a.qasm` and `.b.qasm` circuit files.
It targets **OpenQASM 3.0**, including `stdgates.inc` and `gphase`; a QASM 2-only
reader is not compatible. IR mode writes versioned JSON that the library's
`read_ir` function can reconstruct. Both formats preserve circuit global phase,
while the equivalence labels allow one global phase for the entire operator.

For each case, read these manifest fields:

- `expected_relation`: `equivalent` or `non_equivalent`.
- `files.a` and `files.b`: circuit paths relative to the dataset directory.
- `transform` and `negative_mutation`: construction details and label evidence.
- `seeds`: master, prefix, suffix, and transformation seeds, as decimal strings.
- `sharing.raw` and `sharing.normalized`: gate counts and common-prefix metrics.
- `verification.status`: numerical verification outcome.

The common-prefix fraction divides matching initial gates by the shorter sequence
length. Normalization can change order and length, so requested, raw, and
normalized sharing need not match. A negative pair can have a fraction of 1 when
the shorter circuit is a complete prefix of the other. Output files retain the
raw circuits; normalized statistics come from copies.

## 5. Check verification and common errors

| Verification status | Meaning |
| --- | --- |
| `passed` | Full-operator label and normalization checks passed, at N <= 5 |
| `skipped_size_limit` | N > 5; structural checks and construction proofs still apply |
| `skipped_by_request` | The run used `--verify none`; structural checks still apply |

Generation supports 1-32 qubits, but dense operators require memory proportional
to `4^N`, so numerical verification stops at five qubits. Large batches also retain
all generated circuits in memory; size them to available RAM.

- **Executable not found:** check the platform/configuration path in step 2 and
  confirm that the build completed.
- **Compiler error in `<filesystem>`:** use a toolchain with a complete C++17
  standard library. Some older Windows MinGW distributions are incomplete.
- **Destination is not empty:** choose another output directory. Files are never
  overwritten automatically.
- **Invalid argument:** check option spelling, integer ranges, and fractions with
  `--help`. Success returns exit code 0; failures return nonzero.
- **Only `manifest.jsonl.partial` exists:** the export did not finish. Inspect the
  reported error and use a fresh destination; the complete manifest is published last.

Generated data under `generated/` stays local and Git-ignored. For exact formats,
all option semantics, and library usage, see the [generator reference](test_case_generator.md).
