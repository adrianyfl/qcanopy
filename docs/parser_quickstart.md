# Parser quick start

The parser reads the QASM 3.0 files written by `qcanopy_generate --format qasm`
and returns each circuit as a `Circuit` or a `CircuitDAG`, ready for the rewrite
passes, the simulator, and the equivalence checker. This page covers the
command-line check, the C++ calls, the accepted input, and how parsed circuits
relate to the generator's IR files.

## 1. Build

The parser builds with the rest of the repository:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

It needs only the C++17 standard library. It has been built and tested with
GCC 13 on Linux (WSL); other toolchains are untested.

## 2. Generate a dataset, then parse it

```sh
./build/src/generator/qcanopy_generate --format qasm --output generated/quickstart-qasm --seed 759

./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/case-000000.a.qasm
./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/*.qasm
```

`qcanopy_parse FILE...` prints one line per file, and a summary when given
several:

```text
generated/quickstart-qasm/cases/case-000000.a.qasm: qubits=4 gates=100 global_phase=3.926990816987241
...
parsed 64 of 64 files
```

A file that fails is reported on stderr and the remaining files are still parsed.

| Exit status | Meaning |
| --- | --- |
| `0` | Every file parsed |
| `1` | At least one file could not be read or parsed |
| `2` | No file was given |

On Windows the executable sits under `src\parser` in the same layout as
`qcanopy_generate`; see step 2 of the [generator quick start](generator_quickstart.md#2-generate-test-cases).

## 3. Call the parser from C++

Include `qcanopy/parser.hpp` and link `qcanopy::parser`.

| Function | Input | Returns |
| --- | --- | --- |
| `parse_qasm(text)` | source text | `Circuit` |
| `parse_qasm_file(path)` | file path | `Circuit` |
| `parse_qasm_dag(text)` | source text | `CircuitDAG` |
| `parse_qasm_file_dag(path)` | file path | `CircuitDAG` |

Gates are kept in file order. For `cx` the first operand is the control, as in
`ir.hpp`. The returned circuit always passes `validate()`.

This program reads one generated pair and checks it for equivalence:

```cpp
#include "qcanopy/equiv.hpp"
#include "qcanopy/parser.hpp"
#include "qcanopy/stimulus.hpp"

#include <cstdio>
#include <exception>

int main(int argc, char** argv) {
    if (argc != 3) { std::fprintf(stderr, "usage: check_pair A.qasm B.qasm\n"); return 2; }
    try {
        const qcanopy::Circuit a = qcanopy::parse_qasm_file(argv[1]);
        const qcanopy::Circuit b = qcanopy::parse_qasm_file(argv[2]);

        const qcanopy::Stimulus inputs = qcanopy::haar_random_states(a.n, 8, 759);
        const qcanopy::EquivResult r = qcanopy::EquivalenceChecker().check(a, b, inputs);
        std::printf("%s\n", r.equivalent ? "equivalent" : "non_equivalent");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    return 0;
}
```

```cmake
add_executable(check_pair check_pair.cpp)
target_link_libraries(check_pair PRIVATE qcanopy::parser qcanopy::equiv)
```

Run on the default dataset above, its verdict matched the manifest's
`expected_relation` for all 32 pairs. For the rewrite passes, read each file with
`parse_qasm_file_dag` and pass the result to `normalize()` and `peel()`.

## 4. Accepted input

The parser accepts what the generator writes, plus comments and arithmetic:

- The version header `OPENQASM 3.0;` (or `3`), required and first.
- `include "stdgates.inc";`, required before any gate other than `U`.
- One register, `qubit[n] name;`, with 1 to 32 qubits.
- `gphase(a);`, any number of times.
- The gates `id h x y z s sdg t tdg rx ry rz p U cx cz swap`, applied to
  `name[i]`.
- Arguments built from numbers, `pi`, parentheses, unary `-`, and `+ - * /`.
  All arithmetic is real, so `1/2` is 0.5.
- `//` and `/* */` comments.

Everything else is rejected with an error, including:

- QASM 2.0 files (`OPENQASM 2.0;`, `qreg`, `qelib1.inc`).
- Any other gate, such as `sx` or `ccx`, and `gate` definitions.
- A second register, an unsized `qubit q;`, or more than 32 qubits.
- `measure`, `reset`, `barrier`, and classical declarations.
- `**`, function calls, and constants other than `pi`.

## 5. Parsed circuits and the generator's IR files

A circuit parsed from a `.qasm` file has the same unitary, global phase included,
as the generator's circuit, but not always the same gate list. The QASM writer
has no `sx` or `sxdg`, so those two gates come back in a different form:

| In the generator's circuit | In the `.qasm` file and the parsed circuit |
| --- | --- |
| `sx` | `rx(pi/2)`, with `pi/4` added to the global phase |
| `sxdg` | `rx(-pi/2)`, with `-pi/4` added to the global phase |
| every other gate | unchanged |

Three consequences:

- **Gate lists differ from the IR dataset.** For the default dataset (seed 759)
  generated in both formats, none of the 64 parsed circuits had the same raw gate
  list as its `.ir.json` counterpart.
- **`normalize()` removes the difference.** After `normalize()`, all 64 had the
  same gate list as their IR counterpart. Compare a parsed circuit with an IR
  circuit only after normalizing both.
- **`global_phase` is not the file's `gphase` value.** The file's `U` carries a
  phase the IR's `U` does not, so the parser adds `theta/2` per `U` gate. The
  printed phase also differs from the IR file's by `pi/4` per `sx` and `-pi/4`
  per `sxdg`.

## 6. Errors

Text that is not accepted throws `qcanopy::ParseError`, a
`std::invalid_argument` with `line()` and `column()`. Its message reads
`line:column: message`, or `file:line:column: message` when the source was a
file. A file that cannot be opened throws `std::runtime_error` instead.

| Message | Cause |
| --- | --- |
| `unsupported version '2.0'; only 3.0 is supported` | A QASM 2.0 file |
| `gate 'h' is not defined; include "stdgates.inc" first` | Missing include line |
| `unknown gate or statement 'ccx'` | A gate or statement outside the list in step 4 |
| `expected ';' but found 'q'` | Usually a missing comma between operands |
| `qubit index 4 is out of range for register 'q' of size 4` | Indices start at 0 |
| `unknown constant 'tau'` | Only `pi` is defined |
| `cannot open <path>` | The file does not exist or is unreadable |

## 7. Tests

| Test | Covers |
| --- | --- |
| `test_lexer` | Tokens, comments, positions |
| `test_parser` | Every gate, phases, arithmetic, and rejected inputs with their error positions |
| `test_parser_roundtrip` | Generator output written as QASM and read back, for every gate set, both angle modes, and 1 to 32 qubits; full unitaries are compared up to 5 qubits |
