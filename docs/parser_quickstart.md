# Parser quick start

The parser reads OpenQASM circuit files and returns each circuit as a `Circuit`
or a `CircuitDAG`, ready for the rewrite passes, the simulator, and the
equivalence checker. It reads two dialects:

- **OpenQASM 3.0**, as written by `qcanopy_generate --format qasm`.
- **OpenQASM 2.0**, as used by benchmark sets such as Quokka#.

This page covers the command-line check, the C++ calls, the accepted input, how
parsed circuits relate to their source, and errors.

## 1. Build

The parser builds with the rest of the repository:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQCANOPY_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

It needs only the C++17 standard library. It has been built and tested with
GCC 13 on Linux (WSL); other toolchains are untested.

## 2. Parse files from the command line

Generated OpenQASM 3.0 files:

```sh
./build/src/generator/qcanopy_generate --format qasm --output generated/quickstart-qasm --seed 759

./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/case-000000.a.qasm
./build/src/parser/qcanopy_parse generated/quickstart-qasm/cases/*.qasm
```

OpenQASM 2.0 benchmark files are parsed the same way; the file's version header
selects the dialect:

```sh
./build/src/parser/qcanopy_parse path/to/benchmark/*.qasm
```

`qcanopy_parse FILE...` prints one line per file, and a summary when given
several:

```text
generated/quickstart-qasm/cases/case-000000.a.qasm: qubits=4 gates=100 global_phase=3.926990816987241
...
parsed 64 of 64 files
```

A file that fails is reported on stderr with its position and the reason, and
the remaining files are still parsed:

```text
ccx_gate_0.qasm:4:10: expected ';' but found 'q'
```

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

This program reads one pair of files and checks it for equivalence:

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

Run on the default generated dataset above, its verdict matched the manifest's
`expected_relation` for all 32 pairs. For the rewrite passes, read each file with
`parse_qasm_file_dag` and pass the result to `normalize()` and `peel()`.

## 4. Accepted input

The version header is required, comes first, and selects the dialect.

| | OpenQASM 3.0 | OpenQASM 2.0 |
| --- | --- | --- |
| Header | `OPENQASM 3.0;` (or `3`) | `OPENQASM 2.0;` |
| Library | `include "stdgates.inc";` | `include "qelib1.inc";` |
| Registers | One `qubit[n] name;` | Any number of `qreg name[n];`; `creg name[n];` is ignored |
| Built-in gates | `U` | `U`, `CX` |
| Library gates | `id h x y z s sdg t tdg rx ry rz p cx cz swap` | All of `qelib1.inc`: `u1 u2 u3 id x y z h s sdg t tdg rx ry rz cx cz cy ch ccx crz cu1 cu3` |
| Other statements | `gphase(a);` | `gate` definitions |

In both dialects:

- A library gate needs the include before it; built-in gates do not.
- A gate is applied to `name[i]` operands separated by commas, all different.
- A circuit has 1 to 32 qubits in total.
- Arguments are built from numbers, `pi`, parentheses, unary `-`, and `+ - * /`.
  All arithmetic is real, so `1/2` is 0.5. (OpenQASM 3 divides two integers as
  integers, giving 0; that rule is not implemented.)
- `//` and `/* */` comments are allowed.

With several `qreg`s, qubits are numbered end to end in declaration order: after
`qreg a[2]; qreg b[3];`, `b[0]` is qubit 2 of the circuit.

A 2.0 gate definition has the form

```text
gate name(p, ...) a, b, ... { g(e, ...) a, b; ... }
```

Its body may call `U`, `CX`, and gates defined earlier, on its own arguments.
Each call is expanded in place, so the circuit holds only one- and two-qubit
gates. A name can be defined only once.

Everything else is rejected with an error, including:

- More than 32 qubits.
- `measure`, `reset`, `barrier`, `if`, `opaque`, and whole-register operands such
  as `h q;`.
- Gates that are in neither the library nor the file, such as `mcx` and `csdg`.
- Operands written without commas, such as `cx q[0] q[1];`.
- `**`, `^`, function calls, and constants other than `pi`.
- In 3.0: a second register, an unsized `qubit q;`, `gate` definitions, and
  gates outside the list above, such as `sx` and `ccx`.

## 5. How a parsed circuit relates to its source

### OpenQASM 3.0 files from the generator

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

### OpenQASM 2.0 files

- **Gate counts can exceed the statement count.** `id x y z h s sdg t tdg rx rz
  cx cz` become one IR gate each. The other library gates and the file's own
  gates are expanded from their definitions; `ccx` becomes 15 gates.
- **`rz` is read as the IR's `RZ`,** as Qiskit reads it. `qelib1.inc` itself
  defines `rz` as `u1`, which differs by a global phase.
- **`ch` carries a global phase** of `e^{i pi/4}` relative to a controlled-H,
  from its `qelib1.inc` definition.
- **`cu3` has one line added to its `qelib1.inc` definition,**
  `u1((lambda+phi)/2) c;`, as in Qiskit's copy of the library. Without it the
  gate is not a controlled-U3.

Checked against Qiskit's OpenQASM 2 loader on 977 Quokka# files and on 600
generated programs with gate definitions, the simulated output states were
identical, global phase included.

## 6. Errors

Text that is not accepted throws `qcanopy::ParseError`, a
`std::invalid_argument` with `line()` and `column()`. Its message reads
`line:column: message`, or `file:line:column: message` when the source was a
file. A file that cannot be opened throws `std::runtime_error` instead. Parsing
stops at the first error in a file.

| Message | Cause |
| --- | --- |
| `unsupported version '4.0'; only 2.0 and 3.0 are supported` | A header other than 2.0 or 3.0 |
| `gate 'h' is not defined; include "qelib1.inc" first` | Missing include line |
| `unknown gate or statement 'csdg'` | A gate or statement the dialect does not have |
| `circuit would have more than 32 qubits` | The registers exceed the IR's limit |
| `expected ';' but found 'q'` | Usually a missing comma between operands |
| `qubit index 4 is out of range for register 'q' of size 4` | Indices start at 0 |
| `gate 'cx' needs distinct qubits` | The same qubit used twice in one gate |
| `gate 'h' is already defined` | A `gate` definition reusing an existing name |
| `unknown constant 'tau'` | Only `pi` is defined |
| `cannot open <path>` | The file does not exist or is unreadable |

## 7. Tests

| Test | Covers |
| --- | --- |
| `test_lexer` | Tokens, comments, positions |
| `test_parser` | Both dialects: every gate, phases, arithmetic, registers, `qelib1.inc` gates and gate definitions against their matrices, and rejected inputs with their error positions |
| `test_parser_roundtrip` | Generator output written as OpenQASM 3.0 and read back, for every gate set, both angle modes, and 1 to 32 qubits; full unitaries are compared up to 5 qubits |
