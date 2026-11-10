#!/usr/bin/env python3
"""Optional independent small-case export validation (requires openqasm3[parser]).

Parse QASM with the upstream OpenQASM parser, interpret this generator's straight-line
gate subset independently, and compare entire operators, including global phase,
against matching JSON IR. This is validation tooling, not a production QASM reader.
"""
import argparse
import cmath
import json
import math
from pathlib import Path

import openqasm3
from openqasm3 import ast


def expression(node):
    if isinstance(node, (ast.FloatLiteral, ast.IntegerLiteral)):
        return node.value
    if isinstance(node, ast.UnaryExpression) and node.op.name == "-":
        return -expression(node.expression)
    raise ValueError(f"Unsupported expression in generated QASM: {node!r}")


def operand(node):
    assert node.name.name == "q" and len(node.indices) == 1 and len(node.indices[0]) == 1
    return int(expression(node.indices[0][0]))


def read_qasm(path):
    program = openqasm3.parse(path.read_text(encoding="utf-8"))
    assert program.version == "3.0"
    n, phase, gates = None, 0.0, []
    for node in program.statements:
        if isinstance(node, ast.Include):
            assert node.filename == "stdgates.inc"
        elif isinstance(node, ast.QubitDeclaration):
            assert n is None and node.qubit.name == "q"
            n = int(expression(node.size))
        elif isinstance(node, ast.QuantumPhase):
            assert not node.modifiers and not node.qubits
            phase += expression(node.argument)
        elif isinstance(node, ast.QuantumGate):
            assert not node.modifiers
            gates.append({"name": node.name.name, "params": [expression(a) for a in node.arguments],
                          "qubits": [operand(q) for q in node.qubits]})
        else:
            raise ValueError(f"Unexpected generated statement: {node!r}")
    assert n is not None
    return {"n": n, "global_phase": phase, "gates": gates}


def matrix(g):
    name, params = g["name"], g["params"]
    exp = lambda t: cmath.exp(1j*t)
    r = math.sqrt(0.5)
    fixed = {
        "id": [1, 0, 0, 1], "h": [r, r, r, -r], "x": [0, 1, 1, 0],
        "y": [0, -1j, 1j, 0], "z": [1, 0, 0, -1],
        "s": [1, 0, 0, 1j], "sdg": [1, 0, 0, -1j],
        "t": [1, 0, 0, exp(math.pi/4)], "tdg": [1, 0, 0, exp(-math.pi/4)],
        "sx": [.5+.5j, .5-.5j, .5-.5j, .5+.5j],
        "sxdg": [.5-.5j, .5+.5j, .5+.5j, .5-.5j],
        "cx": [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
        "cz": [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1],
        "swap": [1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1],
    }
    if name in fixed:
        assert not params
        return fixed[name]
    if name in ("unitary_1q", "unitary_2q"):
        return [complex(*z) for z in g["matrix"]]
    if name in ("u", "U"):
        theta, phi, lam = params
        c, s = math.cos(theta/2), math.sin(theta/2)
        # OpenQASM 3's built-in U has exp(i*theta/2) relative to this IR's U.
        phase = exp(theta/2) if name == "U" else 1
        return [phase*x for x in [c, -exp(lam)*s, exp(phi)*s, exp(phi+lam)*c]]
    theta, = params
    c, s = math.cos(theta/2), math.sin(theta/2)
    return {
        "rx": [c, -1j*s, -1j*s, c], "ry": [c, -s, s, c],
        "rz": [exp(-theta/2), 0, 0, exp(theta/2)], "p": [1, 0, 0, exp(theta)],
    }[name]


def unitary(c):
    n = c["n"]
    if not 1 <= n <= 5:
        raise ValueError("Independent dense validation is limited to 1..5 qubits")
    dim = 1 << n
    # Row-major operator; each row is a separate list. Apply gates to every column.
    u = [[cmath.exp(1j*c["global_phase"]) if row == col else 0j for col in range(dim)] for row in range(dim)]
    for g in c["gates"]:
        qs, m = g["qubits"], matrix(g)
        local = 1 << len(qs)
        assert len(m) == local*local and len(set(qs)) == len(qs)
        mask = sum(1 << q for q in qs)
        for base in range(dim):
            if base & mask:
                continue
            ix = [base | sum(((j >> k) & 1) << q for k, q in enumerate(qs)) for j in range(local)]
            old = [u[i][:] for i in ix]
            for row, index in enumerate(ix):
                u[index] = [sum(m[row*local+k]*old[k][col] for k in range(local)) for col in range(dim)]
    return [z for row in u for z in row]


def same(a, b, ignore_phase=False):
    assert len(a) == len(b)
    phase = 1
    if ignore_phase:
        k = max(range(len(b)), key=lambda i: abs(b[i]))
        phase = a[k]/b[k]
        if abs(abs(phase)-1) > 1e-8:
            return False
        phase /= abs(phase)
    return max(abs(x-phase*y) for x, y in zip(a, b)) < 1e-8


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ir", type=Path, required=True)
    p.add_argument("--qasm", type=Path, required=True)
    args = p.parse_args()
    ir = [json.loads(line) for line in (args.ir/"manifest.jsonl").read_text().splitlines()]
    qasm = [json.loads(line) for line in (args.qasm/"manifest.jsonl").read_text().splitlines()]
    assert len(ir) == len(qasm) and ir[0]["config"] == qasm[0]["config"]
    circuits, kinds = 0, set()
    for a, b in zip(ir[1:], qasm[1:]):
        assert {k: v for k, v in a.items() if k != "files"} == {k: v for k, v in b.items() if k != "files"}
        operators = []
        for side in ("a", "b"):
            ci = json.loads((args.ir/a["files"][side]).read_text())
            cq = read_qasm(args.qasm/b["files"][side])
            assert ci["n"] == cq["n"]
            ui, uq = unitary(ci), unitary(cq)
            assert same(ui, uq), f"Export semantics differ: {a['case_id']} {side}"
            operators.append(ui)
            kinds.update(g["name"] for g in ci["gates"])
            circuits += 1
        assert same(*operators, ignore_phase=True) == (a["expected_relation"] == "equivalent"), a["case_id"]
    print(f"Parsed and independently verified {circuits} QASM circuits, {len(ir)-1} labels, matching cross-format metadata.")
    print("Gate coverage: " + ", ".join(sorted(kinds)))


if __name__ == "__main__":
    main()
