#!/usr/bin/env python3
"""Re-read exported files and verify every class using exact, symbolic matrices."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from prefix_classifier import ClassLibrary
from exact import ZERO, ONE, add, mul, conj, gate_matrix, unitary, compare, from_dag


def verify_gate(op, params):
    matrix = gate_matrix(op, params)
    for a in range(len(matrix)):
        for b in range(len(matrix)):
            value = ZERO
            for k in range(len(matrix)):
                value = add(value, mul(conj(matrix[k][a]), matrix[k][b]))
            if value != (ONE if a == b else ZERO):
                raise AssertionError(f"nonunitary gate definition: {op} {params}")


def verify(path):
    data = json.loads(path.read_text())
    ClassLibrary.from_dict(data)
    stats = Counter(); gates = set()
    for c in data["classes"]:
        if len(c["members"]) < 2: raise AssertionError("singleton class")
        reference = None
        for m in c["members"]:
            steps = from_dag(m["dag"])
            for op, _, params in steps:
                if (op, params) not in gates:
                    verify_gate(op, params); gates.add((op, params))
            matrix = unitary(c["num_qubits"], steps)
            if reference is None: reference = matrix
            result = compare(reference, matrix)
            if result == "different" or data["equivalence"] == "exact_unitary" and result != "exact":
                raise AssertionError(f"failed {path.name}:{c['id']}:{m['id']}")
            stats[result] += 1
        stats["classes"] += 1
    return {"file": path.name, "equivalence": data["equivalence"],
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "member_checks": dict(stats), "unitary_gate_instances_checked": len(gates)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    if not args.library:
        for item in json.loads((ROOT / "sources/manifest.json").read_text()):
            path = ROOT / "sources" / item["file"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
                raise AssertionError("source changed: " + str(path))
        report = json.loads((ROOT / "report.json").read_text())
        for name, expected in report["artifacts"].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
                raise AssertionError("artifact changed: " + name)
    paths = [args.library] if args.library else [ROOT / n for n in
             ("classes.phase.json", "classes.exact.json", "classes.phase.small.json")]
    results = [verify(p) for p in paths]
    output = {"method": "exact cyclotomic Laurent-polynomial matrix comparison and gate unitarity",
              "numeric_tolerance": None, "results": results, "elapsed_seconds": time.perf_counter() - start}
    if not args.library: (ROOT / "verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
