#!/usr/bin/env python3
"""Reproducible source import, exact algebraic validation, and class construction."""
import ast
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import heapq
from itertools import permutations
import json
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from prefix_classifier import ClassLibrary
from exact import VARIABLES, compare, gate_matrix, unitary


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def hashed(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def variables(steps):
    return sorted({p[1:] for _, _, ps in steps for p in ps if isinstance(p, str) and p.startswith("$")})


def normalize_gate(op, qs, params=()):
    op = {"i": "id", "phase": "p", "cphase": "cp", "u1": "p", "u3": "u"}.get(op.lower(), op.lower())
    params = tuple(params)
    if op == "u2": op, params = "u", ("pi/2",) + params
    if op in ("rx1", "rx3", "ry1", "ry3"):
        op, params = op[:2], ("pi/2" if op[-1] == "1" else "-pi/2",)
    gate_matrix(op, params)  # Fail closed on unknown gate conventions/parameters.
    return op, tuple(qs), params


def structure(steps):
    """Canonical topological order: eliminate only independent gate-list order."""
    prev, successors, degrees = {}, defaultdict(set), [0] * len(steps)
    labels = [canonical(g) for g in steps]
    for i, (_, qs, _) in enumerate(steps):
        parents = {prev[q] for q in qs if q in prev}
        degrees[i] = len(parents)
        for j in parents: successors[j].add(i)
        for q in qs: prev[q] = i
    ready = [(labels[i], i) for i, d in enumerate(degrees) if not d]
    heapq.heapify(ready); order = []
    while ready:
        _, i = heapq.heappop(ready); order.append(steps[i])
        for j in successors[i]:
            degrees[j] -= 1
            if not degrees[j]: heapq.heappush(ready, (labels[j], j))
    return tuple(order)


def boundary(steps, reverse=False):
    seen, result = set(), {}
    indices = range(len(steps)-1, -1, -1) if reverse else range(len(steps))
    for i in indices:
        gate = steps[i]
        if not seen.intersection(gate[1]): result[gate] = i
        seen.update(gate[1])
    return result


def minimal_pair(a, b):
    a, b = list(a), list(b)
    removed = 0
    while True:
        changed = False
        for reverse in (False, True):
            aa, bb = boundary(a, reverse), boundary(b, reverse)
            common = aa.keys() & bb.keys()
            if common:
                g = min(common, key=canonical)
                a.pop(aa[g]); b.pop(bb[g]); removed += 1; changed = True
                break
        if not changed: break
    if not a and not b: return None
    wires = sorted({q for g in a + b for q in g[1]})
    params = sorted(set(variables(a)) | set(variables(b)))
    if set(variables(a)) != set(variables(b)):
        raise ValueError("parameter_disappears_on_one_side")
    if len(wires) > 3 or len(params) > 3:
        raise ValueError("width_or_parameter_count")
    best = None
    # Simultaneously rename both sides, never independently rename their ports.
    for qorder in permutations(range(len(wires))):
        qm = dict(zip(wires, qorder))
        for porder in permutations(VARIABLES[:len(params)]):
            pm = dict(zip(params, porder))
            pair = []
            for side in (a, b):
                seq = [(op, tuple(qm[q] for q in qs),
                        tuple("$" + pm[p[1:]] if isinstance(p, str) and p.startswith("$") else p for p in ps))
                       for op, qs, ps in side]
                used = {q for _, qs, _ in seq for q in qs}
                # v1 has no empty members or unused ports. Explicit identities
                # carry unchanged wires without altering the operator.
                seq += [("id", (q,), ()) for q in range(len(wires)) if q not in used]
                pair.append(structure(seq))
            pair.sort(key=canonical)
            key = canonical(pair)
            if best is None or key < best[0]: best = key, tuple(pair)
    return len(wires), best[1][0], best[1][1], removed


def dag(n, steps):
    nodes, edges, last = [], [], {}
    for i, (op, qs, ps) in enumerate(steps):
        gid = f"g{i}"
        nodes.append({"id": gid, "op": op, "qubits": list(qs),
                      "params": [{"var": p[1:]} if isinstance(p, str) and p.startswith("$") else p for p in ps]})
        for q in qs:
            if q in last: edges.append({"source": last[q], "target": gid, "qubit": q})
            last[q] = gid
    return {"schema": "qdag/1", "num_qubits": n, "global_phase": 0, "nodes": nodes, "edges": edges}


def expression(member):
    nodes = member["dag"]["nodes"]
    return " → ".join(g["op"] + ("(" + ",".join(p["var"] if isinstance(p, dict) else str(p)
                                                   for p in g["params"]) + ")" if g["params"] else "") +
                      "[" + ",".join("q" + str(q) for q in g["qubits"]) + "]" for g in nodes)


def renamed(steps, qorder, pmap):
    return structure(tuple((op, tuple(qorder[q] for q in qs),
                            tuple("$" + pmap[p[1:]] if isinstance(p, str) and p.startswith("$") else p for p in ps))
                           for op, qs, ps in steps))


def split_args(text):
    level, start, result = 0, 0, []
    for i, ch in enumerate(text):
        if ch in "([{": level += 1
        elif ch in ")]}": level -= 1
        elif ch == "," and level == 0:
            result.append(text[start:i].strip()); start = i + 1
    if text[start:].strip(): result.append(text[start:].strip())
    return result


def rust_array(text):
    text = text.strip()
    if not text.startswith("&[") or not text.endswith("]"):
        raise ValueError("unsupported_rust_array")
    return split_args(text[2:-1])


def arithmetic(text):
    """Read PI-only rational arithmetic without evaluating downloaded code."""
    def visit(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return Fraction(str(node.value)), False
        if isinstance(node, ast.Name) and node.id == "PI": return Fraction(1), True
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            a, pi = visit(node.operand); return -a, pi
        if isinstance(node, ast.BinOp):
            a, ap = visit(node.left); b, bp = visit(node.right)
            if isinstance(node.op, ast.Div) and not bp: return a / b, ap
            if isinstance(node.op, ast.Mult) and not (ap and bp): return a * b, ap or bp
            if isinstance(node.op, (ast.Add, ast.Sub)) and ap == bp:
                return (a + b if isinstance(node.op, ast.Add) else a - b), ap
        raise ValueError("unsupported_parameter_expression")
    value, pi = visit(ast.parse(text.strip(), mode="eval").body)
    if not pi:
        if value == 0: return 0
        raise ValueError("non_pi_numeric_angle")
    if not value: return 0
    numerator = ("-" if value < 0 else "") + (str(abs(value.numerator)) + "*" if abs(value.numerator) != 1 else "")
    return numerator + "pi" + ("/" + str(value.denominator) if value.denominator != 1 else "")


def rust_param(text):
    m = re.fullmatch(r"Param::ParameterExpression\((\w+)\.clone\(\)\)", text.strip())
    if m:
        name = {"theta": "theta", "phi": "phi", "lam": "lam", "lambda": "lam"}.get(m[1])
        if name is None: raise ValueError("parameter_arithmetic_not_supported_by_interface")
        return "$" + name
    text = text.strip()
    if text.startswith("Param::Float(") and text.endswith(")"): text = text[len("Param::Float("):-1]
    if text.endswith(".into()"): text = text[:-7]
    try: return arithmetic(text)
    except (SyntaxError, ValueError): raise ValueError("parameter_arithmetic_not_supported_by_interface")


def quartz_relations(root, excluded, stats):
    for filename in ("quartz--eccset--Clifford_T_5_3_complete_ECC_set.json",
                     "quartz--eccset--Nam_5_3_complete_ECC_set.json",
                     "quartz--experiment--ecc_set--rigetti_4_ecc.json"):
        data = json.loads((root / filename).read_text())[1]
        stats["source_quartz_classes"] += len(data)
        for cid, circuits in data.items():
            parsed = []
            for i, (meta, ops) in enumerate(circuits):
                try:
                    seq = []
                    for op, outs, inputs in ops:
                        if op in ("add", "neg", "mult", "pi"): continue
                        qs = [int(q[1:]) for q in inputs if q.startswith("Q")]
                        if outs != ["Q" + str(q) for q in qs]: raise ValueError("unexpected_output_ports")
                        params = []
                        for p in inputs:
                            if p.startswith("Q"): continue
                            if p not in ("P0", "P1"): raise ValueError("parameter_arithmetic_not_supported_by_interface")
                            params.append("$" + VARIABLES[int(p[1:])])
                        seq.append(normalize_gate(op, qs, params))
                    parsed.append((i, int(meta[0]), tuple(seq)))
                except ValueError as exc:
                    excluded.append({"source": filename, "class": cid, "member": i, "reason": str(exc)})
            if not parsed: continue
            # The upstream class is verified projectively. Every exported
            # relation is checked again under our explicit gate conventions.
            ref_i, n, ref = parsed[0]
            for i, width, seq in parsed[1:]:
                if width != n:
                    excluded.append({"source": filename, "class": cid, "reason": "width_mismatch"}); continue
                yield n, ref, seq, {"source": filename, "class": cid, "members": [ref_i, i],
                                   "upstream": "Quartz ECC; local exact algebraic re-verification"}


def qiskit_relations(root, excluded, stats):
    filename = "qiskit--crates--transpiler--src--standard_equivalence_library.rs"
    raw = (root / filename).read_text()
    text = re.sub(r"//[^\n]*", lambda m: " " * len(m[0]), raw)
    for match in re.finditer(r"\bcreate_standard_equivalence\(", text):
        start, depth, pos = match.end(), 1, match.end()
        while depth and pos < len(text):
            if text[pos] == "(": depth += 1
            elif text[pos] == ")": depth -= 1
            pos += 1
        line = raw.count("\n", 0, match.start()) + 1
        stats["source_qiskit_calls"] += 1
        try:
            args = split_args(text[start:pos-1])
            if len(args) != 5: raise ValueError("nonliteral_equivalence_call")
            name = re.fullmatch(r"StandardGate::(\w+)", args[0])
            if not name: raise ValueError("nonliteral_gate")
            ps = tuple(rust_param(p) for p in rust_array(args[1]))
            rhs = []
            for instruction in rust_array(args[2]):
                if not instruction.startswith("("): raise ValueError("nonliteral_instruction")
                op, qs, params = split_args(instruction[1:-1])
                op = re.fullmatch(r"StandardGate::(\w+)", op)
                if not op: raise ValueError("nonliteral_gate")
                wires = []
                for q in rust_array(qs):
                    qm = re.fullmatch(r"Qubit\((\d+)\)", q)
                    if not qm: raise ValueError("nonliteral_qubit")
                    wires.append(int(qm[1]))
                rhs.append(normalize_gate(op[1], wires, tuple(rust_param(p) for p in rust_array(params))))
            n = max((q for g in rhs for q in g[1]), default=0) + 1
            lhs = (normalize_gate(name[1], tuple(range(n)), ps),)
            if len(gate_matrix(lhs[0][0], lhs[0][2])) != 1 << n: raise ValueError("gate_width_mismatch")
            yield n, lhs, tuple(rhs), {"source": filename, "line": line,
                                      "declared_global_phase": args[3],
                                      "upstream": "Qiskit standard equivalence; local exact algebraic verification"}
        except (ValueError, IndexError, ZeroDivisionError) as exc:
            excluded.append({"source": filename, "line": line, "reason": str(exc)})


def sqir_relations():
    filename = "SQIR--SQIR--Equivalences.v"
    def g(op, *qs): return normalize_gate(op, qs)
    for name, n, a, b in [
        ("X_X_id", 1, [g("x",0),g("x",0)], [g("id",0)]),
        ("H_H_id", 1, [g("h",0),g("h",0)], [g("id",0)]),
        ("CNOT_CNOT_id", 2, [g("cx",0,1),g("cx",0,1)], [g("id",0),g("id",1)]),
        ("H_comm_Z", 1, [g("h",0),g("z",0)], [g("x",0),g("h",0)]),
        ("X_comm_CNOT_control", 2, [g("x",0),g("cx",0,1)], [g("cx",0,1),g("x",0),g("x",1)]),
        ("X_comm_CNOT_target", 2, [g("x",1),g("cx",0,1)], [g("cx",0,1),g("x",1)]),
        ("CNOT_comm_CNOT_sharing_target", 3, [g("cx",0,2),g("cx",1,2)], [g("cx",1,2),g("cx",0,2)]),
        ("CNOT_comm_CNOT_sharing_control", 3, [g("cx",0,1),g("cx",0,2)], [g("cx",0,2),g("cx",0,1)]),
        ("H_swaps_CNOT", 2, [g("h",0),g("h",1),g("cx",0,1)], [g("cx",1,0),g("h",0),g("h",1)]),
    ]:
        yield n, tuple(a), tuple(b), {"source": filename, "lemma": name,
                                    "upstream": "SQIR Coq lemma; translated and independently re-verified"}


def main():
    started = time.perf_counter(); stats = Counter(); excluded = []; relations = {}
    sources = ROOT / "sources"
    manifest = json.loads((sources / "manifest.json").read_text())
    for item in manifest:
        if hashlib.sha256((sources / item["file"]).read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError("source checksum mismatch: " + item["file"])
    streams = (quartz_relations(sources, excluded, stats), qiskit_relations(sources, excluded, stats), sqir_relations())
    for stream in streams:
        for n, a, b, source in stream:
            stats["input_relations"] += 1
            verdict = compare(unitary(n, a), unitary(n, b))
            if verdict == "different":
                excluded.append({**source, "reason": "not_projectively_equal_under_exported_gate_conventions"})
                continue
            stats["accepted_raw_" + verdict] += 1
            try: minimal = minimal_pair(a, b)
            except ValueError as exc:
                excluded.append({**source, "reason": str(exc)}); continue
            if minimal is None:
                stats["trivial_after_independent_order_or_context_removal"] += 1; continue
            n, a, b, removed = minimal
            stats["common_boundary_gates_removed"] += removed
            if a == b:
                stats["trivial_after_independent_order_or_context_removal"] += 1; continue
            verdict = compare(unitary(n, a), unitary(n, b))
            if verdict == "different": raise AssertionError("normalization changed equivalence")
            key = canonical((n, a, b))
            if key in relations:
                stats["duplicate_relations_removed"] += 1
            else:
                relations[key] = {"n": n, "a": a, "b": b, "equivalence": verdict, "sources": []}
            relations[key]["sources"].append(source)
    print(f"Imported {stats['input_relations']} relations; {len(relations)} normalized unique relations.", flush=True)
    # Merge by actual projective semantics, not by source labels or hashes alone.
    groups, buckets = [], defaultdict(list)
    for r in sorted(relations.values(), key=lambda r: (r["n"], canonical(r["a"]), canonical(r["b"]))):
        n, a, b = r["n"], r["a"], r["b"]
        params = variables(a)
        matrix = unitary(n, a)
        bucket = n, tuple(params), tuple(bool(x) for x in matrix)
        group = None
        # A class's wire/parameter names are local. Compare whole classes modulo
        # one consistent renaming of BOTH sides; do not manufacture all relabelings.
        for qorder in permutations(range(n)):
            if group is not None: break
            for porder in permutations(params):
                pm = dict(zip(params, porder))
                aa = renamed(a, qorder, pm)
                mm = unitary(n, aa)
                probe = n, tuple(params), tuple(bool(x) for x in mm)
                group = next((groups[i] for i in buckets[probe]
                              if compare(mm, groups[i]["matrix"]) != "different"), None)
                if group is not None:
                    if aa != a: stats["class_wire_or_parameter_renamings"] += 1
                    a, b = aa, renamed(b, qorder, pm)
                    break
        if group is None:
            group = {"n": n, "parameters": variables(a), "matrix": matrix, "members": {}, "sources": []}
            buckets[bucket].append(len(groups)); groups.append(group)
        for seq in (a, b):
            group["members"][canonical(seq)] = seq
        group["sources"].extend(r["sources"])
    provenance = {"schema": "verified-rule-provenance/1", "classes": {}, "sources_manifest": "sources/manifest.json"}
    all_rules = []
    for mode in ("up_to_global_phase", "exact_unitary"):
        output = []
        for group in groups:
            partitions = defaultdict(list)
            for seq in group["members"].values():
                partitions[unitary(group["n"], seq) if mode == "exact_unitary" else "all"].append(seq)
            for members in partitions.values():
                if len(members) < 2: continue
                members.sort(key=lambda s: (len(s), canonical(s)))
                cid = "eq_" + hashed((group["n"], members))[:16]
                out_members = [{"id": f"m{i:03d}_" + hashed(seq)[:10], "dag": dag(group["n"], seq)} for i, seq in enumerate(members)]
                output.append({"id": cid, "num_qubits": group["n"], "parameters": group["parameters"], "members": out_members})
                checks = []
                for i, seq in enumerate(members):
                    equality = compare(unitary(group["n"], members[0]), unitary(group["n"], seq))
                    if equality == "different" or mode == "exact_unitary" and equality != "exact":
                        raise AssertionError("invalid final class")
                    checks.append({"member_id": out_members[i]["id"], "relation_to_representative": equality})
                    if i:
                        all_rules.append({"profile": mode, "class_id": cid,
                                          "lhs_member": out_members[0]["id"], "rhs_member": out_members[i]["id"],
                                          "bidirectional": True, "verification": equality})
                evidence = {canonical(s): s for s in group["sources"]}
                provenance["classes"][mode + ":" + cid] = {"sources": list(evidence.values()), "checks": checks}
        output.sort(key=lambda c: c["id"])
        data = {"schema": "qclass-library/1", "equivalence": mode, "classes": output}
        ClassLibrary.from_dict(data)
        name = "classes.phase.json" if mode == "up_to_global_phase" else "classes.exact.json"
        (ROOT / name).write_text(json.dumps(data, indent=2) + "\n")
        if mode == "up_to_global_phase":
            lines = ["# 去重后的等价类目录", "",
                     "执行顺序从左到右。`≃` 表示精确到全局相位，不是数值近似。参数对所有实数取值成立。",
                     "每个类只列一次；来源与逐成员验证见 provenance.json。", ""]
            for width in (1, 2, 3):
                lines.extend([f"## {width} qubit", ""])
                for c in output:
                    if c["num_qubits"] != width: continue
                    lines.append(f"### {c['id']} ({len(c['members'])} members)")
                    lines.append("")
                    lines.extend("- `" + expression(m) + "`" for m in c["members"])
                    lines.append("")
            (ROOT / "catalog.md").write_text("\n".join(lines) + "\n")
            small = []
            for c in output:
                if c["num_qubits"] > 2: continue
                ms = [m for m in c["members"] if len(m["dag"]["nodes"]) <= 4]
                if len(ms) >= 2: small.append({**c, "members": ms})
            (ROOT / "classes.phase.small.json").write_text(json.dumps({**data, "classes": small}, indent=2) + "\n")
            stats["small_classes"] = len(small)
            stats["small_members"] = sum(len(c["members"]) for c in small)
        stats[mode + "_classes"] = len(output)
        stats[mode + "_members"] = sum(len(c["members"]) for c in output)
        stats[mode + "_independent_bidirectional_rules"] = sum(len(c["members"]) - 1 for c in output)
        stats[mode + "_symbolic_classes"] = sum(bool(c["parameters"]) for c in output)
        print(name, len(output), 'classes;', stats[mode + "_members"], 'members', flush=True)
    (ROOT / "rules.jsonl").write_text("".join(canonical(r) + "\n" for r in all_rules))
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (ROOT / "excluded.jsonl").write_text("".join(canonical(r) + "\n" for r in excluded))
    stats["excluded_records"] = len(excluded)
    report = {"verification": "exact cyclotomic Laurent-polynomial matrices; no floating point or tolerance",
              "statistics": dict(stats), "excluded_reasons": dict(Counter(r["reason"] for r in excluded)),
              "build_seconds": time.perf_counter() - started,
              "artifacts": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in
                            ("classes.phase.json", "classes.phase.small.json", "classes.exact.json", "rules.jsonl", "catalog.md", "provenance.json", "excluded.jsonl")}}
    (ROOT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["statistics"], indent=2))


if __name__ == "__main__":
    main()
