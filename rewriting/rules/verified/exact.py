"""Exact matrices over Q(zeta_16)[e^(i theta/2), e^(i phi/2), e^(i lam/2)]^{+/-1}.

No floating-point arithmetic or tolerances are used for acceptance. Polynomials
are immutable sparse tuples; zeta_16^8 = -1. All supported gates are unitary for
real parameters. Projective equality is checked by exact cross multiplication.
"""
from fractions import Fraction
from functools import lru_cache
import re

ZERO = ()
VARIABLES = ("theta", "phi", "lam")


def polynomial(items):
    terms = {}
    for (p, a, b, c), value in items:
        turns, p = divmod(p, 8)
        value = Fraction(value) * (-1 if turns % 2 else 1)
        key = p, a, b, c
        terms[key] = terms.get(key, 0) + value
    return tuple(sorted((k, Fraction(v)) for k, v in terms.items() if v))


def constant(value=1):
    return polynomial([((0, 0, 0, 0), value)])


ONE = constant()
I = polynomial([((4, 0, 0, 0), 1)])
HALF = constant(Fraction(1, 2))
HSCALE = polynomial([((2, 0, 0, 0), Fraction(1, 2)), ((6, 0, 0, 0), Fraction(-1, 2))])


@lru_cache(maxsize=100000)
def add(a, b):
    if not a: return b
    if not b: return a
    return polynomial(a + b)


@lru_cache(maxsize=100000)
def neg(a):
    return tuple((key, -value) for key, value in a)


@lru_cache(maxsize=200000)
def mul(a, b):
    if not a or not b: return ZERO
    if a == ONE: return b
    if b == ONE: return a
    return polynomial((tuple(x + y for x, y in zip(ka, kb)), va * vb)
                      for ka, va in a for kb, vb in b)


@lru_cache(maxsize=100000)
def conj(a):
    return polynomial((tuple(-x for x in key), value) for key, value in a)


def exp_half(angle):
    if isinstance(angle, str) and angle.startswith("$"):
        exponents = [0, 0, 0, 0]
        exponents[1 + VARIABLES.index(angle[1:])] = 1
        return polynomial([(tuple(exponents), 1)])
    if angle == 0 or angle == "0":
        return ONE
    match = re.fullmatch(r"(-?)(?:(\d+)\*)?pi(?:/(\d+))?", str(angle))
    if not match:
        raise ValueError(f"unsupported exact angle {angle!r}")
    sign, numerator, denominator = match.groups()
    power = Fraction(4 * int(numerator or 1), int(denominator or 1))
    if power.denominator != 1:
        raise ValueError(f"angle outside zeta_16 field: {angle}")
    return polynomial([((int(power) * (-1 if sign else 1), 0, 0, 0), 1)])


def identity(n):
    return tuple(tuple(ONE if r == c else ZERO for c in range(n)) for r in range(n))


@lru_cache(maxsize=None)
def gate_matrix(op, params=()):
    fixed = {
        "id": ((ONE, ZERO), (ZERO, ONE)),
        "x": ((ZERO, ONE), (ONE, ZERO)),
        "y": ((ZERO, neg(I)), (I, ZERO)),
        "z": ((ONE, ZERO), (ZERO, neg(ONE))),
        "h": ((HSCALE, HSCALE), (HSCALE, neg(HSCALE))),
    }
    if op in fixed:
        if params: raise ValueError(f"unexpected params for {op}")
        return fixed[op]
    if op in ("s", "sdg", "t", "tdg"):
        if params: raise ValueError(f"unexpected params for {op}")
        power = {"s": 4, "sdg": -4, "t": 2, "tdg": -2}[op]
        return ((ONE, ZERO), (ZERO, polynomial([((power, 0, 0, 0), 1)])))
    if op in ("sx", "sxdg"):
        if params: raise ValueError(f"unexpected params for {op}")
        a, b = mul(HALF, add(ONE, I)), mul(HALF, add(ONE, neg(I)))
        return ((a, b), (b, a)) if op == "sx" else ((conj(a), conj(b)), (conj(b), conj(a)))
    if op in ("rx", "ry", "rz", "p", "rxx", "ryy", "rzz", "rzx"):
        if len(params) != 1: raise ValueError(f"arity of {op}")
        e = exp_half(params[0]); em = conj(e)
        co = mul(HALF, add(e, em))
        si = mul(neg(I), mul(HALF, add(e, neg(em))))
        if op == "p": return ((ONE, ZERO), (ZERO, mul(e, e)))
        if op == "rz": return ((em, ZERO), (ZERO, e))
        if op == "rx": return ((co, mul(neg(I), si)), (mul(neg(I), si), co))
        if op == "ry": return ((co, neg(si)), (si, co))
        left, right = gate_matrix(op[1]), gate_matrix(op[2])
        # Operand 0 is the least significant bit; RZX has Z on operand 0.
        pauli = tuple(tuple(mul(left[r & 1][c & 1], right[r >> 1][c >> 1])
                            for c in range(4)) for r in range(4))
        return tuple(tuple(add(co if r == c else ZERO, mul(mul(neg(I), si), pauli[r][c]))
                           for c in range(4)) for r in range(4))
    if op == "u":
        if len(params) != 3: raise ValueError("arity of u")
        a = exp_half(params[0]); co = mul(HALF, add(a, conj(a)))
        si = mul(neg(I), mul(HALF, add(a, neg(conj(a)))))
        ph, la = exp_half(params[1]), exp_half(params[2])
        ph, la = mul(ph, ph), mul(la, la)
        return ((co, neg(mul(la, si))), (mul(ph, si), mul(mul(ph, la), co)))
    controls = {"cx": "x", "cy": "y", "cz": "z", "ch": "h", "cs": "s",
                "csdg": "sdg", "csx": "sx", "cp": "p", "crx": "rx", "cry": "ry", "crz": "rz"}
    if op in controls:
        base = gate_matrix(controls[op], params)
        return tuple(tuple((ONE if r == c else ZERO) if not c & 1 else
                           (base[r >> 1][c >> 1] if r & 1 else ZERO)
                           for c in range(4)) for r in range(4))
    if op in ("swap", "iswap", "iswapdg", "ccx", "cswap"):
        if params: raise ValueError(f"unexpected params for {op}")
        size = 8 if op in ("ccx", "cswap") else 4
        result = [[ZERO] * size for _ in range(size)]
        for col in range(size):
            row, factor = col, ONE
            if op == "ccx":
                row = col ^ 4 if col & 3 == 3 else col
            elif op == "cswap":
                row = col ^ 6 if col & 1 and ((col >> 1) ^ (col >> 2)) & 1 else col
            elif col in (1, 2):
                row = col ^ 3
                factor = I if op == "iswap" else neg(I) if op == "iswapdg" else ONE
            result[row][col] = factor
        return tuple(map(tuple, result))
    raise ValueError(f"unsupported gate {op}")


@lru_cache(maxsize=16000)
def unitary(n, steps):
    dim = 1 << n
    columns = []
    for basis in range(dim):
        state = [ONE if i == basis else ZERO for i in range(dim)]
        for op, qs, params in steps:
            matrix = gate_matrix(op, params)
            if len(matrix) != 1 << len(qs) or len(set(qs)) != len(qs):
                raise ValueError(f"invalid gate operands: {op} {qs}")
            out = [ZERO] * dim
            mask = sum(1 << q for q in qs)
            for index, amplitude in enumerate(state):
                if not amplitude: continue
                col = sum(((index >> q) & 1) << port for port, q in enumerate(qs))
                for row in range(len(matrix)):
                    value = matrix[row][col]
                    if not value: continue
                    dest = (index & ~mask) | sum(((row >> port) & 1) << q for port, q in enumerate(qs))
                    out[dest] = add(out[dest], mul(value, amplitude))
            state = out
        columns.append(state)
    return tuple(columns[c][r] for r in range(dim) for c in range(dim))


def compare(a, b):
    """Return exact/projective/different. Inputs are complete unitary matrices.

    Cross products prove scalar proportionality over the Laurent polynomial
    fraction field. Analytic continuation covers isolated pivot zeros; unitarity
    ensures the proportionality scalar has magnitude one for real parameters.
    """
    if a == b: return "exact"
    if len(a) != len(b): return "different"
    pivot = next((i for i, x in enumerate(a) if x), None)
    if pivot is None or not b[pivot]: return "different"
    if all(mul(x, b[pivot]) == mul(y, a[pivot]) for x, y in zip(a, b)):
        return "projective"
    return "different"


def from_dag(data):
    from prefix_classifier import CircuitDAG
    dag = CircuitDAG.from_dict(data, pattern=True)
    if dag.global_phase != 0: raise ValueError("explicit member phase must be zero")
    return tuple((dag.by_id[g].op, dag.by_id[g].qubits,
                  tuple("$" + p.name if hasattr(p, "name") else p for p in dag.by_id[g].params))
                 for g in dag.order)
