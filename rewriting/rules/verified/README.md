# Verified local equivalence library

Default: `classes.phase.json`, with **1,356 classes / 3,516 members**, exact up to
one overall global phase. No distance tolerance or approximate resynthesis is
used. Each member is a `qdag/1` graph in a `qclass-library/1` catalog.

Other profiles: `classes.exact.json` (1,235 classes / 3,201 members, exact unitary)
and `classes.phase.small.json` (172 classes / 458 members, at most two qubits and
four gates per member). They are alternative profiles, not disjoint rule sets.

## Sources

- [Quartz](https://github.com/quantum-compiler/quartz), commit
  `2a5892efd6dab9ce763510abb131a2615e24205c`: Clifford+T, Nam, and Rigetti ECCs.
- [SQIR / VOQC](https://github.com/inQWIRE/SQIR), commit
  `134fb942265091baa3e3972f78f72da41827351f`: adapted basic equivalence lemmas.
- [Qiskit](https://github.com/Qiskit/qiskit), commit
  `da047b1fc0742316bfc2d09fc4d8e5407fa3e544`: supported standard decompositions.

Upstream validation differs across projects. Every retained local relation was
also checked with our exact symbolic matrix implementation; this is not a claim
that we reran upstream Coq proofs. Snapshots and upstream licenses are in
`sources/`; URLs and SHA256 hashes are in `sources/manifest.json`. Per-class source
records are in `provenance.json`, with counts in `report.json` and original local
verification results in `verification.json`.

## Reproduce

From the QCanopy repository root, using Python 3.9+ and no third-party packages:

```sh
# Check source/artifact hashes and all three exported profiles.
python3 rewriting/rules/verified/verify.py
# Rebuild the deduplicated profiles from the frozen source snapshots, offline.
python3 rewriting/rules/verified/build.py
```

`exact.py` compares full matrices using rational coefficients and cyclotomic
Laurent polynomials, including symbolic parameters. It checks gate unitarity and
exact equality or proportionality of matrices. No random-state test or floating
point tolerance determines catalog admission. `catalog.md` lists the classes;
`excluded.jsonl` records unsupported source relations.

Gate conventions: q0 is the low bit; CX ports are `[control, target]`; U uses the
Qiskit U(theta, phi, lambda) convention. Constants such as `"pi/2"` are exact
strings, not automatically equal to floating-point approximations. Templates can
bind parameters, but the fixed matcher does not solve general parameter
arithmetic. The separate `inverse.identity/v1` schema supports its declared exact
inverse relations without adding infinite catalog members.

The phase profile does not retain phase corrections needed to reproduce identical
complex amplitudes. Do not use its class tag alone as a phase-aware simulator
replacement. The catalog is finite and not claimed complete.
