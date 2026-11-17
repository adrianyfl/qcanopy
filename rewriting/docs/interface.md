# Interface contracts

Input `qdag/1` uses stable gate IDs and **immediate per-wire dependency edges**:

```json
{
  "schema": "qdag/1", "num_qubits": 2, "global_phase": 0,
  "nodes": [
    {"id": "a", "op": "h", "qubits": [0], "params": []},
    {"id": "b", "op": "cx", "qubits": [0, 1], "params": []}
  ],
  "edges": [{"source": "a", "target": "b", "qubit": 0}]
}
```

Operands are ordered: `cx` uses `[control, target]`. Qubit 0 is the least
significant bit. Parameters are finite numbers or exact strings; fixed-library
matching does not round numeric angles to symbolic `pi` expressions. Inputs must
be unitary circuits; measurement, reset, and classical control are unsupported.

| Output | Contents |
|---|---|
| `remove-swaps`: `qswap-removal/1` | `core_dag`, `output_wire_to_core`, removed gate IDs, original fingerprint |
| `classify`: `qrewrite-analysis/1` | The SWAP-removal result plus `classification` below |
| `pair`: `qrewrite-pair/1` | Both SWAP-removal results, layout policy, proof trace, and residual DAGs |

The single-circuit `classification` contains:

- `matches`: class/member IDs, matched core gate IDs, parameter bindings, ordered
  qubit map, and prefix/suffix flags.
- `regions`, `links`, `composition_groups`: occurrence-backed combinations and
  groups with the same composition signature; `annotations` maps gates to matches
  and regions. Groups can overlap; they are candidates, not a selected reuse plan.
- `rewrite_candidates` **only with ON**: original/alternative window DAGs,
  mappings, and two-step rule witnesses. These proposals are separate from the
  original-DAG tags and are not automatically applied.
- `unmatched_gate_ids`, `stop_reason`, `complete_for_bounds`, and search statistics.
  ON has its own completion status under `overlap`; inspect both statuses.

Keep the mapping with the core: `output_wire_to_core[q] = k` means original output
wire **q** is core wire **k**. For example, `[1, 0]` exchanges the two output bits.
`SwapRemovalResult.output_index_to_core(i)` converts a little-endian basis index;
`original_amplitude[i] = core_amplitude[output_index_to_core(i)]`.
For pair checking, cores are compared only when their final mappings agree;
otherwise the original DAGs, including SWAPs, are checked.

Pair verdicts are `equivalent_up_to_global_phase` or `inconclusive`.
Inconclusive does not establish non-equivalence. Matching a single DAG is not a
whole-circuit equivalence verdict; matching every gate is insufficient by itself.


Pair proof gate IDs refer to the cores only if `layout.cores_used` is true; otherwise
they refer to the original inputs. Low-level matcher and SWAP CLIs remain available
as `python3 -m prefix_classifier` and `python3 -m prefix_classifier.swap_cli` when
run from the `rewriting/` directory.

## Identity rules within rewriting

`inverse.identity/v1` represents parametric identity equivalences, including
supported gate/inverse pairs and reversed inverse groups. The specialized matcher
is an implementation detail of equivalence-rule recognition, not a separate
required pipeline step or CLI command. Pair proofs consume its certificates in
both OFF and ON. Single-circuit classification does not yet include these
parametric matches in `matches`, `annotations`, or `rewrite_candidates`; ordinary
identity templates already present in the fixed library can still match there.

## Search limits

`analyze_dag` accepts `limits=Limits(...)` and `overlap_limits=OverlapLimits(...)`;
only supply overlap limits when overlap is enabled. The CLI exposes base and
overlap pass time limits (`--time-limit-seconds`, `--overlap-time-limit-seconds`).
`pair` also exposes `--base-seconds`, `--overlap-seconds`, and
`--max-residual-gates`. Defaults are bounded searches, not exhaustive saturation.

