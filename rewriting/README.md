# DAG Rewriting Interface

Standalone Python 3.9+ module; standard library only. It identifies equivalent
boundary regions without changing the input DAG or expanding the fixed class
library. No simulator, CUDA, C++ build, or package installation is required.

## Workflow

```text
Input DAG → remove SWAPs and retain mapping → rewriting analysis → results
                                             ├─ fixed equivalence templates
                                             ├─ parametric identity relations
                                             └─ optional overlap search
```

Exact inverse relations are identity equivalences within this rewriting model,
not a separate user-facing optimization service. Current routing of these rules
is described below; not every rule type is yet exposed by both analysis paths.

## Run from the QCanopy repository root

```sh
# Remove explicit SWAPs and retain the output permutation.
python3 rewriting/run.py remove-swaps --dag rewriting/examples/swaps.json \
  -o generated/rewrite/swaps.json

# Full pipeline: remove SWAPs, then classify the core DAG. Overlap is OFF.
python3 rewriting/run.py classify --dag rewriting/examples/circuit.json \
  -o generated/rewrite/off.json

# ON: a small fixture that demonstrates partial-overlap candidates.
python3 rewriting/run.py classify --dag rewriting/examples/overlap-circuit.json \
  --classes rewriting/examples/overlap-classes.json \
  --enable-overlap -o generated/rewrite/on.json

# Compare two circuits: common-boundary cancellation and residual proof.
python3 rewriting/run.py pair --left rewriting/examples/pair-a.json \
  --right rewriting/examples/pair-b.json --enable-overlap \
  -o generated/rewrite/pair.json
```

The default is `rules/verified/classes.phase.json` (1,356 classes / 3,516 members).
Use `--classes PATH` to supply another library. Paths to inputs/outputs are relative
to the working directory; the default library path is relative to this module.
Omit `-o` to print JSON. Exit codes: **0** completed classification / proved pair;
**2** partial classification or inconclusive pair; **1** input/IO error.
A completed classification does **not** mean all gates matched.

For generator JSON IR v1, add `--input-format qcanopy-ir` to any command. This
adapter constructs per-wire DAG edges and preserves parameters; it does not
parse QASM. Generic matrix gates are retained as opaque, exact entry tokens,
without inferring named-gate equivalence or inverse relations. There are no
in-process bindings to the C++ `qcanopy::CircuitDAG` yet: exchange JSON or use the
Python DAG directly. This module does not change the existing C++ rewrite passes.

## Python API

From the repository root:

```python
from rewriting import (
    load_dag, load_default_library, analyze_dag, compare_dags,
    remove_swaps, classify, Limits, OverlapLimits,
)

library = load_default_library()  # Load once and reuse across circuits.
dag = load_dag("rewriting/examples/circuit.json")
result = analyze_dag(dag, library, enable_overlap=True)
core = result["swap_removal"]["core_dag"]         # JSON-compatible DAG
mapping = result["swap_removal"]["output_wire_to_core"]
labels = result["classification"]

# If the caller already removed SWAPs, call the matcher directly.
prepared = remove_swaps(dag)
labels = classify(prepared.core_dag, library, enable_overlap=False)

# Optional pair checker: use equivalence rules to prove the remaining difference.
a = load_dag("rewriting/examples/pair-a.json")
b = load_dag("rewriting/examples/pair-b.json")
proof = compare_dags(a, b, library, enable_overlap=True)
print(proof["proof"]["verdict"])
```

## Output

- **Single circuit:** `swap_removal` contains the core DAG and final output-wire
  mapping; `classification` contains class/member matches, gate tags, links,
  boundary regions, and shared-composition groups. ON adds proof-backed local
  `rewrite_candidates`, separately from original-DAG tags.
- **Circuit pair:** both SWAP-removal results, a mapping policy, proof steps,
  residual DAGs, and verdict `equivalent_up_to_global_phase` or `inconclusive`.
  Inconclusive does not mean non-equivalent. If output mappings differ, the
  checker retains the original SWAPs rather than comparing cores unsafely.

Keep the mapping with the core: `output_wire_to_core[q] = k` means original output
wire q is core wire k. Qubit 0 is the least significant bit. Groups can overlap;
matching every gate alone is not a whole-circuit equivalence proof.
See [JSON schemas, mapping, and output fields](docs/interface.md).

## How rewriting works

- **Fixed equivalence templates:** match member DAGs with ordered ports,
  dependencies, and exact parameter bindings; link compatible occurrences into
  larger boundary regions. Class/member definitions stay fixed while occurrence
  and link records grow within budgets.
- **Parametric identity relations:** recognize `G G† = I` and reversed inverse
  groups using `inverse.identity/v1`. A specialized matcher checks exact parameter
  relations and per-wire adjacency, avoiding enumeration of every angle or group.
  These are equivalence-rule proofs, not an independent optimization model.
- **Optional overlap:** temporarily substitute a member in a small window, then
  check a second rule that partially overlaps the replacement. OFF skips this
  search; ON enables it. Arbitrary-depth nesting and full saturation are not implemented.

**Current integration:** the pair checker consumes identity-rule certificates in
both OFF and ON, alongside structural boundary cancellation and fixed-library
proofs. Single-circuit classification does not yet expose these parametric
identity matches in its tags/candidates. Unifying that output remains pending.

## Current scope

**Simulator integration remains future work:** no cost-based extraction, reuse
scheduling, GPU execution, or measured simulation speedup. The default library
allows global phase differences without tracking relative phase factors; class
matches therefore do not promise bit-identical complex state vectors. Suffix
operator matches also do not establish reusable zero-input states.

## Files and checks

| Path | Purpose |
|---|---|
| `run.py`, `prefix_classifier/pipeline.py` | CLI and composed team-facing APIs |
| `prefix_classifier/model.py`, `qcanopy.py` | Validated DAG/library schemas and generator IR adapter |
| `prefix_classifier/swaps.py` | SWAP elision and output mapping |
| `prefix_classifier/engine.py`, `inverse.py`, `inverse_schemas.json`, `overlap.py` | Equivalence-rule matching: fixed templates, parametric identities, and optional overlap |
| `prefix_classifier/pair.py` | Consume rule proofs to cancel common boundaries and check residuals |
| `rules/verified/` | Fixed libraries, provenance, verification tools, and source licenses |
| `examples/`, `tests/` | Small inputs and correctness/integration tests |

```sh
(cd rewriting && python3 -m unittest discover -s tests -v)
# Optional: independently recheck every member of the default library.
python3 rewriting/rules/verified/verify.py \
  --library rewriting/rules/verified/classes.phase.json
```

See [rule sources and verification](rules/verified/README.md). Benchmark datasets,
profiling outputs, plots, and machine-specific build files are intentionally not
part of this module.
