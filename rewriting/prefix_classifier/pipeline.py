"""Team-facing DAG pipelines; the matcher and input DAGs remain unchanged."""
from pathlib import Path

from .engine import classify
from .model import load_library
from .pair import prove_pair
from .swaps import remove_swaps


DEFAULT_LIBRARY = Path(__file__).resolve().parents[1] / "rules/verified/classes.phase.json"


def load_default_library():
    """Load once and reuse for a batch; no network or simulator is needed."""
    return load_library(DEFAULT_LIBRARY)


def analyze_dag(dag, library, limits=None, *, enable_overlap=False, overlap_limits=None):
    """Remove SWAPs, then return core-DAG classification plus output mapping.

    Class labels, gate IDs, cuts and overlap candidates refer to core_dag.
    This does not apply overlap candidates or extract a simulator reuse plan.
    """
    prepared = remove_swaps(dag)
    result = classify(prepared.core_dag, library, limits,
                      enable_overlap=enable_overlap, overlap_limits=overlap_limits)
    return {"schema": "qrewrite-analysis/1", "swap_removal": prepared.to_dict(),
            "classification": result}


def compare_dags(left, right, library, limits=None, overlap_limits=None, *,
                 enable_overlap=False, base_seconds=6, overlap_seconds=8,
                 max_alternatives=32, max_residual_gates=128, enable_inverse=True):
    """Prove a pair modulo global phase, retaining SWAPs if layouts differ.

    The underlying prove_pair API assumes common physical wire coordinates.
    Two SWAP-free cores alone are insufficient when their final mappings differ.
    """
    a, b = remove_swaps(left), remove_swaps(right)
    use_cores = a.output_wire_to_core == b.output_wire_to_core
    proof = prove_pair(a.core_dag if use_cores else left,
                       b.core_dag if use_cores else right, library, limits, overlap_limits,
                       enable_overlap=enable_overlap, base_seconds=base_seconds,
                       overlap_seconds=overlap_seconds, max_alternatives=max_alternatives,
                       max_residual_gates=max_residual_gates, enable_inverse=enable_inverse)
    return {"schema": "qrewrite-pair/1",
            "swap_removal": {"left": a.to_dict(), "right": b.to_dict()},
            "layout": {"cores_used": use_cores,
                       "policy": "Compare cores only with equal final mappings; otherwise compare original DAGs."},
            "proof": proof}
