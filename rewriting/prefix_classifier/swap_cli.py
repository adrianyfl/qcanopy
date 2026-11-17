"""Standalone DAG SWAP-removal CLI; classification remains a separate API."""
import argparse
import json
from pathlib import Path
import sys

from . import InputError, load_dag, remove_swaps


def main(argv=None):
    parser = argparse.ArgumentParser(description="Remove explicit SWAP gates from a qdag/1 DAG.")
    parser.add_argument("--dag", required=True, type=Path)
    parser.add_argument("--output", "-o", type=Path,
                        help="qswap-removal/1 JSON with core DAG and mapping (default: stdout)")
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.resolve() == args.dag.resolve():
            raise InputError("output must not overwrite an input")
        result = remove_swaps(load_dag(args.dag))
        encoded = json.dumps(result.to_dict(), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.write_text(encoded, encoding="utf-8")
        else:
            print(encoded, end="")
        return 0
    except (InputError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
