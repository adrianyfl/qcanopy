import argparse
import json
from pathlib import Path
import sys

from . import InputError, Limits, OverlapLimits, classify, load_dag, load_library


def main(argv=None):
    parser = argparse.ArgumentParser(description="Match a fixed class library and link boundary occurrences.")
    parser.add_argument("--dag", required=True, type=Path)
    parser.add_argument("--classes", required=True, type=Path)
    parser.add_argument("--output", "-o", type=Path, help="result JSON (default: stdout)")
    parser.add_argument("--progress", action="store_true", help="round progress on stderr")
    parser.add_argument("--enable-overlap", action="store_true",
                        help="opt in to bounded partial-RHS overlap search")
    for name, value in Limits().__dict__.items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(value), default=value)
    for name, value in OverlapLimits().__dict__.items():
        parser.add_argument("--overlap-" + name.replace("_", "-"),
                            type=type(value), default=None)
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.resolve() in (args.dag.resolve(), args.classes.resolve()):
            raise InputError("output must not overwrite an input")
        dag, library = load_dag(args.dag), load_library(args.classes)
        limits = Limits(**{k: getattr(args, k) for k in Limits().__dict__})
        overrides = {k: getattr(args, "overlap_" + k) for k in OverlapLimits().__dict__
                     if getattr(args, "overlap_" + k) is not None}
        if overrides and not args.enable_overlap:
            raise InputError("--overlap-* options require --enable-overlap")
        callback = (lambda row: print(json.dumps(row), file=sys.stderr)) if args.progress else None
        result = classify(dag, library, limits, on_round=callback,
                          enable_overlap=args.enable_overlap,
                          overlap_limits=OverlapLimits(**overrides) if args.enable_overlap else None)
        encoded = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.write_text(encoded, encoding="utf-8")
        else:
            print(encoded, end="")
        return 0 if (result["complete_for_bounds"] and
                     result.get("overlap", {}).get("complete_for_bounds", True)) else 2
    except (InputError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
