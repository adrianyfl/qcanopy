#!/usr/bin/env python3
"""Standalone JSON CLI. Run from any directory; all defaults are module-relative."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

from prefix_classifier import InputError, Limits, OverlapLimits, load_dag, load_library, remove_swaps
from prefix_classifier.pipeline import DEFAULT_LIBRARY, analyze_dag, compare_dags
from prefix_classifier.qcanopy import parse_qcanopy_ir


def positive_seconds(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("expected positive finite seconds")
    return number


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("expected positive integer")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description="SWAP removal, boundary classification, and pair proof.")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("remove-swaps", "classify", "pair"):
        cmd = commands.add_parser(name)
        if name == "pair":
            cmd.add_argument("--left", type=Path, required=True)
            cmd.add_argument("--right", type=Path, required=True)
        else:
            cmd.add_argument("--dag", type=Path, required=True)
        cmd.add_argument("--input-format", choices=("dag", "qcanopy-ir"), default="dag")
        cmd.add_argument("--output", "-o", type=Path, help="JSON output (default: stdout)")
        if name != "remove-swaps":
            cmd.add_argument("--classes", type=Path, default=DEFAULT_LIBRARY)
            cmd.add_argument("--enable-overlap", action="store_true", help="OFF by default")
            cmd.add_argument("--time-limit-seconds", type=positive_seconds, default=5.0,
                             help="base classification pass limit")
            cmd.add_argument("--overlap-time-limit-seconds", type=positive_seconds)
        if name == "pair":
            cmd.add_argument("--base-seconds", type=positive_seconds, default=6.0)
            cmd.add_argument("--overlap-seconds", type=positive_seconds, default=8.0)
            cmd.add_argument("--max-residual-gates", type=positive_int, default=128)
    args = parser.parse_args(argv)
    try:
        inputs = [args.left, args.right] if args.command == "pair" else [args.dag]
        protected = inputs + ([args.classes] if hasattr(args, "classes") else [])
        if args.output and args.output.resolve() in {p.resolve() for p in protected}:
            raise InputError("output must not overwrite an input or class library")
        dags = [load_dag(p) if args.input_format == "dag" else
                parse_qcanopy_ir(json.loads(p.read_text(encoding="utf-8")), allow_opaque_matrices=True)
                for p in inputs]
        exit_code = 0
        if args.command == "remove-swaps":
            result = remove_swaps(dags[0]).to_dict()
        else:
            if args.overlap_time_limit_seconds is not None and not args.enable_overlap:
                raise InputError("--overlap-time-limit-seconds requires --enable-overlap")
            library = load_library(args.classes)
            limits = replace(Limits(), time_limit_seconds=args.time_limit_seconds)
            overlap_limits = (replace(OverlapLimits(), time_limit_seconds=args.overlap_time_limit_seconds)
                              if args.overlap_time_limit_seconds is not None else None)
            if args.command == "classify":
                result = analyze_dag(dags[0], library, limits, enable_overlap=args.enable_overlap,
                                     overlap_limits=overlap_limits)
                c = result["classification"]
                exit_code = 0 if (c["complete_for_bounds"] and
                    c.get("overlap", {}).get("complete_for_bounds", True)) else 2
            else:
                result = compare_dags(dags[0], dags[1], library, limits, overlap_limits,
                                      enable_overlap=args.enable_overlap, base_seconds=args.base_seconds,
                                      overlap_seconds=args.overlap_seconds,
                                      max_residual_gates=args.max_residual_gates)
                exit_code = 0 if result["proof"]["verdict"] == "equivalent_up_to_global_phase" else 2
        encoded = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(encoded, encoding="utf-8")
        else:
            print(encoded, end="")
        return exit_code
    except (InputError, OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
