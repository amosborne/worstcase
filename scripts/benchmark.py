"""Bounded scalar DAG benchmark; run with `poetry run python scripts/benchmark.py`.

Times are observations, never pass/fail thresholds. Optional backend imports and
one warmup are excluded. Use --profile for cumulative cProfile evaluation costs.
"""

import argparse
import cProfile
import platform
import pstats
import statistics
from time import perf_counter

from worstcase import param


def unit_for(name):
    if name == "scalar":
        return 1.0
    if name == "pint":
        import pint

        return pint.UnitRegistry().volt
    if name == "unyt":
        import unyt

        return unyt.V
    if name == "astropy":
        from astropy import units

        return units.V
    raise ValueError(name)


def graph(shape, unit, depth):
    primitives = [
        param.bytol_rel(float(i + 1) * unit, 0.01, tag=f"p{i}") for i in range(4)
    ]
    expression = sum(primitives)
    if shape == "deep":
        for _ in range(depth):
            expression = expression + primitives[0]
    elif shape == "shared":
        for _ in range(depth):
            expression = (expression + expression) / 2
    return expression


def measured(call, repeat):
    call()  # Exclude first-use caches and name-inference warmup.
    samples = []
    result = None
    for _ in range(repeat):
        start = perf_counter()
        result = call()
        samples.append(perf_counter() - start)
    return statistics.median(samples), result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backends",
        nargs="+",
        default=["scalar", "pint", "unyt"],
        choices=["scalar", "pint", "unyt", "astropy"],
    )
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--depth", type=int, default=32)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()
    if not (1 <= args.samples <= 10000 and 1 <= args.depth <= 256):
        parser.error("samples must be 1..10000 and depth must be 1..256")
    if not 1 <= args.repeat <= 10:
        parser.error("repeat must be 1..10")

    print(f"Python {platform.python_version()} on {platform.system()}")
    print(f"median of {args.repeat}; depth={args.depth}; MC n={args.samples}")
    print("backend shape phase milliseconds evaluations")
    for backend in args.backends:
        unit = unit_for(backend)
        for shape in ("small", "deep", "shared"):
            seconds, expression = measured(
                lambda: graph(shape, unit, args.depth), args.repeat
            )
            print(f"{backend} {shape} construct {seconds * 1000:.3f} -")
            for method in ("ev", "rss", "mc"):
                options = {"n": args.samples, "seed": 42} if method == "mc" else {}

                def call():
                    return expression.derive(method=method, **options)

                seconds, result = measured(call, args.repeat)
                print(
                    f"{backend} {shape} {method} {seconds * 1000:.3f} "
                    f"{result.n_evaluations}"
                )
            if args.profile and shape == "shared":
                profiler = cProfile.Profile()
                profiler.runcall(call)
                print(f"Cumulative MC profile: {backend}/{shape}")
                pstats.Stats(profiler).strip_dirs().sort_stats(
                    "cumulative"
                ).print_stats(18)


if __name__ == "__main__":
    main()
