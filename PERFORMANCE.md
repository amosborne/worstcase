# Performance investigation

Run the bounded benchmark with Poetry:

```sh
poetry run python scripts/benchmark.py --repeat 3 --samples 128 --depth 32
poetry run python scripts/benchmark.py --profile
```

It separates construction from EV, RSS, and MC evaluation, reports root evaluation
counts, and includes plain scalars, Pint, and unyt (Astropy is optional). Timings are
observations, not test thresholds. Shared graphs retain four independent inputs;
increasing graph depth must not increase the number of corner combinations.

## Findings and changes

- Profiling the initial expression implementation found assignment inference dominated
  construction. A 32-level shared graph (68 operations) took about 285 ms; direct
  `executing.Source` lookup reduced the focused construction probe to about 1.1 ms.
  Avoiding `varname`'s full-stack inspection also removes a runtime dependency.
  Names remain best-effort; stack walking is bounded and frames are not retained.
- Copying unyt quantities with `deepcopy` copied whole unit registries. Native
  shallow quantity copying owns the numeric storage while sharing registry metadata.
  Mutation tests protect both parameter and result ownership.
- Known native scalar arithmetic reads operands without redundant defensive copies.
  User callbacks and custom backend subclasses still receive copies. The focused
  shared-MC probe improved from 36.77 to 23.14 ms for scalars and 351.40 to 240.43 ms
  for Pint. This is not a claim of uniform speedup across backends or workloads.
- MC lower bounds and widths are computed once per evaluation, rather than copied
  and subtracted for every sample.
- All methods traverse a compiled postorder for a run. Shared nodes execute once
  per input assignment, without graph rebuilding, intermediate range collapse, or
  result caching across user calls.

Measurements were made locally on Windows with Python 3.12.8, NumPy 2.4.6,
Pint 0.25.3, and unyt 3.1.0, using medians and excluding import/warmup costs.
Use repeated measurements on the same machine to assess future changes. Native
quantity arithmetic and finite-value validation remain substantial costs and are
retained for correctness.
