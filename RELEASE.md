# 1.0 release checklist

Release development targets `release/1.0`. The final squash merge into `master`
requires Alex's explicit approval of the reviewed commit. Keep the release PR
in draft and do not enable automatic merging.

- Integrate PR #4's corrected optional unit support, retaining contributor credit.
- Build scalar expressions, explicit evaluation, function wrapping, formatting,
  inferred names, and dependency-preserving sensitivity studies.
- Verify numerical regressions and all five optional unit backends.
- Execute README examples, the current-sensing notebook, and installed artifacts.
- Run an adversarial review covering simplicity, single source of truth,
  ownership, correctness, Python idioms, tests, ergonomics, integrations,
  performance, reproducibility, API boundaries, packaging, and documentation.
- Fix release-blocking findings and repeat affected checks.
- Present the final diff, migration guide, validation evidence, review findings,
  and exact commit for approval before merging, tagging, or publishing.

Deferred: optimization solver, array-valued uncertainty propagation, correlated
distributions, and general NumPy integration. EV enumerates corners; MC samples
uniform ranges; RSS is an endpoint sensitivity approximation. None proves a
global bound for arbitrary nonlinear functions.
