# 1.0 release checklist

Release development targets `release/1.0`. The final squash merge into `master`
requires Alex's explicit approval of the reviewed commit. Keep the release PR
in draft and do not enable automatic merging.

- [x] Integrate PR #4's corrected optional unit support, retaining contributor credit.
- [x] Build scalar expressions, explicit evaluation, function wrapping, formatting,
  inferred names, and dependency-preserving sensitivity studies.
- [x] Verify numerical regressions and all five optional unit backends.
- [x] Execute README examples, the current-sensing notebook, and installed artifacts.
- [x] Profile construction and evaluation and implement measured improvements.
- [x] Run an adversarial review covering simplicity, single source of truth,
  ownership, correctness, Python idioms, tests, ergonomics, integrations,
  performance, reproducibility, API boundaries, packaging, and documentation.
- [x] Fix release-blocking findings and repeat affected checks.
- [ ] Alex explicitly approves the exact final release-branch commit.
- [ ] Squash merge into master, retaining contributor credit.
- [ ] Verify the merged commit before tagging or publishing.

Local final-code validation: 144 tests pass, 6 explicit backend-specific skips,
98% statement coverage. The strict dependency audit reports no known
vulnerabilities; all 48 current default-branch Dependabot alerts are reconciled
in [REVIEW.md](REVIEW.md). GitHub will close applicable alerts only after the
approved merge reaches master.

See [MIGRATION.md](MIGRATION.md), [PERFORMANCE.md](PERFORMANCE.md), and
[REVIEW.md](REVIEW.md). The release candidate is `1.0.0rc1`, built and tested from
source; it has not been tagged or published.

The final squash commit must include the original PR contributor:

```text
Co-authored-by: Reid Sox-Harris <rharris@kepler.space>
```

Issues #5 (format specifications) and #6 (inferred names) are implemented in this
release. Leave them open until the approved default-branch merge closes them.

Deferred: optimization solver, array-valued uncertainty propagation, correlated
distributions, and general NumPy integration. EV enumerates corners; MC samples
uniform ranges; RSS is an endpoint sensitivity approximation. None proves a
global bound for arbitrary nonlinear functions.
