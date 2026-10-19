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
- [x] Add the manual publishing workflow and protect the GitHub `pypi` environment.
- [ ] Alex explicitly approves the exact final release-branch commit.
- [ ] Squash merge into master, retaining contributor credit.
- [ ] Verify the merged commit before tagging or publishing.
- [ ] Configure the PyPI trusted publisher described below.
- [ ] Run the separate publishing workflow and approve its environment gate.

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

## Publishing through GitHub

Poetry owns dependency resolution and builds. The official PyPA publishing action
uploads those verified artifacts using PyPI trusted publishing and generates
attestations. No long-lived PyPI token or literal `poetry publish` is needed.
Ordinary CI never publishes; neither does merging or pushing a tag.

One-time setup, while signed in as a maintainer of the existing PyPI project:

1. Open [worstcase publishing settings](https://pypi.org/manage/project/worstcase/settings/publishing/).
2. Add a GitHub trusted publisher with owner `amosborne`, repository `worstcase`,
   workflow filename `publish.yml`, and environment `pypi`.
3. Confirm the GitHub `pypi` environment requires approval by `amosborne`, allows
   only the `master` branch, and disallows administrator bypass. Self-approval is
   allowed because this project has one maintainer. Do not add a PyPI token secret.

For each release, after explicit approval and the squash merge:

1. Choose the version in `pyproject.toml` before the merge. It is currently
   `1.0.0rc1`; a final `1.0.0` requires a version change and fresh review/checks.
   The README targets the stable release with `pip install worstcase`. That
   command normally skips prereleases: publishing only `1.0.0rc1` will not make
   it select the new API over an existing stable 0.x release. Use
   `pip install --pre worstcase` when deliberately testing a published candidate;
   publish stable `1.0.0` for the README's ordinary installation path.
2. Wait for successful **push CI on the exact master commit**. Create and push
   its matching tag, such as `v1.0.0rc1`, only after release approval.
3. In GitHub Actions, manually run **Publish to PyPI**, select branch `master`,
   and enter that tag. The tag must equal `v` plus Poetry's version and point to
   the exact selected commit. Branch runs other than master cannot publish.
4. Inspect the build job and its wheel/sdist artifact, then approve deployment
   to `pypi`. The upload job downloads only this run's verified artifacts; it
   does not check out source or execute the build with publishing credentials.
5. Verify the version and provenance on PyPI. Uploads are irreversible; an
   existing version is not silently skipped or overwritten. Fixes need a new
   version and another reviewed release.

The publishing workflow cannot run from the release branch. PyPI configuration
and a real upload are deliberately untested until maintainer setup and release
approval. See [PyPI's trusted-publishing guide](https://docs.pypi.org/trusted-publishers/using-a-publisher/).
