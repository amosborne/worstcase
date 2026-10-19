# 1.0 adversarial review

Reviewed code: [88f4d7a](https://github.com/amosborne/worstcase/commit/88f4d7ab26f98e3d1f71f8979344a3e0f29db1c6).

Reviewed the cumulative implementation after the initial feature work, including
the original code, PR #4 integration, the expression redesign, packaging, and the
profiling changes. This is a code review with executable counterexamples, not a
claim that arbitrary numerical models have guaranteed bounds. The performance
investigation also received a separate subagent profiling pass.

## Findings and disposition

| Severity | Finding | Resolution and evidence |
|---|---|---|
| High | Pint can enclose an expression as a quantity magnitude and forward `.derive()` to it, losing outer units | Use Pint's documented upcast-type registration for concrete expression types. Test quantity-left arithmetic for every backend, plus Pint imported before and after primitive construction. No library methods are replaced. |
| High | Magnitude extraction in RSS can turn a 3–5 V range into a kilovolt upper bound | Normalize outputs to nominal units and perform RSS in native quantity arithmetic. Scaled-output regression runs across compatible backends. |
| High | Backend approximate equality can freeze small forallpeople tolerances | Compare normalized numeric bound magnitudes exactly to identify varying inputs. A 1e-12 relative-tolerance regression verifies three EV evaluations and distinct bounds. |
| Medium | Deep copying quantities can duplicate entire unit registries | Copy native quantity storage, share registry metadata, and test mutations of input quantities, returned values, provenance values, callbacks, and custom subclasses. |
| Medium | NumPy's minimum signed integer overflows under `abs(nom)` | Apply the validated floating tolerance before taking absolute deviation; verify ordered bounds for `np.int64(-(2**63))`. |
| Medium | Concrete decorated calls copied arbitrary user objects | Detect symbolic arguments before capturing data; verify a concrete object with a deliberately failing copy operation passes through unchanged. |
| Medium | Equality-based LaTeX symmetry can hide visible high-precision differences | Compare formatted deviations at the requested precision. Verify asymmetric output at `.14gL` on every backend. |
| Medium | Compound python-quantities units used internal text in LaTeX | Delegate unit rendering to its native dimensionality LaTeX formatter; verify compound units and handcalcs rendering. |
| Medium | Strict auditing failed on the unpublished editable project | Audit exact installed dependencies, excluding only this project; fail on unresolved dependency audits. Do not dismiss or ignore advisory IDs. |
| Low | Redundant node keyword plumbing, display overrides, and callback copies obscured ownership | Keep keyword/template binding in the function wrapper, use the display mixin consistently, and isolate operands once at the callback boundary. Mutation regressions remain enabled. |
| Low | Duplicate push/PR CI runs wasted execution | Limit push runs to master, use PR checks for development, and cancel superseded runs. |
| Low | A prerelease installation command could install the older published API | Document installing this unpublished candidate from its checkout. |

Pint integration follows its documented extension point:
https://pint.readthedocs.io/en/stable/user/numpy.html#technical-commentary

## Review angles

| Angle | Assessment |
|---|---|
| Simplicity | One postorder traversal and evaluator serve all methods. No runtime NetworkX/pyDOE, legacy method-specific decorators, or bundled notebook server. Runtime dependencies are NumPy and executing. |
| Single source of truth | Unit normalization and validation live in the adapter boundary; source parameter identity drives every solver and sensitivity study. Function binding has one implementation. |
| Ownership | Nodes are structurally immutable. Quantity storage and result/provenance access are isolated. Registry metadata is deliberately shared. No cross-call result cache or global RNG state. |
| Correctness | Shared variables cancel; equal-looking independent parameters remain independent. Inputs and evaluated outputs must be finite scalar values with compatible dimensions. Failure notes identify the operation and endpoint/sample assignment. |
| Python idioms | Native operators, functools.wraps, normal signature/default/positional-only behavior, numeric format specifications, NumPy ufunc dispatch, and explicit exceptions. Name inference is bounded, best-effort, and does not retain frames. |
| Test quality | Analytical expectations, original circuit values, observed sampling strata, mutation counterexamples, callback counts, and real backend behavior supplement ordinary examples. No probabilistic reruns or time-threshold assertions. Inapplicable backend cases are explicit skips. |
| Ergonomics | Expressions are readable algebra; bounds require explicit evaluation. Results display without rerunning. Wrapped functions handle internal arrays/solvers and static configuration. Names have an explicit override and a safe empty fallback. |
| Integrations | All five backends receive common conformance tests. handcalcs renders actual results. Optional backends are not imported by the core; installed-artifact tests run without them. Raw quantity-left arithmetic either preserves the expression or remains a documented unsupported operation. |
| Performance | Profile-backed inference and copying changes, cached MC ranges, shared-node evaluation counts, bounded repeatable benchmarks, and no repeated graph construction within a solver. See PERFORMANCE.md. |
| Reproducibility | Stable traversal/parameter order and local seeded generators. EV/MC extrema record input assignments; RSS deliberately has no fictitious joint bound assignment. |
| API boundaries | Reject unsupported arrays, symbolic truth tests, result-as-expression arithmetic, incompatible backends/registries, and nonfinite outputs. User callbacks must be deterministic and meaningful for the input domain. |
| Packaging/docs | Poetry lock, clean wheel/sdist installs, executed README and notebook examples, migration guidance, supported-platform CI, and a strict dependency audit. No deployment or publishing workflow. |

## Validation evidence

Local final-code regression run on Python 3.12.8/Windows: **144 passed, 6 skipped,
98% statement coverage**. Skips concern physical dimensions on plain scalars,
backends without retained scaled magnitudes, and immutable quantities that cannot
be mutated. They do not skip mandatory backend imports.

The final hardening commit passed Linux and Windows on Python 3.11–3.14, notebook
execution, quality checks, artifact installation, and dependency auditing in
[CI run 39183817729](https://github.com/amosborne/worstcase/actions/runs/39183817729).
Both rebuilt artifacts also passed local isolated installs without unit packages.
Final integration status is tracked in the draft release PR.

## Dependabot reconciliation

All 48 open default-branch alerts were compared with the release Poetry lock:

| Package | Alerts | Release disposition |
|---|---:|---|
| anyio | 2 | Removed with obsolete frontend/server stack |
| black | 1 | Updated to 26.10.1 |
| bleach | 3 | Removed |
| jupyter-server | 1 | Removed |
| jupyterlab | 9 | Removed |
| mistune | 10 | Removed |
| pillow | 13 | Removed |
| soupsieve | 2 | Removed |
| tornado | 7 | Updated to 6.5.10 |

The final Poetry environment, including development and examples dependencies,
reported **no known vulnerabilities** through the strict audit script. GitHub's
default-branch alerts remain open until the approved release reaches master;
none has been dismissed merely to hide it.

## Remaining limits and release gate

No known blocking review finding remains after the corrections and passing
regressions. EV is corner analysis, MC is a sampled envelope, and RSS is an endpoint
approximation. None proves global bounds for arbitrary nonlinear functions.
Quantity-left syntax is not uniformly supported across libraries; explicit
`param.constant(...)` is the portable form. Use linear units for uncertainty models.

Optimization, uncertain arrays, correlated distributions, and broader NumPy
interoperability remain deferred. Source/AST availability can prevent label
inference. Arbitrary callback side effects and global/registry mutation are outside
the evaluation contract.

The release PR stays draft. Alex must explicitly approve the exact final commit
before squash merging into master; changes after approval require renewed approval.
Preserve contributor credit in the squash commit. Do not tag or publish before the
release checks and approval gate are satisfied.
