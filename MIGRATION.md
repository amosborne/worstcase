# Migrating to 1.0

1.0 replaces method-specific decorators and partial binding with immutable scalar
expressions and explicit evaluation. This is intentionally a breaking release.

| Before | After |
|---|---|
| `from worstcase import unit` | Create a registry/environment in your preferred unit library |
| `param.bytol(nom, tol, rel=True)` | `param.bytol_rel(nom, tol)` |
| `param.bytol(nom, tol, rel=False)` | `param.bytol_abs(nom, tol)` |
| `@derive.byev(a=A)` | Write an expression or use `@derive` and pass `A` explicitly |
| `C()` / `C.lb` | `result = expression.derive(); result.lb` |
| `@derive.bymc(..., n=1000)` | `expression.derive(method="mc", n=1000, seed=42)` |
| `@derive.byrss(...)` | `expression.derive(method="rss")` |
| `sigfig=3` | `format(result, ".3g")` |
| Mandatory explicit tags | Inferred labels with optional `tag` override |
| `C(tag="study").ss(A)` | `expression.ss(A, tag="study")` |
| Pint-only `ito` and `apply` methods | Convert concrete quantities in your own code or a wrapped function |

Expressions have no implicit bounds. Their `.nom` remains available for nominal
references and calibration formulas. Result display and access are stable snapshots.
If you need to freeze an evaluated uncertainty as a new independent input, construct
`param.byrange(result.nom, result.lb, result.ub)` explicitly; this discards original
dependencies. Prefer reusing the expression when those dependencies matter.

Choose one solver at the final expression rather than mixing methods inside a graph.
Source primitives are identity-based, not value-equal. Shared primitives remain
shared even when nested through several expressions or wrapped function calls.

Result provenance is keyed by original primitive objects. RSS bound assignments are
`None`; Monte Carlo provenance now records the observed extrema (including nominal).
Input validation raises `TypeError`/`ValueError` instead of removable assertions.

The `examples` Poetry group no longer installs JupyterLab, its formatter extension,
or Matplotlib. The repository needs notebook execution, not a bundled web frontend;
removing these dependencies also eliminates their obsolete vulnerable transitive
stack. Install a maintained frontend separately if desired.
