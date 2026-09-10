# Closure fixtures

Each directory is an isolated downstream application used to verify one profile
in `architecture/closure-profiles.toml` (ADR 0036).

## Rules

- Every fixture has an EMPTY `[workspace]` table. It must be its own workspace
  root, otherwise workspace feature unification masks what a real consumer
  resolves and the measurement becomes meaningless.
- Depend on leaf packages by relative path with `default-features = false`.
- `Cargo.lock` is gitignored per fixture; the checker resolves with `--offline`.
- Exercise real behaviour, not just symbol names. A fixture that only mentions a
  type can pass while the package is unusable.
- Adding a dependency here changes a declared compatibility promise. Update the
  profile deliberately and regenerate the closure docs.

## Verify

```bash
cargo xtask architecture closure check
cargo xtask architecture closure explain <profile>
bash scripts/probe_closure_gate.sh
```

## Build output is never tracked

Each consumer is a real cargo project, so running one creates a local
`target/`. Those artifacts must never enter git.

The root `.gitignore` had `/target`, which is anchored and matches only
the repo-root directory. Nested `tests/consumers/*/target/` was therefore
NOT ignored, and 1427 artifact files (70 MB) were committed by accident.

The rule is now `**/target/`, which matches at any depth. Verify with:

```bash
git check-ignore -v tests/consumers/2d-curves/target
```

A `.gitignore` rule does NOT apply to a file already tracked, so fixing
the pattern alone changes nothing: the path must be untracked first with
`git rm -r --cached <dir>` (index only -- leaves files on disk).

