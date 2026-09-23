# architecture/

Machine-checked architecture declarations.

`closure-profiles.toml` declares minimal downstream dependency closures. Each
profile names an isolated consumer fixture under `tests/consumers/`, the exact
internal packages that must be present, and the packages that must be absent.

Verify with:

```bash
cargo xtask architecture closure check
cargo xtask architecture closure explain <profile>
```

A closure change is an API change. Update `expected_internal` deliberately and
record the reason in an ADR — never to silence a failing gate.

## Capability ledger (what to build next)

`capability-ledger.toml` grades every geometry capability found in OCCT and
CGAL against Axiolid (implemented, narrow or absent), with evidence paths,
reference packages at pinned commits, and the tracking issue.
`reference-packages.toml` lists every package path in those pinned trees.

```bash
cargo xtask gaps            # ready work by priority, then blocked work
cargo xtask gaps show 119   # or a row id (B10) or issue key
cargo xtask gaps list --area brep --open
```

A landing that changes a row's level edits the row in the same commit. The
gate (`cargo xtask gaps check`) rejects evidence that no longer resolves,
reference paths outside the pinned trees, and dangling issue keys. It does
not judge whether a grade is true; review does.

Re-pinning OCCT or CGAL means regenerating `reference-packages.toml` from the
new tree and re-checking every row's reference paths. The OCCT/CGAL sources
are not needed for anything else: the ledger carries what the comparison
learned from them.
