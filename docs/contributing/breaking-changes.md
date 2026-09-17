# Breaking-change policy

A policy nobody can check is not in force. This page states what may change
at each version step, defines exactly what counts as public surface, and
names the gate that enfor...[truncated]
## What counts as public surface

Every crate published to crates.io. The facade re-exports leaf crates, so
an item is public if a consumer can name it through ANY published path,
not only through `axiolid`. Removing an item from a leaf crate breaks the
facade re-export too.

Not public surface:

- crates with `publish = false` (`xtask`, `axiolid-benchmark`, `axiolid-oracle`)
- `#[doc(hidden)]` items
- private modules, and anything reachable only through them
- `Cargo.toml` feature names are public: a consumer writes them down

## What may change at each step

The workspace shares one version, so the rules apply to the workspace as
a whole. Pre-1.0, Cargo treats the MINOR field as the major: `0.2.x` and
`0.3.0` are incompatible, which is why a breaking change today lands as
a minor bump and not a patch.

| Change | Patch `0.2.0 -> 0.2.1` | Minor `0.2.x -> 0.3.0` | Major `0.x -> 1.0` |
| --- | --- | --- | --- |
| Add an item, variant behind `#[non_exhaustive]`, or feature | yes | yes | yes |
| Fix behaviour without changing a signature | yes | yes | yes |
| Remove or rename a public item | no | yes | yes |
| Change a signature, trait bound, or public field type | no | yes | yes |
| Add a required trait method | no | yes | yes |
| Add a variant to an enum NOT marked `#[non_exhaustive]` | no | yes | yes |
| Change an array length in a public constant | no | yes | yes |
| Tighten a refusal so previously-accepted input now errors | no | yes | yes |

The array-length row is not hypothetical. `capability_ids::ALL` was
`[CapabilityId; 9]`, and adding a tenth capability changed the type to
`[CapabilityId; 10]` -- a breaking change caused by registering a value.
It is now `&[CapabilityId]`, so the vocabulary grows additively.
A public collection whose LENGTH is part of its type makes every
addition breaking; prefer a slice.

## The gate

`scripts/check-semver.py`, wired into `scripts/gate.sh`, runs
[`cargo-semver-checks`](https://github.com/obi1kenobi/cargo-semver-checks)
against the last published release of every publishable crate. It
compares the working tree to what is actually on crates.io, so it
cannot be fooled by a stale baseline committed in-tree.

The rule it enforces:

> A breaking change is allowed only when the workspace version has
> already been bumped past the published one in a way that admits it.

So a breaking change with no bump fails; the same change with the
minor bump applied passes.

Run it directly:

```bash
python3 scripts/check-semver.py          # gate mode
python3 scripts/check-semver.py --explain  # show each crate's baseline
```

A crate with no published release yet is skipped: there is no baseline
to break. It joins the gate automatically on its first publish.

## When a breaking change is the right answer

This policy does not forbid breaking changes; it forbids SILENT ones.
Fixing a contract that cannot express a correct answer is worth a
bump. What is not acceptable is a consumer discovering the break at
compile time after a patch upgrade.

Record the change in `docs/CHANGELOG.md` under a `Changed` or
`Removed` heading, naming the replacement. A removal with no stated
migration is an unfinished change.
