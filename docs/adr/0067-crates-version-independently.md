# Crates version independently

Status: accepted

## Context

All 55 publishable crates shared one version through `version.workspace = true`.
The workspace exists in many small crates deliberately, so a consumer can take
`axiolid-core` without dragging in NURBS evaluation or mesh booleans — pay for
what you use.

Shared versioning defeats that at release time. Any change, however local,
moved every crate's version, so every crate had to be published again. The
publish script verifies each crate by packaging it twice — once to build the
archive, once to prove `cargo publish` reproduces identical bytes — and each
packaging compiles the crate. Across the whole workspace that is a multi-hour
release, most of it spent re-verifying crates whose source did not change.

The practical effect was that small fixes waited for a release window instead
of reaching consumers when they were ready.

## Decision

Each crate carries its own literal version and moves independently.

- `version.workspace = true` is replaced by an explicit `version` per manifest.
  All crates start at the version they already had, so adopting this changed
  nothing; versions diverge from the next change onward.
- `[workspace.dependencies]` already pinned an explicit `version` next to each
  `path`, so dependents already state a requirement. Nothing there changed.
- A patch or minor fix publishes only the crates that actually changed.
- Breaking changes still move the whole graph together, on purpose. Two majors
  of a foundation crate can coexist in one dependency tree, and the resulting
  type mismatch is confusing to diagnose, so a major is a workspace-wide event.
- The publish script checks the registry before doing any packaging work.
  Registry versions are immutable, so an already-published version has nothing
  to upload and nothing a byte comparison could change.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Keep shared versions, parallelise publishing | crates.io index propagation is the serialising constraint, and dependents must observe a dependency before they publish. Parallelism cannot remove that ordering, and it would still republish unchanged crates. |
| Keep shared versions, skip unchanged crates | The version is what identifies a release; if every crate's version moved, every crate genuinely does need publishing. The skip only becomes correct once versions can stay still. |
| Independent versions with no path deps | Local development would resolve against the registry, so a change could not be tested across crates before release. |
| Release-please or cargo-smart-release | Worth revisiting for changelog assembly, but both want to own version selection and commit conventions; adopting one is a larger change than the bottleneck requires. |

## Consequences

**Positive**

- A one-crate fix publishes one crate.
- A patch to a foundation crate reaches dependents through their existing
  caret requirement without republishing them.
- Version numbers start carrying information: a crate at `0.3.4` has changed
  four times since `0.3.0`, where before it only tracked the workspace.

**Negative / costs**

- A breaking foundation bump means editing every dependent's requirement, which
  shared versioning did implicitly.
- Reviewers can no longer assume two crates in the tree share a version.
- Release notes must be assembled from several sources rather than read from
  one file.

**Follow-ups / risks to watch**

- Nothing yet detects that a crate's source changed but its version did not.
  Until that exists, choosing which crates to bump is a manual judgement, and a
  missed bump ships a fix nobody can depend on.

## Relation to existing code

- `crates/**/Cargo.toml` — literal `version` per manifest.
- `crates/**/CHANGELOG.md` — one per publishable crate, Keep-a-Changelog.
- `Cargo.toml` — `[workspace.package] version` now seeds new crates only;
  `[workspace.dependencies]` still carries the requirement dependents use.
- `scripts/publish-workspace.py` — `version_already_published` gates the
  expensive verification.
- `scripts/prepare-crate-release.py` — bumps one crate's version and rolls its
  own changelog; refuses a bump its own caret compatibility would break.
- `scripts/assemble-crate-changelogs.py` — renders `docs/reference/changelog.md`
  from every crate's dated releases; `--check` gates drift in `scripts/gate.sh`.
