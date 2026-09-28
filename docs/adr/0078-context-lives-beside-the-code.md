# 0078 — Context lives beside the code; open work lives in issues

- **Status:** Accepted
- **Date:** 2026-09-28
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

The repository kept its contributor context in a tree of progressive
`AGENTS.md` files, one per ownership directory and crate (about ninety),
plus a `PLAN.md` beside many crates, `docs/plans/`, and root-level agent
session plans (`PLAN-*.md`, `*-FINDING.md`). Each was meant to add only
what its directory needed on top of its parent.

In practice they drifted. Most crate files repeated the same boilerplate
("follow the parent", "public values derive `Debug` and `Clone`",
dependency lists that `cargo xtask architecture check` already enforces).
The statements that were unique were often stale: kernel#25 found six of
seven unchecked `PLAN.md` boxes describing work that was already done,
and `scripts/check-roadmap-freshness.py` grew a special case to stop
`PLAN.md` from recording status. A reader could not tell which of the
three or four files on a path was authoritative, and crates.io showed
the workspace README for every crate.

openbimrs/ifc made the same move in its PR #179 and has run on it since.

## Decision

Every piece of context has exactly one home, chosen by what it is:

| Content | Home |
| --- | --- |
| Rules every contributor must know before touching the repository | the root `AGENTS.md`, the only file of that name (at most 120 lines) |
| A durable decision and its reasoning | an ADR in `docs/adr/` |
| An invariant, pitfall or reason tied to specific code | the `//!` or `///` doc of that code |
| A rule that can be checked | a test, `cargo xtask architecture check` or `cargo xtask context check` |
| Open work, gaps, next steps | a GitHub issue; a code marker is written `TODO(#N)` |
| What a crate is for and what it deliberately does not do, when nothing above holds it | the crate's `README.md` (its crates.io page, at most 150 lines) |

No nested `AGENTS.md` or `CLAUDE.md`, no `PLAN.md`, no `docs/plans/` and
no checked-in session plans. The repository root holds only the workspace
manifest, lockfile, toolchain pin, licence, `README.md` and `AGENTS.md`.
The machine-checked declarations that used to sit in a root
`architecture/` directory (closure profiles, capability ledger, reference
packages, semver exceptions) move to `docs/architecture/`, beside the maps
generated from them. A crate README repeats nothing a test, ADR
or module doc already says. It has no checkboxes and makes no claim
that goes stale at the next release, such as "not yet published".

`cargo xtask context check` enforces this and runs in `scripts/gate.sh`.
`scripts/probe_context_gate.sh` shows each rule can fail.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Keep nested `AGENTS.md`, lint for drift | The drift is semantic (a true-looking statement about code that changed), which no lint catches. Module docs sit next to the code and get reviewed with it. |
| Keep `PLAN.md` without status (the kernel#25 rule) | A plan without status still goes stale. It duplicates the issue tracker, which already has status, discussion and cross-links. |
| One `AGENTS.md` per crate, no parents | This still duplicates the README and module docs, and agents that read the nearest file still miss the root rules. |

## Consequences

**Positive**

- One place to look for each kind of fact. A reviewer sees a stale
  invariant in the diff that changes the code it describes.
- crates.io shows each crate's own page.
- Open work is visible, triaged and linked on the project board.

**Negative**

- Agents no longer get directory-specific context just by opening a
  directory. They must read the crate README and the module docs they
  touch, and the root `AGENTS.md` says so.
- Branches that edit a nested `AGENTS.md` or `PLAN.md` conflict with the
  deletion. Move the edit into the new home when rebasing.
