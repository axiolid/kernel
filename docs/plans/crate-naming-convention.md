# Crate naming convention, enforced

Status: done (ADR 0064)

## The rule

Derived from `[package.metadata.axiolid]` `role` + `domain`, which the
architecture model already parses. Names are checked, not conventions
remembered.

1. `role = contract.operation` -> `axiolid-<domain>-contract`
2. `role = provider.*`         -> `axiolid-<domain>-<engine>`, engine non-empty
3. anything else               -> must NOT end in `-contract`, and must NOT
   equal a contract's name with `-contract` stripped

`<domain>` is the `domain` field with `.` replaced by `-`.

Rule 3's second half is the one that catches the real defect: the
unsuffixed capability name reads as THE implementation of that
capability. Rule 2 means providers always carry an engine suffix, so
the unsuffixed name is reserved and belongs to nobody.

## Why derive from metadata rather than a hand-written list

ADR 0004 chose machine-checkable metadata over naming convention. That
was right, and it is why naming drifted: nothing read the names. The
fix is not to abandon the metadata but to DERIVE the expected name from
it, so the two cannot disagree silently.

## Current state against the rule

Compliant, no change: `axiolid-mesh-boolean-contract`,
`axiolid-mesh-section-contract`, `axiolid-curve-evaluate-contract`,
`axiolid-mesh-compile-contract`, `axiolid-mesh-boolean-boolmesh`,
`axiolid-pointcloud-reconstruction-sdf`.

### Fixed here (metadata only, not breaking)

- `axiolid-pointcloud-reconstruction-contract` domain
  `operation.pointcloud-reconstruction` -> `pointcloud.reconstruction`.
  The `operation.` prefix restated the role; stripping it makes the
  derived name match the actual name, so the crate becomes compliant
  WITHOUT a rename. Its provider already used `pointcloud.reconstruction`,
  so this also makes contract and provider agree on their domain.
- `axiolid-tessellation-contract` domain `operation.tessellate` ->
  `tessellate`. Does not make it compliant, but reduces the violation to
  exactly one thing: noun vs verb in the crate name.

### Grandfathered (renames are BREAKING; all three are published)

Confirmed on the sparse index: 0.1.0 and 0.2.0 both live for each.

| Crate | Should be | Why deferred |
| --- | --- | --- |
| `axiolid-tessellation-contract` | `axiolid-tessellate-contract` | published; 26 refs |
| `axiolid-exact-compile-contract` | `axiolid-brep-compile-contract` | published; domain is `brep.compile` |
| `axiolid-mesh-compile` | `axiolid-graph-compile` | published; 48 refs; impersonates the contract |

Deferred to 0.3.0, the next breaking release. The exception list is
CLOSED: a new violation fails the gate. An exception that has been fixed
ALSO fails, so the list cannot rot into a permanent amnesty.

