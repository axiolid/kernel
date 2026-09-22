# Ownership replaces read-only/mutable type pairs

Status: accepted

## Context

A capability audit compared this kernel against an interface catalogue from a
managed-language geometry API, which pairs every value type with a mutable
twin:

```
Vector3d    / MVector3d
Triangle2d  / MTriangle2d
AABB3d      / MAABB3d
Polygon2d   / MPolygon2d
```

Sixteen such pairs appear across its 2D and 3D primitives. The audit asked
whether the same pairs should exist here, and the primitives added in
`c8dc050` and `ae894d8` made the question concrete: `Triangle2`, `Triangle3`,
`Aabb2`, `Polygon2`, `Polygon3`, `Rectangle2`, `Rectangle3`, and `Box3` all
have an obvious `M*` counterpart in that catalogue and none of them here.

The pattern exists because Java and C# cannot express "you may read this but
not write it" in a type. A method returning `Triangle2d` is returning a
reference to a mutable object; the read-only interface is the only way to
withhold the setters. The split is a workaround for a missing language
feature, not a modelling decision.

Rust has that feature. `&Triangle2` withholds mutation, `&mut Triangle2`
grants it, and the compiler enforces the distinction at every call site. The
aliasing rules go further than the interface pair can: a `&` borrow guarantees
*nobody* is mutating the value, which a read-only interface cannot promise
because another holder of the concrete type may still write through it.

## Decision

We will not mirror value types into read-only and mutable variants. Each
primitive is one type, and access is governed by `&` versus `&mut`.

Where an invariant must survive mutation, we will express it with a private
field and a constructor or accessor pair on the single type, not with a second
type. Where a type is a plain aggregate with no invariant — which is every
primitive in `axiolid_core::primitives2` and `primitives3` — fields stay
public and the borrow checker is the whole mechanism.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Mirror the catalogue: add `MTriangle2` and friends | Doubles the type count and every conversion between them, to re-implement a guarantee the compiler already makes. The mutable twin would be structurally identical, so the pair carries no information. |
| Read-only *traits* (`Triangle2View`) implemented by the concrete types | Buys generic abstraction over storage, which no caller has asked for. Costs a trait bound on every signature and blocks field access. Revisit only if a second storage layout appears. |
| Immutable types plus `with_*` builders returning copies | Reasonable for small `Copy` types, but `Polygon2`/`Polygon3` own a `Vec`, so every edit would clone the boundary. Silent O(n) cost on what looks like a setter. |

## Consequences

**Positive**

- One type per concept, so there is no conversion surface between a value and
  its mutable twin, and no question about which one an API should accept.
- The guarantee is stronger than the interface pair provides: `&` proves no
  aliased writer exists, where a read-only interface only proves *this* holder
  will not write.
- Signatures document intent without extra vocabulary: `&Triangle3` reads as
  input, `&mut Triangle3` as in-place edit.

**Negative / costs**

- Consumers porting from a catalogue-shaped API will not find the names they
  expect and must learn the borrow forms instead. This ADR is the answer to
  that question when it recurs.
- Language bindings that cross the FFI boundary (C ABI, future Python or
  JavaScript packages) must re-introduce the distinction themselves, because
  the host languages lack it. That is a binding-layer concern, deliberately
  not pushed into the kernel's type system.

**Follow-ups / risks to watch**

- If a primitive later gains an invariant that mutation could break — a
  polygon that must stay simple, say — the fix is a private field plus
  validated accessors on that type, not a second type. Watch for the mirror
  pattern reappearing under a different name.

## Relation to existing code

- `crates/foundation/core/src/primitives2.rs` — `Aabb2`, `Rectangle2`,
  `Triangle2`, `Polygon2`, all single types with public fields.
- `crates/foundation/core/src/primitives3.rs` — `Rectangle3`, `Box3`,
  `Polygon3`, `Triangle3`, likewise.
- `crates/foundation/core/src/bounds.rs` — `Aabb` predates the audit and
  already follows this shape: `extend` and `union` take `&mut self` rather
  than living on a separate mutable type.
