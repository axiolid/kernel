# Roadmap

::: warning This page does not track status
Every work item lives on GitHub. This page explains **why the work is ordered
the way it is** — the part that does not change when a checkbox gets ticked.

- **What is planned, in progress, or done** → [project board](https://github.com/orgs/axiolid/projects/1)
- **What a version must satisfy** → [milestones](https://github.com/axiolid/kernel/milestones)
- **What exists today** → [capabilities](./capabilities.md)
- **How to propose something** → [where things go](./contributing/where-things-go.md)

If this page and GitHub disagree, GitHub wins. Deliberately, there is not a
single item checklist below — a checklist here would go stale the first time
someone closed an issue.
:::

Axiolid is an early geometry kernel **striving to become a multipurpose, exact
B-rep kernel** — for CAD construction, for rule checking over building models,
and for the analysis those applications need. This page orders that work; it is
not a promise of dates, and listing something here is not a claim that it
exists.

## The ordering principle

**Build the meat of the kernel first, shape it for performance later.**

Correctness and exact geometry come before speed. We keep following good
practice — layering, feature isolation, provider seams, measured oracles — so
that optimization stays possible when its turn comes. We do not trade
capability for benchmarks now.

This is why performance work is *parked* rather than *declined*: the ordering
is intentional, and each parked item carries a written unblock condition. See
[ADR 0013](./adr/0013-deferred-performance-techniques.md).

## Why the milestones are shaped this way

Milestones are **capability gates, not dates**. Each one exists because the
next cannot be honestly attempted without it.

| Milestone | Why it comes when it does |
|---|---|
| [Exact B-rep survives operations](https://github.com/axiolid/kernel/milestone/1) | The first exact graph path must preserve analytic identity instead of routing through the discrete cache: a cylinder that enters a supported extrusion operation leaves as a cylinder. The focused exact compiler/cache establishes that invariant; later milestones expand the supported operation families ([ADR 0020](./adr/0020-exact-brep-kernel-model.md)). |
| [Intersection and inversion](https://github.com/axiolid/kernel/milestone/2) | Exact booleans, section curves, offsets, and fillets all reduce to intersection and inversion. Attempting them before the exact representation holds would build on sand. |
| [Downstream integration](https://github.com/axiolid/kernel/milestone/8) | A kernel nobody can consume is a research artefact. Versioned Rust and native integration surfaces, reproducible builds, and executable consumer probes come once there is a real capability worth integrating against. |
| [Trustworthy discrete geometry](https://github.com/axiolid/kernel/milestone/3) | The mesh path stays supported for callers who explicitly want discrete results — and it is the differential oracle for the exact path, so it must be trustworthy *after* there is an exact path to check against. |
| [Compiled geometry and plans](https://github.com/axiolid/kernel/milestone/4) | Reproducible operation plans need stable exact semantics underneath. Freezing a plan format over shifting geometry would bake in the wrong contract. |
| [Measurement and validity](https://github.com/axiolid/kernel/milestone/9) | Measurement is the oracle every later capability is checked against: a simplification that preserved shape, a repair that fixed a defect, and an offset that honoured a thickness are all unverifiable without volume and area. Defect diagnosis lands with it, because a repair that cannot name what it changed is not auditable. |
| [Mesh processing breadth](https://github.com/axiolid/kernel/milestone/10) | Repair, decimation, and hull. Each is only trustworthy once measurement and validity can prove what it preserved, so this cannot come first — a decimator with no error bound is damage, not simplification. |
| [Solid modelling breadth](https://github.com/axiolid/kernel/milestone/11) | General exact boolean, offset/shelling, and fillet. The heaviest geometry in the sequence and the most dependent on the validity foundation beneath it: these operations produce the self-intersections that measurement and validity detects and mesh processing repairs. |
| [Mesh-side breadth without exact loss](https://github.com/axiolid/kernel/milestone/12) | Closes the functional gaps against mesh-first kernels: per-vertex attributes that survive a boolean, component decomposition, refinement and smoothing, implicit level-set input, and the standard mesh queries. Scheduled after solid modelling breadth because these are additive breadth, not foundations -- and constrained by it: a mesh capability may not be bought by degrading an exact one. Refinement is the clearest case for doing this here rather than adopting a mesh library, since we keep the analytic surface and can evaluate it instead of interpolating between triangles. |
| [v1.0: Stable public API](https://github.com/axiolid/kernel/milestone/5) | An API is only worth stabilising once the capability surface behind it is real. Committing earlier would freeze scaffolding. |
| [Backlog: accepted, unscheduled](https://github.com/axiolid/kernel/milestone/6) | Not a capability gate and deliberately undated: work that is accepted and wanted but not committed to a release. Pull from here when a gate above is clear, rather than widening a gate in flight. |

Minor versions can be inserted between any two of these at any time. The
sequence is an ordering constraint, not a fixed count.

## What this kernel refuses to do

These are not backlog items. They are boundaries, and they do not expire.

- **Not a file parser.** IFC, STEP, and CAD source semantics belong in adapter
  projects. Representing imported data is not the same as parsing a format.
- **No C++ in the dependency graph.** Including OpenCascade. Pure Rust is a
  hard constraint, not a preference.
- **A type is not an implementation.** A feature flag, provider seam, or
  vocabulary type is never evidence that an algorithm is production-ready.
  Capability claims need an executable implementation and a test.
- **No performance claim without same-harness evidence.** Benchmarks currently
  exist for `axiolid-reference` and `axiolid-mesh-boolean-boolmesh` only. There
  are no broad performance claims to defend, and none should be made.

## Release discipline

Semantic versioning: patch for compatible fixes, minor for compatible
capability additions during `0.x`, and a documented breaking-change policy
before `1.0`. The [changelog](./CHANGELOG.md) records user-visible changes.

Releases are deliberate, reviewed changes rather than automatic tags. Automating
that is itself tracked work, not a claim about today.
