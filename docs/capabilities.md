# Capabilities and status

Axiolid is **striving to be a multipurpose, exact B-rep geometry kernel** —
usable for CAD construction, for rule checking over building models, and for
the analytical work those applications need. It is explicitly *not* intended to
be another tessellation pipeline. Exact geometry is the model; tessellation is
one output you can ask for, with a tolerance you choose. See
[ADR 0020](/adr/0020-exact-brep-kernel-model).

That is the **ambition**. This page separates it from the **status**, which is
deliberately conservative:

- **Implemented** — a focused crate or provider has executable behavior and tests.
- **Represented** — the neutral type vocabulary or contract exists.
- **Seam** — an extension boundary exists; not a claim that Axiolid bundles a
  production implementation.
- **Planned** — accepted as in-scope by an ADR, not yet built.

Where the two differ, the status column wins. An accepted ambition is not
evidence that an algorithm exists.

::: tip Looking for a quick answer?
The [support map](/support-map) is the scannable version: one row per
capability, red/yellow/green, including what is deliberately **out** of scope
and which Cargo feature turns each thing on. This page is the detailed
prose account.
:::

## Where the implementation stands against the ambition

Graph execution now has **separate exact and discrete paths**. The reference exact
compiler memoises supported graph extrusion roots as `ExactBRep`, so a sharp
filled/hollow rectangle remains planar and a positive-axis filled circle remains
a cylinder through compilation. Its cache cannot contain `TriMesh`; unsupported
node, operation, profile, and oblique-circle families fail with typed
`UnsupportedInput` diagnostics. The reference mesh compiler remains an explicit,
tolerance-driven discrete capability with a `TriMesh` cache.

This is a focused exact slice, not general exact graph evaluation. Treat exactness
claims as applying only to the supported extrusion families and the independently
certified affine trace arrangement described below.

## Core and reference work

| Capability | Status | Evidence / boundary |
| --- | --- | --- |
| Scalar values, frames, transforms, bounds, intervals, tolerance | Implemented | `axiolid-core` |
| Mesh values, polygon/triangle utilities and views | Implemented | `axiolid-mesh` |
| Robust-orientation and in-circle / in-sphere predicate reference paths | Implemented | `axiolid-predicates` with degeneracy and filter tests; re-exported by `axiolid-reference`. Depends only on `axiolid-core` and `axiolid-guarantees`, so a consumer gets certified signs without curves, surfaces, or meshes ([ADR 0036](/adr/0036-use-case-specific-compilation-closures)) |
| Scalar geometry generation | Implemented broad discrete path plus focused exact construction | `axiolid-construct` builds exact sharp filled/hollow rectangle extrusions, positive-axis filled-circle extrusions, and the certified one-owner affine trace arrangement. Other profile/path families remain discrete or fail closed; see [ADR 0024](/adr/0024-exact-brep-result-contracts) and [ADR 0029](/adr/0029-certified-trace-topology-integration) |
| Scalar graph compilation | Implemented separate mesh and focused exact paths | `axiolid-mesh-compile`; `ReferenceExactCompiler` owns a per-batch `NodeId -> ExactBRep` cache and supports only the exact extrusion families above. `ReferenceMeshCompiler` remains the explicit tolerance-driven discrete path. Unsupported exact roots refuse instead of tessellating; see [ADR 0020](/adr/0020-exact-brep-kernel-model) |
| Polygon triangulation | Implemented provider | `axiolid-tessellation-contract` adopts Earcut under the contract in [ADR 0015](/adr/0015-adopt-earcut-polygon-triangulation) |
| Mesh Boolean execution | Optional provider | `axiolid-mesh-boolean-boolmesh`; limited to its mesh contract and tests |
| Mesh plane section | Implemented scalar oracle | `axiolid-mesh-section-contract::MeshPlaneSection` plus `axiolid-reference::ScalarSection`; exact binary64 plane-side classification, source-topology stitching of outward-oriented closed solids, closed plane-local contours, explicit limits/cancellation/scratch preflight, and input-mesh provenance. Coplanar surface overlap and open/non-manifold results fail closed; no analytic exactness or region nesting is claimed. See [ADR 0033](/adr/0033-mesh-plane-section-contract) |

## Geometry representation

| Capability | Status | Notes |
| --- | --- | --- |
| Primitive solids and half-spaces | Represented | `axiolid-primitive` owns neutral values and validation |
| Profiles / contours | Represented with validation | `axiolid-profile` |
| Curves and surfaces | Represented | `axiolid-curve`, `axiolid-surface`; representation alone is not an evaluator claim |
| Polynomial/rational B-spline evaluation | Implemented scalar oracle | `axiolid-evaluate` (re-exported by `axiolid-reference`); validates compact knot/control/weight data and exposes homogeneous point, first-, and second-derivative curve/surface jets |
| Linear intersection: line/line and segment/segment in 2D | Implemented | `axiolid-linear-intersection`; certified parallel/coincident/crossing classification, exact endpoint-contact parameters, collinear overlap spans, and typed refusals naming the operand at fault. 3D linear and curve/surface intersection are not implemented |
| Analytic curve evaluation and adaptive flattening | Implemented | `axiolid-evaluate` (re-exported by `axiolid-reference`); native-parameter evaluation with a chord-error budget, fail-closed on depth exhaustion and unbisectable intervals |
| Point→parameter inversion, elementary surfaces | Implemented | `axiolid-evaluate` (re-exported by `axiolid-reference`); plane/cylinder/cone/sphere/torus with residual validation, refusing degenerate axis/apex/pole points |
| General NURBS differential analysis | Implemented reference path | `axiolid-nurbs`; curve tangents/curvature and surface fundamental forms, normals, Gaussian/mean/principal curvature; see [ADR 0022](/adr/0022-general-nurbs-kernel-capability) |
| NURBS inverse queries | Implemented reference paths plus bounded certified slices | `axiolid-nurbs`; deterministic budgeted multistart curve/surface candidates, globally certified clamped curve projection and curve-pair distance, and globally certified open clamped and explicit cyclic-periodic surface projection over the full rectangular or quotient domain. The surface certificate retains every possible global-minimizer box, requires distance and parameter resolution, and refuses trims or neutral closed axes lacking the explicit cyclic schema. Closest-point INVERSION builds on that certificate and additionally proves the answer is unique -- the minimizer cover must be one connected, localized region -- so poles, seams, and self-touching patches yield an explicit ambiguity refusal carrying the rival boxes rather than an arbitrary representative; see [ADR 0025](/adr/0025-certified-nurbs-subdivision-oracle) and [ADR 0030](/adr/0030-globally-certified-surface-projection) |
| Planar clamped NURBS curve/curve roots | Implemented bounded reference slice | `axiolid-nurbs`; exact-sign single-span line/point predicates with a distinct zero-length `PointContact`, contractive transverse polynomial/rational Bézier boxes via strict-interior Krawczyk proof, explicit native-parameter resolution, localized structural overlap/endpoint tangency, compact parameter-only DFS work items, hard work ceilings, per-box contact classification with boundary ownership and deduplication, interval-certified tangency proof over the whole box, and an explicit `BoundaryCrossing` class for transverse roots on a shared cell edge that strict-interior isolation cannot name. General seam/trim ownership and curved overlap tracing remain unimplemented; see [ADR 0026](/adr/0026-certified-planar-nurbs-root-isolation) |
| Clamped, internally continuous NURBS curve/surface roots | Implemented bounded reference slice | `axiolid-nurbs`; accepts internal knot multiplicity `1..=degree` and treats valid full-multiplicity internal knots as unsupported by this certified query; outward tensor rational-Bézier refinement, conservative native-span curve/surface derivative hulls, strict-interior 3×3 Krawczyk existence/uniqueness proofs, explicit three-parameter resolution, compact parameter-only DFS work items, shared hard work ceilings, and explicit unresolved singular/tangential/boundary boxes. General contact classification, periodic/seam/trim ownership, discontinuous span joins, surface/surface tracing, and topology remain unimplemented; see [ADR 0027](/adr/0027-certified-nurbs-curve-surface-root-isolation) |
| NURBS exact transformations | Implemented | `axiolid-nurbs`; homogeneous curve knot insertion/reversal/split/Bézier decomposition and surface U/V insertion/reversal preserve the represented shape |
| Verified periodic curve views | Implemented bounded semantics | `axiolid-nurbs`; `PeriodicCurve2`/`PeriodicCurve3` require declared closure plus verified positional seam continuity, report C0/C1/C2 continuity, wrap only explicit view evaluation, canonicalize interior insertion/split parameters, revalidate insertion, and return neutral edited/open curves. They do not define periodic control-net or surface topology; see [ADR 0031](/adr/0031-verified-periodic-curve-views) |
| Explicit periodic B-spline surfaces | Implemented fixed-topology semantics | Cyclic U/V/UV schema, wrapped jets, alias-safe edits, and certified projection; see [ADR 0032](/adr/0032-explicit-periodic-bspline-surfaces) |
| Trimmed curved-face tessellation | Implemented bounded reference path | `axiolid-mesh-compile`; endpoint-inclusive pcurve boundaries and holes, guarded structured-grid/Earcut seeds, topological seam reuse, explicit outer/bound orientation, face-level conforming support-surface refinement, periodic face charts, and fail-closed input/per-edge/per-face/aggregate/depth limits; see [ADR 0019](/adr/0019-validate-and-refine-nurbs-on-the-scalar-read-path) |
| Topology / B-rep vocabulary | Represented | `axiolid-topology` provides generic typed-role topology; `axiolid-brep` provides strict owned analytic catalogs and required native trim spans |
| Ray/triangle-mesh narrow phase | Implemented | `axiolid-ray-mesh`; nearest hit as parametric distance in direction units, triangle index, barycentric coordinates, and a certified `orient3d` front/back/coplanar side. Double-sided by default with deterministic lowest-index tie-breaking, typed refusals for zero-area triangles and out-of-range indices, and `nearest_hit_among` for composing with the `axiolid-spatial` BVH broad phase. It does not decide what a ray means: sampling patterns, entity identity, and obstruction policy stay with the caller |
| Independent mapped-3D verification of intersection/inversion | Implemented verification tool | `axiolid-oracle` (workspace-internal, not published in the facade); maps claimed parameter boxes back through `axiolid-evaluate` and measures the model-space deviation, sharing no subdivision or interval machinery with `axiolid-nurbs`. It refutes overstated global minimum distances soundly and reports the measured 3D deviation; a small sampled deviation is a witness, never a proof of global minimality; see [ADR 0037](/adr/0037-mapped-3d-verification-oracle) |
| Exact B-rep result contract | Implemented contract with focused constructors | `axiolid-brep::ExactBRep` refuses missing supports/spans and generic topology failures. `axiolid-construct` builds supported exact extrusions and the certified one-owner affine trace arrangement; see [ADR 0024](/adr/0024-exact-brep-result-contracts) and [ADR 0029](/adr/0029-certified-trace-topology-integration) |
| Exact B-rep operation results | Focused extrusion, revolution, boolean, chamfer and affine-arrangement slices; general operations planned | Sharp filled/hollow rectangles and positive-axis filled circles retain analytic planes/cylinders, closed trims, pcurves, and native spans through exact graph compilation. Full-turn revolution of a sharp filled rectangle about an offset axis parallel to the profile's local y yields an exact annular tube. Boolean over two COAXIAL prisms is exact via the planar overlay, including a wall with an interior opening; stepped results (unions with differing spans, tools shorter than the subject) and disconnected results refuse by name. Constant-distance chamfer on one straight vertical edge of an extruded sharp rectangle is exact. One bounded surface-pair slice also returns analytic supports and a residual certificate. Partial-turn revolutions, profiles crossing the revolution axis, fillets, non-coaxial or non-prism boolean operands, other profiles, general sweeps, instances, authored B-rep graph nodes, and curved intersections remain unsupported |
| Single-span affine NURBS surface/surface traces | Implemented bounded reference slice with focused topology handoff | `axiolid-nurbs` proves affine transverse bounded traces. `axiolid-construct` splits the uniquely boundary-owned rectangular patch and attaches the same edge as an embedded pcurve on the containing face. Dual-boundary ownership, corners, mixed ownership, curved tracing, loops, tangency/coincidence, and multispan stitching remain unresolved. Constructed intersection CURVES are emitted only where the construction is exact -- a degree-1 segment for the affine family, with a deviation bound valid over the whole curve -- and curve/surface queries return isolated crossing points rather than a manufactured extent; every other case refuses by name and proven `Disjoint` is distinguished from `Unresolved`. See [ADR 0038](/adr/0038-constructed-intersection-curves) and [ADR 0028](/adr/0028-certified-affine-surface-surface-tracing) and [ADR 0029](/adr/0029-certified-trace-topology-integration) |
| Immutable shared geometry DAG | Implemented structural model | `axiolid-model` uses typed IDs and backward references |
| Sweeps / extrusions / revolutions / lofts | Broad discrete reference path plus focused exact extrusion, revolution and straight sweep | `axiolid-construct` exactly constructs sharp filled/hollow rectangle and positive-axis filled-circle extrusions, full-turn revolutions of a sharp filled rectangle about an offset parallel axis, and fixed-reference sweeps along a STRAIGHT directrix (which delegate to the extrusion they are). Curved directrices, partial turns, other extrusion profiles, general sweeps, and lofts remain discrete or refuse exact requests. The L2 result boundary keeps exact B-rep and tolerance-bearing tessellation requests distinct and never substitutes a mesh for exact output |
| Mass properties | Implemented for meshes and planar exact B-reps | `axiolid-measure`: `MeshMeasure` implements `Measure<TriMesh>` (area, signed volume, volume centroid, second moments), and `exact_properties` measures an `ExactBRep` WITHOUT tessellating it, so a planar-faced prism measures at machine precision rather than at tessellation fidelity. Verified against closed forms, the parallel-axis theorem, and an exact-vs-mesh differential. Curved exact faces are refused by name; open shells are refused rather than assigned a plausible volume. Exact support is behind the `exact` feature so a mesh-only consumer does not acquire B-rep geometry |
| Mesh defect diagnosis | Implemented: located defects, not counts | Coplanar overlap decided exactly via `orient2d`; repair is separate |
| Convex hull (3D) | Implemented for point sets | `axiolid-construct`: `convex_hull` builds a closed, outward-oriented `TriMesh` by incremental insertion. Every visibility decision goes through certified `orient3d`; there is no epsilon in the module. Degenerate input is typed rather than approximated: fewer than four points, all-collinear, and all-coplanar are distinguishable refusals, so a caller can fall back to the 2D hull in `axiolid-reference`. Verified against closed-form volumes and by re-deciding containment for every input point against every face |
| Exact boolean (planar-faced solids) | Implemented for general planar-faced operands | `axiolid-construct`: `boolean_polyhedra_exact` accepts arbitrary planar-faced solids, convex or non-convex, in any orientation, for union, intersection and difference. Faces are split only against ORIGINAL input planes, never derived ones, so vertices stay one intersection from input data. Containment uses exact ray-crossing parity, which a convex-only "inside every plane" test gets wrong for concave regions. Coplanar contact is classified by normal agreement so a shared face is never duplicated. Verified by the volume identity |A u B| + |A n B| = |A| + |B|, a boolmesh differential, and v0.7 diagnosis. Curved operands (including v0.6 revolution output) are REFUSED by name, not tessellated |
| Solid offset and shelling | Implemented for planar-faced solids | `axiolid-construct`: `offset_solid` (miter, outward and inward) and `shell_solid`. Vertices are moved to the meeting point of their incident pushed planes, so concave edges are handled and topology is preserved. An offset that closes the gap between opposing walls is REFUSED, not emitted. NOT sphere offset, NOT variable-distance, NOT Minkowski, and no face removal to open a shell |
| Constant-radius fillet | Implemented for one straight edge of a prism | `axiolid-construct`: `fillet_extruded_profile` emits a real `Surface::Cylinder` blend tangent to both adjacent walls, not a many-segment chamfer. Tangency follows from placing the axis on the internal bisector and is verified by measuring axis-to-wall distance. Refuses radii reaching past a neighbouring corner, variable radius, edge loops, and curved edges |
| Mesh decimation | Implemented, deviation-bounded | `axiolid-decimate`: edge-collapse reduction to a triangle budget or a maximum deviation. The deviation is measured and reported per run, not estimated, and accumulates per vertex so repeated collapses cannot drift past the bound. Collapses that would create a non-manifold edge are refused and counted. Output is deterministic. NOT quadric-based: cost is edge length and placement is the midpoint, which is predictable but not optimal for a given budget; see the crate PLAN.md, which also records that the normal-inversion guard lacks a mutation probe |
| Connected-component decomposition | Implemented: `decompose` / `compose` | Connectivity is by shared vertex index, not geometric proximity; coincident-but-distinct vertices stay separate bodies until welded |
| Clearance / min gap | Implemented for triangle meshes | BVH-accelerated; reports clash as 0.0 and out-of-range as none |
| Point containment / winding number | Implemented, exact | Certified `orient3d` parity; signed, so inside-out shells are visible |
| Ray casting | Implemented for triangle meshes | Exact membership; only the hit parameter is f64 |
| Genus | Implemented for closed two-manifolds | Refuses meshes with boundary or non-manifold edges |
| Boolean construction arithmetic | f64 constructions; predicates verify only | The production mesh boolean constructs intersection coordinates in `f64`; `axiolid-predicates` is a DEV-dependency of the provider and re-decides orientation in tests rather than driving construction. Measured against Manifold, CGAL (exact constructions) and OCCT on byte-identical operands: error does not accumulate over chained operations (5.62e-16 at n=64, at or better than the exact kernels), and degradation appears only as operands approach coincidence. Published thresholds by relative operand overlap: above 1e-6 invisible, 1e-6 to 1e-12 smooth degradation, below 1e-12 severe. Axiolid never collapses to zero volume on this sweep, where Manifold and OCCT both do. `BooleanEvidence::relative_overlap` reports the measured conditioning so a caller can refuse on its own threshold; Axiolid does not choose one on its behalf; see [ADR 0045](/adr/0045-boolean-construction-arithmetic) |
| Spatial, healing | Focused crates / staged capability | Consult each crate’s `PLAN.md`; do not infer broad CAD coverage |

## Execution and acceleration

| Capability | Status | Boundary |
| --- | --- | --- |
| Portable CPU context | Implemented shell | `axiolid-backend-cpu`; portable defaults and explicit feature tiers |
| Parallel / SIMD | Opt-in context features | They require measurement and differential validation before performance claims |
| GPU graph execution | Contract seam | `axiolid-backend-gpu` provides an API-neutral seam, not a bundled GPU algorithm suite |
| Native CUDA/HIP | Planned out-of-tree providers | See [ADR 0011](/adr/0011-native-accelerator-backends-out-of-tree) |

Performance work is deliberately deferred; see
[ADR 0013](/adr/0013-deferred-performance-techniques) and the
[roadmap](./ROADMAP.md). Correctness and exactness come first, and the
architecture is kept clean so optimization stays possible later.

## In scope, not yet built

Accepted by [ADR 0020](/adr/0020-exact-brep-kernel-model) as work the kernel
intends to do. These broader exact/certified forms do not exist today:

- General exact B-rep results for unsupported graph/operation families (including
  booleans, revolutions, lofts, instances, authored B-rep roots, and arbitrary
  profile sweeps), plus exact-result propagation beyond the focused extrusion
  and affine-arrangement slices.
- General curved surface/surface intersection-curve tracing, closed loops, and
  tangent/overlap or seam/trim-owned completion beyond the bounded affine trace slice.
- Exact analytic booleans, section curves, offsets, and fillets, which all sit
  downstream of intersection.

## Explicit non-goals today

- Source-format parsing or semantic interpretation.
- A claim of OpenCascade compatibility or replacement coverage.
- Bundled production CUDA, HIP, Metal, Vulkan, or WebGPU compute kernels.
- A global hidden tolerance policy.


## How to evaluate a claim

Use the evidence nearest the implementation:

1. The relevant crate’s public API and tests.
2. The architecture decision that defines the contract.
3. Feature-isolation and layering gates.
4. Benchmark reports for performance statements.

The [research comparison](./research/geometry-kernel-capability-comparison.md) is useful context, but it is not a capability declaration.
