# Support map

What the kernel can do today, what it partly does, and what it does not do yet.

Every row is one capability. The legend is deliberately blunt:

| | Meaning |
| --- | --- |
| 🟢 | **Implemented.** Executable behaviour, tests, and at least one consumer. You can use it now. |
| 🟡 | **Partial.** Something real exists, but it is narrower than the name suggests. The note says exactly how. |
| 🔴 | **Not implemented.** Either deliberately out of scope, or plausible future work. The note says which. |

A 🔴 row is not a promise. It is here so you can tell "we decided against this"
from "nobody has built it yet" without reading the issue tracker.

For the longer, prose-form status — including the exact/discrete split and what
"exact" is claimed to mean — see [Capabilities and status](/capabilities). Where
the two disagree, that page wins.

## How to read the optional column

Axiolid is a kernel, not an application. Almost everything below is behind an
**additive Cargo feature**: you compile the geometry you use and nothing else.
The *Optional* column names the feature, or `core` when it is always present.

```toml
# Points only. Resolves to axiolid + axiolid-core + axiolid-pointcloud.
axiolid = { version = "0.11", default-features = false, features = ["pointcloud"] }
```

---

## Spatial structures

Indices that answer *which candidates*, never *what the intersection is*.

| Capability | Support | Optional | Notes |
| --- | :---: | --- | --- |
| Bounding volume hierarchy (BVH) | 🟢 | `spatial` | `Bvh<K>` over anything with an AABB. Deterministic median split; AABB, ray, overlap-pair, and filtered-nearest queries. Used by clash, ray casting, and healing. |
| Uniform point grid (KNN / radius) | 🟢 | `pointcloud-queries` | `PointIndex`. Exact distances, not broad-phase bounds, with deterministic tie-breaking by point index. Verified against brute force. |
| Edge adjacency over a mesh | 🟢 | `core` | `EdgeAdjacency`: boundary, non-manifold, inconsistent-winding, vertex/triangle neighbours, Euler characteristic. Derived once and shared, so algorithms stop rebuilding it. |
| Octree | 🔴 | — | Not implemented, and not currently justified: the BVH covers object queries and the grid covers point queries. An octree's advantage is sparse volumetric subdivision. It would arrive with a measured workload, not before. |
| k-d tree | 🔴 | — | Not implemented. Would overlap `PointIndex` for the queries we actually make. |
| R-tree / quadtree / BSP | 🔴 | — | Not implemented, no consumer. |
| Morton / Z-order codes | 🔴 | — | Not computed by us. The `boolmesh` provider builds its own internally. |
| Half-edge with mutable connectivity | 🔴 | — | `EdgeAdjacency` is an immutable derived index; every consumer reads adjacency and writes a *new* mesh. In-place connectivity editing would be a separate type, not a field on this one. |

## Geometry representation

| Capability | Support | Optional | Notes |
| --- | :---: | --- | --- |
| Scalars, frames, transforms, bounds, tolerance | 🟢 | `core` | `axiolid-core`. The vocabulary everything else is written in. |
| Triangle and polygon meshes | 🟢 | `mesh` | `TriMesh`, `PolygonMesh`, named attribute channels, zero-copy views. N-gons stay N-gons until you triangulate. |
| Point clouds | 🟢 | `pointcloud` | `PointCloud`: positions plus optional normal, colour, intensity. No topology, no adjacency, no source-format types. |
| Layered scalar fields (2.5D grid) | 🟢 | `field` | `LayeredField`, row-major cells with a validated configuration. |
| Curves and surfaces | 🟢 | `curves`, `surfaces` | Neutral value types. Representation is not an evaluator claim — see the algorithms table. |
| Profiles and contours | 🟢 | `profiles` | Validated closed/open contour values. |
| Primitive solids and half-spaces | 🟢 | `primitives` | Neutral values with validation. |
| Exact B-rep topology | 🟢 | `brep` | `BRep<Curve3, Curve2, Surface>` — faces, edges, vertices with exact underlying geometry. |
| Authored geometry graph | 🟢 | `model` | `GeometryGraph`, the immutable authored modelling representation. |
| Signed distance fields | 🟢 | `field-ops` | Sampled `Fn(Point3) -> Scalar`, the bridge between points, fields, and meshes. |
| Voxel / dense 3D grid | 🟡 | — | Used transiently inside level-set extraction and B-rep compilation. There is no owned, queryable voxel value type. |
| Tetrahedral meshes | 🟡 | — | Tetrahedra are a *sampling decomposition* inside marching tetrahedra. There is no tet-mesh representation you can hold. |
| Per-point scalar fields as a channel | 🔴 | — | Point clouds carry intensity and nothing more. Arbitrary named scalar channels with gradients/histograms would be new work. |
| Gaussian splats | 🔴 | — | Not implemented. Would be a sibling representation, not an extension of `PointCloud`. |

## Topology and validity

| Capability | Support | Optional | Notes |
| --- | :---: | --- | --- |
| Mesh audit (manifoldness, winding, degeneracy) | 🟢 | `mesh` | `audit_mesh` reports `MeshHealth` counts rather than a pass/fail verdict, so a caller decides what is acceptable. |
| Mesh healing and diagnosis | 🟢 | `heal` | Named defects with the offending element, deterministic order. Diagnosis is separate from repair. |
| Connected components | 🟢 | `mesh` | `component_count`, `decompose`, `compose`. |
| Genus / Euler characteristic | 🟢 | `axiolid-inspect` † | Refuses anything that is not a closed two-manifold instead of returning a meaningless integer. |
| Boundary and non-manifold edge extraction | 🟢 | `core` | Via `EdgeAdjacency`. |
| Containment and winding number | 🟢 | `axiolid-inspect` † | Robust point-in-solid via generalised winding number. |
| Convex decomposition | 🟢 | `axiolid-decompose` † | Splits an arbitrary closed solid into convex parts; volume-conserving, every part a closed solid. |
| B-rep validity audit | 🟢 | `topology` | Face/edge/vertex consistency checks. |
| Non-manifold *repair* | 🟡 | `heal` | Defects are detected and named; automatic repair covers a subset. Some defects are reported for the caller to resolve. |
| Topological simplification (feature removal) | 🔴 | — | Not implemented. Plausible future work for model preparation. |

## Algorithms

| Capability | Support | Optional | Notes |
| --- | :---: | --- | --- |
| Exact orientation / in-circle / in-sphere predicates | 🟢 | `predicates` | Certified signs with degeneracy and filter tests. Depends only on `axiolid-core`, so you can get exact predicates without meshes. |
| Polygon triangulation | 🟢 | `tessellation` | Ear clipping through the tessellation contract. Note: **not** Delaunay. |
| Convex hull (3D) | 🟢 | `generate` | Incremental hull decided by exact predicates; refuses degenerate input by name. |
| Mesh boolean (union / difference / intersection) | 🟢 | `dispatch-mesh-boolean` | Through a swappable provider contract, with evidence and typed refusal. |
| Mesh plane section | 🟢 | `mesh-section` | Exact binary64 plane-side classification with source-topology stitching. |
| Minkowski sum and difference | 🟢 | `axiolid-minkowski` † | Convex sums exactly via hull of pairwise vertex sums; general operands via convex decomposition. Difference is erosion, and refuses a non-convex subject rather than returning a too-large result. |
| Offset and shelling | 🟢 | `generate` | Constant-distance offset with miter handling. |
| Mesh decimation | 🟢 | `axiolid-decimate` † | Edge collapse with a bounded, *reported* deviation; cumulative per vertex so repeated collapses cannot drift. |
| Mesh refinement and smoothing | 🟢 | `axiolid-refine` † | Loop-style subdivision and Laplacian smoothing with pinned boundaries. |
| Level-set extraction | 🟢 | `axiolid-levelset` † | Marching **tetrahedra**, chosen over marching cubes because it is manifold by construction. |
| Point-set reconstruction | 🟢 | `pointcloud-provider` | SDF estimated from samples, extracted as a level set. Composed from capabilities the kernel already owns — no third-party numerics adopted. |
| Ray–mesh intersection | 🟢 | `ray-mesh` | Nearest hit and filtered variants over the BVH. |
| Distance and proximity queries | 🟢 | `measure` | Point-to-triangle, segment-to-segment, mesh-to-mesh distance with witness points. |
| Volume, area, second moments | 🟢 | `measure` | Closed-form over closed meshes, plus an exact path for supported B-rep families. |
| NURBS evaluation and differential analysis | 🟢 | `nurbs` | Curve tangents/curvature, surface fundamental forms, Gaussian/mean/principal curvature. |
| Curve flattening | 🟢 | `evaluate` | Adaptive, with a chord-error budget; fails closed rather than exceeding it. |
| 2D overlay / boolean on regions | 🟢 | `overlay` | Sweep-line arrangement with exact predicates. |
| Linear intersection (2D) | 🟢 | `linear-intersection` | Certified classification with typed refusals naming the operand at fault. |
| Delaunay triangulation | 🟡 | `predicates` | The exact `insphere` predicate exists and is tested against near-degenerate cases. **No triangulator is built on it** — polygon triangulation is ear clipping. |
| 3D curve/surface intersection | 🔴 | — | Not implemented. 2D linear intersection is; the 3D analytic case is not. |
| Normals estimation from raw points | 🔴 | — | We *consume* normals but cannot compute them from a bare point set. This is the main gap blocking good reconstruction from unoriented captures. |
| Cloud-to-cloud / cloud-to-mesh distance | 🔴 | — | The KNN index makes this cheap to add; it does not exist as an operation yet. |
| Point-set subsampling and outlier filtering | 🔴 | — | Not implemented. Standard first step in any scan pipeline. |
| Registration / ICP alignment | 🔴 | — | Not implemented, and a genuine project rather than a small addition. |
| Voronoi diagrams | 🔴 | — | Not implemented. |
| Mesh parameterisation / UV unwrapping | 🔴 | — | Not implemented. |
| Geodesic distance on surfaces | 🔴 | — | Not implemented. |

† Not re-exported through the `axiolid` facade yet — depend on the leaf crate
directly. This is the [ADR 0036](/adr/0036-use-case-specific-compilation-closures)
pattern: the facade is a convenience, never a required route, and a leaf
consumer should not pay for capabilities it does not use.

## Execution and integration

| Capability | Support | Optional | Notes |
| --- | :---: | --- | --- |
| Provider contracts with typed refusal | 🟢 | `contracts` | Every operation that can fail says *why*, by name. An empty result never stands in for "could not". |
| Conformance suites | 🟢 | `contracts` | Exported as library code, generic over the trait — an out-of-tree provider runs the identical checks. Skips are recorded separately, so a provider cannot reach "conformant" by refusing everything. |
| Runtime provider dispatch with fallback | 🟢 | `dispatch-*` | Priority ordering, device matching, and memory-budget admission *before* dispatch. Fallback happens only for `Unsupported`/`Unavailable`. |
| Conformance-gated registration | 🟢 | `dispatch-*` | A non-conformant provider is rejected at registration with its failing report, not discovered later by a caller receiving wrong geometry. |
| Cancellation and scratch budgets | 🟢 | `contracts` | Providers declare their granularity honestly, including "never polls". |
| Deterministic output | 🟢 | `core` | Declared per provider, from best-effort through topological to bit-exact. Routing compares strength, so a stronger provider is never rejected for being too good. |
| C ABI | 🟢 | — | Stable v0.4 surface for native consumers. |
| Parallel CPU execution | 🟡 | `parallel` | The seam exists and some paths use it; most algorithms are still serial. Parallelism is a runtime choice, never a build-for-one-host decision. |
| SIMD | 🟡 | `simd` | Selected at runtime behind the same seam. Coverage is partial. |
| GPU execution | 🟡 | `gpu` | Contract and graph-compile seam exist. Not a claim that we bundle GPU kernels. |
| Python bindings | 🔴 | — | Not implemented. |
| WASM target | 🔴 | — | Untested. The pure-Rust core has no obvious blocker, but no claim without evidence. |

## Deliberately out of scope

These are not gaps. They are decisions, and they are what keeps the kernel a
kernel.

| Not in the kernel | Why | Where it belongs |
| --- | --- | --- |
| LAS / LAZ / E57 / PCD / COPC parsing | Format neutrality. Admitting source-format types is the mistake [ADR 0044](/adr/0044-pointcloud-representation-and-reconstruction) exists to prevent. | A sibling ingestion crate outside `crates/`. |
| STEP / IFC parsing | Same rule, longer standing. | `openbimrs/*`. |
| Rendering and visualisation | A kernel computes geometry; it does not draw it. | Application layer. |
| Scan registration UI, cleanup workflows | Application policy built *on* kernel primitives. | Application layer. |
| Statistical testing, Kriging | Analysis of measurements, not geometry. | Application or a stats library. |
| Navigation meshes | Application semantics on top of geometry. | Application layer. |

## Honest limits

- **"Exact" is scoped.** It applies to the certified predicates and the
  supported exact B-rep families, not to every operation. The
  [capabilities page](/capabilities) draws that line precisely.
- **Scale is untested above ~10⁶ points.** `PointIndex` is a uniform grid,
  appropriate for 10⁵–10⁶ scattered points. We have not measured it at scan
  scale and do not claim it.
- **Reconstruction is an estimate.** A point set does not determine a unique
  surface. Results carry evidence — including how many triangles were
  *interpolated* across gaps in the capture — so a caller can tell measured
  surface from invented surface.
- **Partial rows are partial on purpose.** A 🟡 means the note is load-bearing.
  Read it before depending on the row.
