# Changelog

All notable changes to Axiolid are documented in this file.

## [Unreleased]

### Added

- A capability ledger, `architecture/capability-ledger.toml`, grading 95
  geometry capabilities from OCCT and CGAL against Axiolid, with evidence,
  reference packages at pinned commits and tracking issues (#111,
  #118-#139). `cargo xtask gaps` prints ready work by priority, and
  `gaps show` gives one row or issue. `gaps check` runs in the gate and a
  mutation probe proves each of its rules can fail.
- Capability grade `scoped`, for rows deliberately not raised to
  implemented (a designed refusal, or out of scope with no consumer). It
  requires a written rationale, which the gate enforces; nine rows use it.
  Every other open row now has an issue (#140-#159), with blockers mirrored
  from GitHub. The two design decisions they raised are recorded as
  ADR 0068 (exact numbers: `num-bigint` integers under an owned filtered
  layer, benchmarked in `docs/research/exact-arithmetic-bench/`) and
  ADR 0069 (arc booleans: keep `cavalier_contours` for now, exact path
  later).
- `axiolid-exact` 0.1.0 (#154, ADR 0068): exact constructions over `f64`
  input. An outward-rounded interval filter decides almost every case;
  exact big-integer arithmetic (`num-bigint`) runs only when it cannot.
  Includes `(a + b*sqrt(c)) / d` signs and ordering across radicands,
  exact segment crossings and line/circle hits with exact tangency.
- `axiolid-exact`: nested square roots (`Tower`), exact real roots of
  integer polynomials (`IntPoly`, Sturm), and exact conic intersections
  (`Conic`).
- Exact arc booleans (#155, ADR 0070): `arc_overlay` decides every
  topological question exactly and no longer uses `cavalier_contours`.
  Checked against area identities and point membership on 270 scenes with
  shared edges, tangent and concentric circles, and decimal coordinates,
  plus a twelve-fault mutation probe. A bounding-box broad phase keeps the
  cost close to linear in edge count (two overlapping 256-edge rings: 7 ms,
  was 779 ms before it).
- Exact analytic curve intersection (#119): lines, circles and ellipses
  against each other and against planes, cylinders, cones, spheres and tori,
  with hits, tangency and containment decided exactly on `axiolid-exact`
  polynomials (`axiolid-nurbs`). Checked against a sampling oracle on 600
  random scenes and a six-fault mutation probe.
- Sloped cuts through curved prisms (#120): `clip_arc_prism_exact` cuts a
  column with a sloped plane exactly. The wall stays a cylinder, its cut
  edge is an ellipse, and its trim is the new `Curve2::Sinusoid` (ADR 0071).
- Curved coaxial booleans (#120): `boolean_arc_prisms_exact` now returns
  solids with round through-openings and solids that start above the
  ground plane, instead of refusing them.
- Disconnected coaxial booleans (#120): `boolean_prisms_exact_solids` and
  `boolean_arc_prisms_exact_solids` return one solid per separate piece,
  in a fixed order, instead of refusing.
- Stepped coaxial booleans and caps crossed by a cut (#120, ADR 0072): a
  union of prisms with different heights, a notch, counterbore or blind
  pocket, and a sloped plane crossing a column's top or bottom are now
  exact solids. All are built from one exact arrangement of every input
  ring (`ArcArrangement` in `axiolid-overlay`), so faces share vertices
  exactly. Round walls stay cylinders.
- Authored polygon faces (#160): the reference mesh compiler triangulates
  n-gon, concave and holed `PolygonMesh` faces instead of refusing them.
  Non-planar faces are refused by index.
- Surface models (#161): B-reps with shells but no solid compile to a mesh
  flagged `MeshClosure::Surface`; `CompileOutcome::solid_mesh` refuses a
  volume on it.
- Exact ruled quadric sections (#119, ADR 0076): where a quadric crosses a
  cylinder or cone off any shared axis -- a pipe tee, a column piercing a
  dome off-centre, a pipe entering a hopper, an oblique cut through a
  cone, a pipe bend (torus) meeting a wall or a sphere -- the intersection is
  derived exactly as a root branch of a quadratic over the carrier's angle,
  with its loops and branch ends decided by exact root isolation. Checked
  on 120 random cylinder pairs against a dense scan and a 10-fault mutation
  probe.
- Certified B-rep distance (#125, ADR 0074): `boundary_distance` returns an
  interval certain to contain the distance between two exact B-reps'
  boundaries, and `boundary_clearance` compares it with a limit, answering
  Below, Above or Indeterminate instead of rounding a near-limit value into
  a verdict. Checked on closed-form fixtures (walls, caps, poles, a cone
  apex, a torus crown, a plate with a hole) and an 11-fault mutation probe.
- Exact profile breadth (#111): full-turn revolution of circles, hollow
  circles and rectangles, and any section with holes (each hole a toroidal
  cavity carried as a void shell). Composite profiles are unioned exactly
  over one arc arrangement, so members may carry arcs and their own
  openings, and members that do not touch become separate solids of one
  exact B-rep. The mesh compiler tessellates every solid of a B-rep.
  A translated circle profile lowers to an exact contour instead of being
  refused.
- Cavities (#120): a solid's void shells are tessellated, facing into the
  cavity, so the mesh encloses the outer volume less every cavity. A void
  authored facing the other way is turned round; an open void shell is
  refused. Exact coaxial booleans that bury a tool inside the subject now
  return one solid with a void shell instead of refusing.
- Exact mass properties over curved faces (#125, ADR 0073):
  `exact_properties` integrates cylinders, cones, spheres, tori, elliptical
  cylinders, B-spline faces and arc-bounded planar faces over their own
  parameters instead of refusing them. Checked against closed forms (Pappus,
  parallel-axis moments, hemispheres, cones, a half torus) and an 18-fault
  mutation probe.

### Fixed

- Exact revolutions were built inside out (#125): their surface frames were
  left-handed, so every face pointed into the solid. Both audits passed it;
  measuring the result gave `-2 pi R A`. The frames are now right-handed.
- `axiolid-measure`: `exact_properties` honours face orientation, so an
  exact solid that does not touch `z = 0` measures its true volume (a
  raised unit cube measured 7/3).
- The planar overlay accepts U-shapes and combs: two collinear edges that
  do not touch were misreported as a self-intersection.
- Exact surface intersection no longer mistakes an exact tangency,
  parallel or perpendicular case for a nearby float answer (#119): a
  tangent sphere/plane pair used to yield a circle of radius 1e-7.
- Arc-aware booleans (`arc_overlay`, and `boolean_arc_prisms_exact` on top
  of it) no longer depend on drawing units. The backend's thresholds are
  fixed in drawing units, so a 5 um gap survived a union drawn in
  millimetres but vanished in metres. The drawing is now scaled so those
  thresholds sit at the caller's tolerance (ADR 0069).
- macOS SHARED consumers load `libaxiolid_capi.dylib` again
  ([#113](https://github.com/axiolid/kernel/issues/113)). Its install name
  is rpath-relative, and dyld resolves that only through the consumer's
  own `LC_RPATH`. CMake adds one itself when the dylib exists at configure
  time, which is always true for the installed package, but on a fresh
  source-tree build cargo has not produced it yet, so the consumer got
  none and failed at load. The source-tree target now supplies the dylib's
  directory. A macOS-only mutation removes it and, against an empty cargo
  target dir, requires the consumer to fail at load with dyld's
  `no LC_RPATH` reason, not just any failure.
- Native release archives attach to every release, not only `v0.1.1`.
  The package version was pinned at 0.1.1 and the attach job refused any
  other tag, so no release since v0.1.1 carried native archives. The
  archive now takes the workspace version, so tag, crates.io version and
  archive name agree, and the release-set check requires every archive
  to match the tag.

## [0.3.0] - 2026-09-23

A minor bump because it carries a breaking change, and pre-1.0 Cargo
treats the minor field as the major. Nothing outside this workspace
depends on Axiolid yet, so the break costs nothing to take now and
paying it honestly keeps the version an accurate claim.

Note on the gate: at 0.3.0 every crate's bump already admits breakage,
so `scripts/check-semver.py` reports "nothing to check" rather than
verifying anything. That is correct but weak — the gate regains its
teeth at the next patch release against a published 0.3.0 baseline.

### Added

- `axiolid-core` gained the bounded 2D primitives the kernel was missing:
  `Aabb2`, `Rectangle2` (rotatable, with `from_aabb`), `Triangle2`, and
  `Polygon2`. Signed-area accessors carry the winding rather than discarding
  it, and the shoelace sums are rebased to the first vertex so a polygon in
  georeferenced coordinates does not lose precision to cancellation.
- `axiolid-core` gained the bounded 3D primitives the kernel was missing:
  `Rectangle3`, `Box3` (oriented, unlike the axis-aligned `Aabb`), and
  `Polygon3`. `Triangle3` moved here from `axiolid-field-ops`, which
  re-exports it, so the most reusable 3D primitive no longer requires
  depending on field sampling to reach.

### Added

- Every publishable crate now has its own `CHANGELOG.md`, and `docs/reference/changelog.md` assembles them for the docs site ([ADR 0067](/adr/0067-crates-version-independently)). `scripts/prepare-crate-release.py` bumps one crate and rolls its own changelog; it refuses a bump Cargo's own caret rule would treat as breaking, since that is a workspace-wide event handled by `prepare-release.py` instead. `scripts/assemble-crate-changelogs.py --check` gates drift between the two.

### Changed

- Debug native packages build with `line-tables-only` debug info. Full DWARF
  put the static library at 171.7 MB against the verifier's 128 MiB
  per-member budget, of which 141 MB was debug sections; line tables keep
  file and line numbers in backtraces at 74.7 MB. Packaging now checks the
  budget itself and names the oversized member, rather than leaving the
  failure to surface later as an unexplained verification error.
- `overlay::Polygon` gained `outline()` and `has_holes()`, and `Ring`
  converts to and from `axiolid_core::Polygon2`. `project_mesh` still returns
  `overlay::Polygon`: a projection can have holes and `Polygon2` models a
  simple polygon, so the narrowing is offered as an explicit, documented
  request instead of a silent coercion.
- `minimum_area_rectangle` returns the shared `axiolid_core::Rectangle2`
  instead of a local four-corner struct with a cached area. The edge-vector
  form keeps the parallelogram property by construction, where four
  independently stored corners could be edited into a non-parallelogram and a
  cached area could disagree with them. `OrientedRectangle2` remains as an
  alias, and side lengths moved to a `side_lengths` free function, sorted
  shortest-first so the result does not depend on which hull edge the caliper
  stopped on.
- **BREAKING**: `Mat4` is now `glam::DMat4`, a general 4x4 matrix, instead of
  an alias for the affine `Transform3`. The old alias could not represent a
  perspective projection despite its name, so homogeneous division was
  unavailable under the one type whose name implied it. Code that wanted the
  affine meaning should name `Transform3`; the differing API surface makes
  that a compile error rather than a silent behaviour change.

- **BREAKING**: the `axiolid` facade now defaults to no features
  (`default = []`). A consumer names the capability they need and compiles
  only that; the previous default is available in one line as
  `features = ["standard"]`. Building with no capability feature raises a
  `compile_error!` that names the task-to-feature mapping, the bundles, and
  the provider-vs-contract distinction, so the failure explains itself
  instead of surfacing as "method not found" at every call site
  ([#9](https://github.com/axiolid/kernel/issues/9)).

- Planar-mask dilation and erosion decompose the Euclidean
  structuring element into one contiguous span per row, dropping the
  per-cell cost from O(steps^2) to O(steps): 1.8-2.0x faster on
  128x128 and 256x256 masks. Results are bit-identical, gated by a
  180-case differential test against the previous windowed form
  ([#21](https://github.com/axiolid/kernel/issues/21)).

- Mesh booleans allocate 62% less: 43,028 -> 16,194 allocations on a
  6,912-triangle grouped subtraction (6 -> 2 per triangle). Wall clock
  is unchanged within noise; this buys allocator headroom and a lower
  scratch ceiling, not speed.

### Added

- `MeshBooleanRegistry::with_execution` scopes every dispatched
  provider call to a caller-owned CPU pool instead of rayon's
  process-global one, behind `axiolid-dispatch`'s new `parallel`
  feature. The facade's `parallel` now reaches it, so the feature
  bounds real provider work rather than only adding rayon
  ([#109](https://github.com/axiolid/kernel/issues/109)).

- Mesh booleans are byte-reproducible across runs. Two independent
  sources of run-to-run drift are fixed: `HashMap` iteration order in
  `boolean45`, and `Rc::as_ptr` heap-address tie-breaks in the ear-clip
  comparators, which clipped equal-cost ears in a different order each
  run ([#108](https://github.com/axiolid/kernel/issues/108)).

- The threading page no longer documents dispatch-level pool
  scoping as shipped behaviour. `MeshBooleanRegistry::with_execution`
  and a `parallel` feature on `axiolid-dispatch` do not exist; the
  supported mechanism today is `CpuExecution::install`
  ([#109](https://github.com/axiolid/kernel/issues/109)).

- The ray-index cache carries the same `application` gate as the only
  method that uses it, so a facade build with `ray-mesh` and `spatial`
  but without `application` no longer compiles a cache nothing can
  reach. Found by the new `full` closure profile
  ([#11](https://github.com/axiolid/kernel/issues/11)).

- The native STATIC symbol-mutation probe now proves something on
  every platform: the C consumer promotes implicit declarations to
  errors (`/we4013` on MSVC, `-Werror=implicit-function-declaration`
  elsewhere) and links eagerly on Mach-O, so a removed C ABI symbol
  cannot build quietly on macOS or Windows
  ([#56](https://github.com/axiolid/kernel/issues/56)).

- Four closure profiles the design document named but never
  implemented: `core-only`, `linear-data`, `spatial-rule-checker`,
  and `full`, each with an isolated fixture and a mutation probe
  proving its gate can fail
  ([#11](https://github.com/axiolid/kernel/issues/11)). Eleven
  profiles now gate CI, from a 3-package floor to the maximal
  facade build.

### Fixed

- Splitting a solid and re-uniting the pieces returns ONE solid again,
  not two touching shells
  ([#100](https://github.com/axiolid/kernel/issues/100)). Coplanar
  merging asked whether two faces shared PROVENANCE (same operand, same
  index in that operand's coplanar set) rather than whether they share a
  plane, so a seam between two pieces could never weld. It now falls
  back to the face normals when provenance disagrees. The two
  reconstruction tests are no longer `#[ignore]`d.

### Added

- `scripts/closure-bench.py` measures every closure profile:
  resolved package count, cold-build median over N reps, and
  `target/` size, with a measured noise floor so a gap below it is
  reported as not-a-result rather than ranked
  ([#12](https://github.com/axiolid/kernel/issues/12)). ADR 0036
  gains the full seven-profile table.

- The facade documents and tests the exactness guarantee: no public
  entry point turns exact geometry into a mesh unless the caller asked
  and supplied a tolerance
  ([#36](https://github.com/axiolid/kernel/issues/36)). Audited rather
  than assumed - `tests/exact_primary.rs` pins it, and both a removed
  plane count and a weakened `requires_exact_brep` make it fail.

### Added

- `route::shortest_path_within` takes the vertex budget as a parameter.
  What is affordable depends on the caller's deadline, not on the kernel,
  so a consumer no longer has to shrink its supported model size to adopt
  routing ([#92](https://github.com/axiolid/kernel/issues/92)).
  `MAX_VERTICES` remains the default for `shortest_path`.

### Changed

- **Breaking.** `RouteError::TooManyVertices` gained `budget` and
  `lower_bound` fields, and `RouteError` no longer derives `Eq` because it
  now carries a float. A refusal above the budget reports the
  straight-line distance between the endpoints: a proven lower bound on
  every route, since a polyline is at least as long as the line joining
  its ends and obstacles only lengthen it. A refusal and a bound are
  different facts, and the caller can act on the second
  ([#92](https://github.com/axiolid/kernel/issues/92)).

- `Polyline` reports `ArcLength3d` instead of `Unsupported`: its arc
  length is an exact finite sum of segment lengths, so a distance maps
  to a parameter by a running sum and one linear interpolation
  ([#107](https://github.com/axiolid/kernel/issues/107)). `Ellipse` and
  `BSpline` stay refused. A distance landing on a vertex reads the
  outgoing tangent; a zero-length segment is refused rather than
  normalised; a closed polyline counts its wrap segment but does not
  lap past the total length.

### Added

- Breaking-change policy with a mechanical gate. `docs/contributing/breaking-changes.md`
  states what may change at each version step and what counts as public surface;
  `scripts/check-semver.py` runs `cargo-semver-checks` against the newest published
  baseline below the working version and fails when the public API breaks without a
  bump that admits it (#34).
- Capability evidence gate. `scripts/check-capabilities.py` fails when a row in
  `docs/capabilities.md` claims a status without naming a crate or ADR a reader can
  open; nine rows cited nothing and now cite their implementing crate (#35).

### Changed

- Crate `PLAN.md` files no longer record per-item status. Seven unchecked
  boxes described work that was already implemented, and a contributor
  reading one would have built it a second time. The files that carried
  only a title and a stale status line are gone; the rest keep their
  design rationale and invariants. `check-roadmap-freshness.py` now
  rejects checkboxes and progress headings in `PLAN.md`, so the drift
  cannot return silently (kernel#25).

### Fixed

- The `axiolid-capi` cdylib now carries a SONAME (`libaxiolid_capi.so`)
  and a `@rpath` install name on Mach-O. Without one, a consumer's
  `DT_NEEDED` entry recorded whatever path the linker saw, and ELF never
  resolves a `NEEDED` value containing `/` through `RPATH`/`RUNPATH`, so
  a SHARED-linkage downstream build only loaded from its build tree
  ([#70](https://github.com/axiolid/kernel/issues/70)).

### Added

- The gate's isolated-build list is derived from `cargo metadata` instead of
  hand-maintained, so a new publishable crate is covered the moment it exists.
  Fourteen of 52 publishable crates were escaping the check, including
  `axiolid-curve-evaluate-contract`. A mutation probe registers a throwaway
  crate and fails if it is not picked up (kernel#38).

### Fixed

- Boolean vertex duplication no longer reads past the incidence array
  on meshes where a shared vertex is duplicated more than once. The
  12/21 loops used a result-space bound (`nv_12`/`nv_21`, the sum of
  winding magnitudes) to index source-space arrays, so any winding
  magnitude above one overran and aborted the process (#103).

## [0.2.1] - 2026-09-16

### Added

- Crate names are derived from architecture metadata and checked by
  `cargo xtask architecture check`: a `contract.operation` package is
  `axiolid-<domain>-contract`, a `provider.*` package is
  `axiolid-<domain>-<engine>`, and no other package may end in
  `-contract` or take a name a contract reserves (ADR 0064).
- `scripts/probe_naming_gate.sh` mutation-verifies that rule, including
  decoys that must NOT trip it.

### Fixed

- `capability_ids::ALL` is a `&[CapabilityId]` slice rather than a
  fixed-size array. The length was part of the public type, so
  registering a capability was technically a breaking change; additions
  are now additive.
- `axiolid-pointcloud-reconstruction-contract` and
  `axiolid-tessellation-contract` declared a `domain` that disagreed with
  their own names (`operation.*`); corrected to `pointcloud.reconstruction`
  and `tessellation`. Metadata only, no crate renamed (ADR 0064).

### Added

- `CurveMeasure`: curve evaluation takes either a `Distance` or a native
  `Parameter`, so an authored `IfcParameterValue` cannot be mistaken for
  a length. The parameter route answers every curve family, including
  those whose arc length is refused (issue #106, ADR 0063).
- `axiolid-curve-evaluate-contract`: curve evaluation as a named
  capability, so a consumer can request "a point, tangent or frame at a
  distance" without depending on an engine. `axiolid-evaluate` provides
  `ReferenceCurveEvaluator` (issue #106, ADR 0063).
- `CurveEvaluator::frame_at`: an oriented placement frame using a
  reference-up convention, stable across crests, sags and straights
  where a Frenet frame flips or is undefined (ADR 0063).
- `DistanceConvention`: providers state whether a distance is 3D arc
  length or plan distance, and refuse families where no closed-form arc
  length exists rather than returning a native parameter (ADR 0063).

## [0.2.0] - 2026-09-15
### Performance

Eight changes found by profiling the benchmark suite rather than by
inspection. Every figure below is the median of an interleaved A/B on a
pinned core, taken from the commit that made the change.

- Mesh edge adjacency is laid out as CSR instead of a map of vectors:
  genus 1335 ms -> 543 ms (2.46x). The map spent ~32% of the build in
  malloc/free for one `Vec` per edge.
- Edge records are grouped by counting sort rather than comparison
  sort, in both `audit_mesh` and the adjacency build: audit 514 ms ->
  322 ms (37% faster), and genus a further 590 ms -> 256 ms (2.30x)
  once CSR had made comparison sort the dominant cost.
- Face counts are reused from the build instead of being recomputed:
  genus 2108 ms -> 1309 ms (-37.9%).
- Weld caches hash their keys where ordering cannot be observed:
  levelset 79 ms -> 45 ms (43% faster), refine 147 ms -> 47 ms (68%
  faster). Profiling put 66% of levelset and 77% of refine in the
  midpoint weld cache.
- Decomposition prunes faces that cannot hold the worst concavity:
  646 ms -> 13 ms (50x) on the benchmark mesh, with results identical
  across tolerances from 1e-12 to 2.0.
- Repeated ray casts reuse a cached broad phase: 5227 ms -> 1439 ms
  (3.6x) over a mesh sequence. The index costs ~50 ms to build, so a
  single ray against a fresh mesh is slower than a linear scan; the
  break-even is around 22 rays per mesh.
- The ray-cache key was made cheap enough to stop mattering:
  0.52 ms -> 0.14 ms on 40,962 vertices.
- A caller-held `MeshRayIndex` lets an application keep the index
  across meshes: 4359 ms unindexed -> 369 ms via the application cache
  (11.8x) -> 62 ms held by the caller (70.3x).

### Added

- `CurvatureLaw::shifted`: re-write a law in a coordinate starting at `a`,
  in closed form. Makes trimming a natural-equation curve exact rather
  than a refit (ADR 0062).
- `trim_intrinsic3`, `offset_intrinsic3`, `join_intrinsic3`: relations over
  space curves given by curvature and torsion. Trim and join are exact;
  offset is exact for a helix and refused for a varying law, where the
  offset is not an arc-length curve at all (ADR 0062).
- `Curve3::Intrinsic` is dispatched by `evaluate3`, `derivative3` and
  `domain3`, so existing generic machinery -- graph trimming, composite
  stitching, sweep directrix sampling -- works on a torsion curve without
  special-casing it (ADR 0062).

### Fixed

- Quadrature panels now break at curvature and torsion seams. A panel
  straddling a seam left a joined curve wrong by 3.0e-3, because
  Gauss-Legendre assumes a smooth integrand (ADR 0062).
- A piecewise law can be evaluated over a partial span. Previously any
  seam beyond the requested arc length was refused, which made every
  intermediate evaluation of a joined curve fail (ADR 0062).
- `Curve2::Intrinsic` is dispatched by `evaluate2`, `derivative2` and
  `domain2`. The clothoid type has been storable since 0.1.8 but every
  generic 2D consumer refused it by name, so a transition spiral could be
  built and never evaluated or flattened. `domain2` reports the declared
  ARC LENGTH, not the unit interval, so flattening samples the whole
  curve instead of its first metre.


### Added

- `Curve3::Intrinsic`: a space curve given by its natural equations,
  curvature AND torsion as functions of arc length, anchored to a start
  frame. Carries helices and general space spirals as exact values
  (ADR 0061).
- `axiolid-evaluate::frenet`: frame, point and tangent of a space curve by a
  fourth-order Magnus expansion on SO(3). The frame is orthonormal to
  machine precision at any step size because each step is an exponential of
  a skew matrix, and zero torsion reproduces the planar answer exactly.
- `Intrinsic2::turning_variation_bound`: an upper bound on the total
  variation of heading, for quadrature budgeting.

### Fixed

- Arc-length evaluation of an oscillating curvature law budgeted its
  quadrature from SIGNED turning, which is zero over whole periods of a
  zero-mean law: `k(s) = 2 sin(10 s)` over `[0, pi]` was integrated with one
  panel and landed 2.1e-1 from the true endpoint. Budgeting from total
  variation lands it to 4.9e-15 (ADR 0061).

### Added

- `Curve3::Elevated`: a planar layout paired with an `ElevationLaw`, the
  exact composition an alignment centreline is authored as. The plan keeps
  its own exactness -- including a `Curve2::Intrinsic` transition spiral --
  and the vertical profile keeps its own; neither is approximated to pair
  them (ADR 0060, #105).
- `ElevationLaw`: polynomial and piecewise height laws over PLAN distance,
  with `parabolic` and `constant_grade` constructors for the two vertical
  segment kinds that carry most alignment data.
- Arc-length evaluation of intrinsic curves: `intrinsic_point`,
  `intrinsic_tangent`, `elevated_point`, `elevated_tangent`. Heading is
  exact in closed form; position is Gauss-Legendre quadrature of the
  non-elementary integral, matching a Fresnel reference to better than 1e-9.

### Added

- Overlay hole reachability is pinned by test: `overlay` returns a polygon
  carrying a hole for a difference that encloses a void, so `polygon_area`'s
  hole subtraction is live code rather than an unreachable branch. Mutation-
  verified by deleting the subtraction, which the new area assertion catches.

### Added

- Exact full-turn revolution of any profile that lowers to a contour:
  `Section`, `Contour`, `CenterLine` and `Derived` no longer refuse. Segments
  sweep cylinders, cones, planar annuli and tori according to their own
  geometry; volumes are verified against Pappus. See ADR 0059.

### Fixed

- Reversed cap loops on a revolved annulus carried a forward pcurve interval,
  placing the pcurve start diametrically opposite its 3D edge.
- A torus seam was built as a straight ruling, sagging below the surface by
  `r*(1 - cos(sweep/2))`. It is now an arc around the tube.

### Changed

- `construct` crate docs and `AGENTS.md` no longer claim exact generation is
  limited to rectangle and circle extrusion, or that revolution refuses.
- ADR 0053 (contour holes) and ADR 0057 (section taper) carry supersession
  notes pointing at ADR 0058.
### Added

- Arc extrusion with holes: `extrude_arc_rings` builds one cap face per end
  carrying every ring's loop, with ring 0 outer and the rest through-holes.
  `Profile::Contour` with holes now extrudes instead of refusing, with
  genuine cylindrical walls around curved holes. Hole winding is ENFORCED,
  not demanded: a ring handed either way builds the same solid, because
  `orient_arc_ring` reverses vertices, rotates the bulge assignment and
  negates each bulge. Signed area includes each arc's circular-segment term,
  so the sign is the winding for curved rings too.
- Tapered structural sections: declared flange, web and leg slopes build for
  I, AsymmetricI, T, U and L instead of being refused. Measurement showed
  ADR 0057's stated blocker was wrong -- the existing bisector rounding is
  already tangent to an inclined face at any angle (verified to 1e-9 at 5, 8
  and 14 degrees), so a taper is a different corner list, not new geometry.
  Each tapered face pivots about the mid-point of its run so the declared
  thickness stays the MEAN thickness that section tables state; pivoting
  about the tip would silently change the area. A slope beyond 0.9 of a
  quarter turn is refused as leaving no flange. See ADR 0058.

### Added

- `Profile::Section` exact extrusion for all nine parameterised structural
  variants (I, asymmetric I, L, T, U, C, Z, trapezium). Each lowers to a
  corner ring routed through one shared rounding function, so concave root
  fillets and convex toe radii share a single code path.
- Root fillets are built as exact arcs, not dropped: measured on a 0.4 x 0.3
  I-section with an 0.021 root radius they carry 2.40% of the cross-sectional
  area, so discarding them would leave a section whose area, second moment
  and mass are all wrong while still looking like the right shape.
- Tapered flanges, webs and legs are refused by name rather than silently
  built parallel. A declared `Some(0.0)` slope is a parallel flange and is
  accepted; only a non-zero slope is a taper.

### Fixed

- The section corner router placed the arc centre at `r / sin(turn / 2)`
  instead of `r / cos(turn / 2)`. The two agree exactly at a right angle, and
  every rounded corner reachable through `SectionProfile` is a right angle,
  so the error was invisible to all area tests. Found by testing the router
  directly at a 116.565-degree corner, where the arc came out with radius
  0.01748 instead of 0.02.

## [0.1.8] - 2026-09-07

### Added
- Added `CurvatureLaw::Piecewise` to `axiolid-curve`: several curvature laws over one arc-length domain, tiled by interior seams. This is what lets a straight/transition/arc alignment live in a single `Intrinsic2` under one absolute start frame -- decomposing it into separate curves would require an interior start frame whose origin is the position at the seam, a Fresnel-type integral the representation must not compute. `total_turning` sums each piece over its own subinterval in closed form and refuses (`None`) on a malformed law or a seam outside the curve length rather than clamping; `derivative` is per piece and genuinely discontinuous at seams; `is_straight`/`is_constant` stay structural, with constancy requiring pieces that are constant AND mutually equal. Pieces may nest, so a `Composite` transition can sit inside a `Piecewise` alignment.

## [0.1.7] - 2026-09-06

### Added
- Added `CurvatureLaw::Composite` and `Harmonic` to `axiolid-curve`: a polynomial part plus any number of additive sinusoidal terms in one law, so a transition spiral with both a linear ramp and a sine correction -- `k(s) = k0 + (d/L)s - (d/2pi) sin(2 pi s/L)` -- is stored exactly instead of being refused or approximated. `sine_corrected_transition` derives it from the endpoint curvatures and length. The shape is flat and additive rather than a recursive sum, so a given function has one representation, the family stays closed under differentiation and integration, and `is_straight`/`is_constant` stay structural. Existing variants are unchanged

## [0.1.6] - 2026-09-06

### Fixed

- Coaxiality is now decided by a scale-relative predicate, so the same
  shape derives identically in metres, millimetres and kilometres; frame
  axes are normalised on entry, and a degenerate axis is refused.
- The exact planar-faced boolean gained a contact/tangency lattice test sweep: disjoint, contained, identical, face/edge/vertex touching, and an epsilon ladder from 1e-3 to 1e-15 swept as both a positive gap and a negative overlap. Two limits it found are recorded in the crate PLAN.md rather than hidden: long grid-aligned subtraction chains refuse partway through, and a thin overlap below 1e-12 returns an open shell.
- The exact planar-faced boolean now retries a fixed family of ray directions when a containment probe meets a vertex or edge exactly, instead of refusing the whole operation. Containment is direction-independent, so each attempt stays exact. Grid-aligned operands -- repeated axis-aligned subtraction, where operands share vertices in bulk -- are now answered; previously they refused after a few accumulated operations. Exhausting the family still refuses.
- `axiolid-spatial` no longer claims an octree it does not implement; the crate description and docs name only the BVH and uniform point grid, with the octree listed as a structure that *could* implement the same callback API

### Added
- Added `Curve2::Intrinsic` and `CurvatureLaw` to `axiolid-curve`: plane curves given by their natural equation, curvature as a function of arc length anchored to a start frame. Carries clothoid, Bloss, cubic-parabola, and sinusoidal transition spirals exactly, which no parametric variant can. Exact and symbolic only -- `derivative`, `reversed_orientation`, and `total_turning` are closed form; recovering position needs the Fresnel integral, so it belongs to an evaluator that can state a tolerance, never to the representation
- Added `EdgeAdjacency` to `axiolid-mesh`: edge-to-triangle adjacency derived once, with boundary, non-manifold, inconsistent-winding, vertex-neighbour, and Euler-characteristic queries. `genus`, `smooth`, `decompose`, and heal's orientation unification now ask it instead of each rebuilding their own edge map
- Added `axiolid-pointcloud`: a validated point-set value with optional per-point normal, colour, and intensity channels. Representation only — no topology, no algorithms, and no source-format types; LAS/LAZ/E57/PCD/COPC parsing stays outside the kernel per ADR 0044
- Added KNN and radius queries over point sets to `axiolid-spatial` (`PointIndex`). Callback-based so the hot path allocates nothing per hit, with exact distances rather than the BVH's broad-phase bounds, and deterministic tie-breaking by point index
- Added `axiolid-pointcloud-reconstruction-contract`: request, evidence, typed refusal, and an exported conformance suite for turning point sets into surfaces. A missing surface is always a named refusal, never an empty mesh; surface invented across gaps in a capture is counted in `interpolated_triangles`
- Added `axiolid-pointcloud-reconstruction-sdf`, the reference provider: a signed-distance field estimated from samples and extracted with `axiolid-levelset`. Composed from capabilities the kernel already owns, so no third-party numerics are adopted
- Added `pointcloud-reconstruction` dispatch to `axiolid-dispatch`, with conformance-gated registration. A provider's refusal is returned rather than triggering fallback: it is an answer about the data, and falling through would search for a provider willing to guess
- Added additive `pointcloud`, `pointcloud-queries`, `pointcloud-reconstruction`, `dispatch-pointcloud-reconstruction`, and `pointcloud-provider` features to the `axiolid` facade. Default features are unchanged; `pointcloud` alone resolves to exactly `axiolid` + `axiolid-core` + `axiolid-pointcloud`
- Added `axiolid-minkowski`: `minkowski_sum` for convex planar-faced solids (the convex hull of pairwise vertex sums, exact and needing no boolean solver), `minkowski_sum_with` for arbitrary solids via convex decomposition and union through a `MeshBoolean` provider, and `minkowski_difference_with` for erosion. The difference is computed as the intersection of the subject translated by the negated vertices of the tool -- it is `{ x : x + B subset A }`, not a hull of pairwise differences -- and refuses a non-convex subject by name, since the vertex-wise containment test is only sufficient when the subject is convex. `MinkowskiEvidence` reports part counts, pairwise sums, and boolean operations, so the cost of a decomposition is visible rather than implied.
- Added `axiolid-decompose`: `convex_decompose(mesh, strategy, tolerance)` splits an arbitrary closed two-manifold solid into convex parts. `Strategy::Exact` splits until every part is convex; `Strategy::Approximate { max_concavity }` stops at a stated bound. `Fidelity` reports which was achieved, and the approximate variant reports the concavity actually reached alongside the one requested, so an ask can never be mistaken for an outcome. Concavity is measured against the solid's own face planes rather than its convex hull, since a reflex vertex lies ON its hull and would report zero. Cuts are capped by measuring which edges the clipped shell left used once, rather than predicting the cross-section from the input. `convex_decompose_with` accepts a `Splitter`, so a caller may drive the cut through any `MeshBoolean` provider instead of the built-in clipper; the two are independent implementations of the same contract and are tested against each other.
- Added `axiolid-levelset`: `level_set(field, bounds, edge_length, level)` extracts the level set of a scalar field as a closed manifold `TriMesh`. Decomposes each cell into six Kuhn tetrahedra rather than using marching cubes, whose ambiguous face cases are not manifold without the disambiguated MC33 table: a tetrahedron has no ambiguous sign pattern, so watertightness is structural rather than a property of a table. Vertices are welded by position as well as by edge, without which coincident crossings collapse a triangle and tear the surface. A field that never crosses the level is refused by name rather than answered with an empty mesh. Deterministic for identical inputs. The closed guarantee does not currently hold where the level set is exactly tangent to a grid plane; the crate documents the measured boundary of that case. Simulation of simplicity keeps the result watertight even where the level set is exactly tangent to a grid plane.
- Added `project` to `axiolid-evaluate`: closed-form nearest-point parameters on `Plane`, `Cylinder`, `Sphere`, `Cone` and `Torus`. This is the counterpart to `invert`, which names a point already ON a surface and refuses one that is not. A configuration whose nearest point is genuinely ambiguous -- a point on a cylinder's axis, or at a sphere's centre -- is refused by name rather than resolved to an arbitrary member of the tie. `Cone` projects along the slant, not radially. B-spline surfaces stay with the iterative certified projection in `axiolid-nurbs`.
- Added `axiolid-refine`: uniform and edge-length-driven mesh refinement, plus Laplacian smoothing. Surface-aware refinement places introduced vertices ON the analytic surface the mesh was tessellated from, rather than interpolating between existing triangles, so refining a faceted cylinder converges toward the real cylinder. A surface that refuses to place a vertex fails the refinement instead of silently falling back to linear interpolation, which would return more triangles with none of the promised accuracy. `smooth` holds boundary vertices bit-identical by default.
- Added `AttributeChannel` to `axiolid-mesh`: named per-vertex data (`name`, `values`, `width`, `blend`) carried on `TriMesh::attributes`, with a `Blend` policy (`Linear`, `Nearest`, `None`) that is a property of the DATA rather than of any operation. `validate_structure` refuses a channel that does not cover every vertex, declares a zero tuple width, or repeats a name.

### Changed
- Three sorts in the boolean asked for stability their keys make
  unobservable, and `edge_topology` grew its table by reallocation. A
  union of two 81920-triangle icospheres is ~3.7% faster (median 196.7 ms
  to 189.4 ms, non-overlapping runs) with identical output (ADR 0047).
- The BVH broad phase specialises its traversal per query shape and
  stores identity-free `Aabb` nodes instead of `BBox`, cutting the node
  array from 10.0 MiB to 7.5 MiB at 81920 triangles per operand. A union
  executes 8.9% fewer instructions with 16.5% fewer cache misses; the
  wall-clock effect is smaller than this machine's run-to-run noise, so
  it is reported as counters rather than a speedup (ADR 0047).
- The mesh boolean no longer rebuilds a full `Manifold` for its own
  result: `compute_boolean` returns positions and triangles, keeping the
  empty-result signal and the two-manifold check as explicit steps. With
  a `Copy` `BBox`, an unstable Morton sort, and a correctly sized weld
  map, a union of two 81920-triangle icospheres drops from 264.6 ms to
  201.4 ms (1.31x), measured by interleaved A/B runs with identical
  output checksums (ADR 0047).
- The mesh boolean is now absorbed into `crates/providers/mesh/boolmesh` rather than taken as a `boolmesh` crates.io dependency (ADR 0047, superseding ADR 0014). Profiling put 99.6% of a boolean's runtime inside upstream's single `compute_boolean` call, so neither the hot paths nor the known defects -- including the depth-2 Menger sponge panic -- were reachable from axiolid. Behaviour is unchanged and proven so: a differential test runs union, intersection and difference through both the absorbed algorithm and upstream 0.1.9 and requires bit-identical vertices, equal triangle counts, and volumes agreeing to 1e-12. Upstream's copyright headers are preserved; both projects are MPL-2.0, so absorbing adds no new licence obligation.
- `BooleanEvidence` now reports `attribute_fates`: one `AttributeFate` per named channel on the subject (`Preserved`, `Interpolated`, or `Dropped(DropReason)`). A boolean creates vertices along the cut with no preimage in either operand, so attributes could not always survive -- but they were being dropped SILENTLY, leaving a caller to compare the mesh before and after to discover the loss and with no reason for it. `DropReason` separates `NotBlendable` (the data forbids derivation) from `ProviderLimitation` (this backend does not carry it), so a capability gap does not read as a property of the data. `BooleanEvidence` is no longer `Copy` as a result; it remains `Clone`.


## [0.1.5] - 2026-09-05

### Added
- Added `invert2`/`invert3` to `axiolid-evaluate` (re-exported as `axiolid_reference::curve`): the exact point-to-parameter map for lines, circles and ellipses. Curves had `evaluate`, `derivative` and `jet` but no inversion, so a trim stated as a POINT could not be turned into a parameter. Families with no closed-form inversion are refused by name rather than iterated: introducing Newton here would put a tolerance and a convergence failure mode into every consumer of a point trim, and the certified iterative path belongs to a caller that can carry its evidence. A point off the curve is refused with its residual rather than projected onto the nearest parameter.

### Fixed
- `CurveRelation::Trimmed` with `TrimmingPreference::Cartesian` now resolves. Point selectors were validated and stored but never read: `parameter()` returned `None` for Cartesian, so every point-trimmed curve compiled to `trimmed directrix start needs a finite parameter selector`. Formats that can only stated a trim as endpoints -- a three-point arc knows its endpoints, not their parameters -- were representable but not usable. A basis that is itself a curve relation has no analytic curve to invert against, so a point selector there is still refused, now by a message that says why.

## [0.1.4] - 2026-09-05

### Added
- `SurfaceCurve` now records which p-curve belongs to which surface. `associated_geometry: Vec<NodeId>` was an unordered list that accepted a single entry, a swapped pair, or three unrelated nodes equally; it is replaced by `SurfaceSides`, which pairs each surface with its own p-curve and keeps a single-sided curve expressible.
- `MasterRepresentation::ParameterCurve` is split into `ParameterCurveS1` and `ParameterCurveS2`, so a surface curve can name which parametric side governs. Formats that state the pairing explicitly (such as IFC `PCURVE_S1`/`PCURVE_S2`) no longer have that information discarded at the boundary.
- Graph validation refuses a surface curve whose master names the second parametric side when only one side is present, and checks each side's surface and p-curve against their own node kinds instead of a permissive curve-or-surface test.
- Added `bounded_half_space_in_frame` to `axiolid-construct`: builds a bounded half-space with an explicitly authored in-plane boundary frame. `bounded_half_space` derived the boundary's in-plane axes from the clip plane normal alone, fixing only 2 of 3 orientation degrees of freedom and picking the rotation about the normal by an internal `Vec3::X`/`Vec3::Y` heuristic. Consumers whose source format authors an independent boundary placement (IFC `IfcPolygonalBoundedHalfSpace.Position`) can now supply it. The existing entry point keeps the derived default.
- Added `SpaceFrame` to `axiolid-core`: a validated right-handed orthonormal 3D frame owning `to_local`/`to_world`. Surface evaluation, sampled-field config, and section dispatch each carried their own orthonormality rule under a different tolerance policy; all three now delegate to it.
- Added `PlaneFrame` to `axiolid-core`: an in-plane coordinate frame that validates its basis on construction and owns both directions of the map (`project`, `lift`, `signed_distance`, `normal`). `Frame2` and `Frame3` are inert storage, so every consumer previously re-derived the same two mappings and its own validity rule -- four orthonormality checks existed under three different tolerance policies. Fields are private, so a skewed basis is unrepresentable rather than merely rejected, and the normal is derived from the axes rather than stored alongside them where the two could disagree. `from_normal` covers the case where the caller cares about the plane but not which in-plane direction is x.

### Fixed
- `SolidOperation::BoundedHalfSpace.placement` now orients the boundary profile itself, not just the finished mesh. The compiler passed only the clip plane's normal to the construction, so the boundary's in-plane axes were guessed and an authored placement could not express a rotation about that normal; the transform was applied afterwards, too late to fix the profile's orientation.
- `axiolid-evaluate` now refuses a left-handed surface frame. `invert` checked unit length and perpendicularity but never handedness, so a mirrored basis passed validation and silently produced reflected parameters that still round-tripped through their own bad frame.
- `axiolid-evaluate` now validates surface frames against the caller's tolerance instead of a hardcoded `1e-9`.
- `axiolid-project` decided plane orthonormality with the LINEAR tolerance. Orthonormality is a dot product of unit vectors and therefore dimensionless, so the check scaled with the model length unit: under `Tolerance::MILLIMETRE` a basis skewed by up to 0.5 milliradians passed as orthonormal and every projected coordinate was silently wrong, while the same basis was correctly refused under `Tolerance::METRE`. `Plane` is now an alias for `PlaneFrame`, which decides validity with the ANGULAR tolerance. `ProjectionError::InvalidPlane` is no longer produced and is retained only because removing a public variant is breaking.

## [0.1.3] - 2026-09-05

### Fixed
- `mesh_distance` and `proximity_components` now report zero separation for transverse triangle crossings. The pairwise scan sampled vertex/triangle and edge/edge candidates only, so two surfaces crossing edge-through-face reported their nearest non-intersecting feature instead of zero: a genuine interpenetration read as a real gap, which is fail-open for any consumer asking whether two bodies clash. Coplanar overlap was already caught by the edge/edge family, which is why the existing crossing test did not detect this. A segment/triangle family now runs last, so it cannot disturb the documented tie order of the metric candidates and only ever lowers a result to exactly zero.

## [0.1.2] - 2026-09-05

### Added
- Added `axiolid-route`: exact planar shortest path over a visibility graph, with typed unreachable reasons (`StartOutside`, `GoalOutside`, `Disconnected`) rather than an empty path. Barriers are zero-width polylines, so a wall modelled as a line still blocks a route without bounding area. The kernel reports that no route exists under a given envelope; it never reports that a design is non-compliant.
- Added `axiolid-project`: planar projection of meshes onto a plane, producing a polygon set with holes rather than an outline or hull, plus `intersect_prism` for clipping by a vertical prism. Triangles are unioned pairwise into an accumulator because the planar validator rejects self-intersecting input and a raw triangle soup routinely overlaps itself.
- Added mesh/mesh distance and proximity components to `axiolid-measure`, reporting witness points rather than a bare scalar so a caller can show where the minimum occurs.
- Added `axiolid-inspect`: `min_gap` for clearance and clash detection, `winding_number` and `contains` promoted from the exact boolean's private implementation, `ray_cast`, and `genus`. Containment is decided by certified predicates and is scale-free; `genus` refuses any mesh that is not a closed two-manifold rather than returning a meaningless integer.
- Added `decompose` and `compose` to `axiolid-mesh`, splitting a mesh into connected components and recombining them. Two providers previously counted components and discarded the partition; `component_count` is now a caller of the shared implementation and `boolmesh`'s private union-find is removed. Component order follows first appearance in the input, and a single-body mesh is returned unchanged rather than reindexed.
- Added `offset_solid` and `shell_solid`: constant-distance miter offset and hollowing of planar-faced solids. Offsetting vertices rather than faces handles concave edges, and an offset that closes the gap between opposing walls is refused rather than emitted as a collapsed solid.
- Added `fillet_extruded_profile`: a constant-radius fillet on one straight prism edge, producing a genuine cylindrical blend face tangent to both neighbours. v0.6 refused this by name rather than approximate it with a segmented chamfer.
- Added `boolean_polyhedra_exact`, an exact boolean over general planar-faced solids. v0.6 handled coaxial prisms only; this accepts convex and non-convex operands in any orientation. Containment is exact ray-crossing parity rather than a convex-only plane test, and coplanar contact is resolved by normal agreement so a shared face is kept exactly once. Curved operands are refused by name rather than approximated.
- Added `axiolid-decimate`: edge-collapse mesh decimation to a triangle budget or a maximum deviation, reporting the deviation it actually introduced and refusing collapses that would damage the mesh. A budget still honours the caller's tolerance as a ceiling, so reducing to N triangles cannot licence arbitrary damage. Output is deterministic across runs.
- Added `convex_hull` in `axiolid-construct`: the convex hull of a 3D point set as a closed, outward-oriented mesh. Every face-visibility decision uses certified `orient3d` rather than a tolerance, and degenerate input is refused by name (too few points, collinear, coplanar) so a caller can act on knowing which dimension its input collapsed into. Interior and duplicate points are absorbed without producing degenerate faces.
- Added the `OrientOutward` repair to `axiolid-heal`, flipping a closed shell that encloses negative volume. `UnifyOrientation` makes neighbouring faces agree but leaves the absolute sense as its seed found it, so a consistently-wound inside-out shell passed every structural audit; that is the defect that makes `Difference` behave as `Union` and return a larger mesh with no error. Opt-in like every repair, reported when applied, and a no-op on a shell that is already outward.
- Added `BooleanEvidence::relative_overlap`, reporting how well conditioned the operand configuration was. The production boolean constructs in f64, so accuracy degrades as operands approach coincidence; previously a caller got a badly wrong result with no signal. Reported as evidence rather than policy: Axiolid measures, the caller decides its own refusal threshold. Recorded in [ADR 0045](adr/0045-boolean-construction-arithmetic.md) with the measured cross-kernel comparison and published thresholds.
- Added mesh defect diagnosis to `axiolid-heal`. `Diagnosis`, `Defect` and `DefectKind` were complete vocabulary that nothing produced, and `blocks_boolean` was a predicate over data that never arrived; `diagnose` now produces them from measured evidence, carrying locations where the audit knows them and counts where it does not, rather than inventing an index. `self_intersections` adds a detector Axiolid had no equivalent of: self-intersection is the defect class that most reliably yields a plausible-looking wrong boolean, because a self-intersecting mesh passes every structural audit the kernel has. Crossing is decided through certified `orient3d` on the intersection line of the two triangle planes, never by constructing an intersection point in floating point. Plane-side rejection alone is insufficient and the unit cube proves it: opposite faces each straddle the other plane while missing entirely, which an earlier implementation reported as 8 false pairs on the canonical valid solid. Adjacency is decided on vertex indices before any arithmetic, so triangles sharing an edge or a vertex are never reported as intersecting; without that rule every closed mesh reports as broken. `self_intersections_brute_force` is kept in production code as the reference the BVH-accelerated path is checked against, since an accelerator that changes the answer is a bug. Coplanar overlapping triangles are reported conservatively as intersecting rather than silently dropped, because over-reporting is visible to a caller and a missed self-intersection is not.
- Added the mass-properties provider that `MassProperties` and `Measure<T>` were declared for but never had. `MeshMeasure` implements `Measure<TriMesh>`, `second_moments` computes the `second_moment_diagonal` field that was declared and never populated, and `exact_properties` measures an `ExactBRep` without tessellating it -- so a planar-faced solid measures at machine precision instead of at tessellation fidelity. Previously, measuring an exact solid meant approximating it first and then reporting the approximation's volume as the solid's; v0.6's revolution and chamfer tests hand-rolled divergence sums for exactly this reason. Curved exact faces are refused by name with a message pointing at the mesh path, and open shells are refused rather than assigned the finite-but-meaningless divergence sum of an open surface. Verified against closed forms, the parallel-axis theorem, and an exact-vs-mesh differential that agrees to 1e-9 across two implementations sharing no code.
- Added the `exact` feature to `axiolid-measure`, gating exact-B-rep measurement. It is off by default because `exact_properties` needs `axiolid-brep`, `axiolid-surface` and `axiolid-topology`, which the `mesh-rule-checker` closure profile forbids: a discrete rule checker must not acquire exact B-rep geometry merely to measure a mesh.
- Added an opt-in analytic subtraction path to the `boolmesh` provider: `BoolmeshBoolean::subtract_boxes_analytic` cuts axis-aligned box tools out of an axis-aligned box subject in closed form, bypassing the general mesh boolean for the dominant IFC case (a wall with rectangular openings). Measured ~25x faster than the general solver at 64 openings in the cross-kernel harness. The path is never auto-dispatched: the caller selects it and handles an `Ok(None)` refusal, so output topology never varies for reasons the caller cannot see. It returns `None` rather than an approximate answer whenever an operand is not an axis-aligned box, no tool meets the subject, or the induced grid would exceed the caller's cell budget. Operands are recognised structurally — index count, corner lattice, and per-face-plane triangle distribution — not by bounding box, since every mesh has one of those and a sphere would otherwise be cut as though it were the cube around it. Results are watertight by construction and byte-reproducible across processes, because vertex identity is an ordered map rather than a randomly-seeded `HashMap`.
- Added `BooleanEvidence::analytic_path`, recording whether a closed-form path or the general solver produced a result. It is sticky through `absorb`, so a composed operation cannot report itself as a pure general-solver product. This is the same distinction `sub_operations` already draws between a composed and a primitive result: the two paths produce different (equally valid) triangulations, and a caller reproducing or diffing results needs to know which machinery ran.
- Added `MeshBoolean::determinism`, so a provider declares the reproducibility level it can actually honour. It defaults to the weakest level (`Determinism::BestEffort`), matching the fail-safe shape of `scratch_requirement` and `cancellation_granularity`: overstating is the dangerous direction, because `Plan::admit` refuses a step whose guarantee is weaker than the caller requested and that refusal is only sound if the declaration is honest. The conformance suite gained `determinism_declaration_is_verifiable`, which refuses a `Determinism::Bitwise` claim outright — cross-process byte equality is not observable from a single-process suite, so it must be earned by construction rather than asserted.
- Added an opt-in `parallel` feature to the `boolmesh` provider, wiring upstream's `rayon` multi-threading, which was previously unreachable because the dependency was declared without a `features` key. Off by default: enabling threads changes performance characteristics and pulls in rayon, so a caller should choose it. It does not weaken the declared guarantee — upstream dedups vertices through a randomly-seeded `HashMap`, so the general path's vertex ordering already varies between processes single-threaded, and `Determinism::Topological` (stable connectivity, unstable ordering) is the honest ceiling at any thread count. Callers needing byte reproducibility use `subtract_boxes_analytic`, whose ordered vertex identity earns it.
- Added `Plan` and `PlanStep` to `axiolid-contracts`: a reproducible record of the options an operation ran under, which provider served each step, and what guarantee that provider actually delivered. `Plan::admit` refuses a step whose provider guarantees weaker determinism than the caller requested, using ordering rather than equality so a stronger provider still satisfies a weaker request. Before this, `Determinism` was declared in the contract — including `Bitwise`, documented as the only level supporting cross-machine artifact comparison — but no provider anywhere in the workspace read it, so requesting a level a backend could not meet returned best-effort output with no signal. Provenance records the guarantee delivered, not the one requested, so a plan is evidence of what happened rather than a restatement of intent. Plans are in-process artifacts and deliberately not serialised: a stored plan is a compatibility promise, and freezing a wire format before the API stabilises would commit to a shape the kernel has not finished learning.
- Added exact degree elevation, knot removal, degree reduction, curve interpolation, and lofting to `axiolid-nurbs`. Elevation is exact and always succeeds; it blends in homogeneous coordinates so rational curves elevate correctly, and a polynomial input stays polynomial rather than acquiring a vector of ones. Knot removal and degree reduction are lossy, so each computes a candidate, measures the deviation against the original by sampling, and refuses when it exceeds the caller's tolerance — returning the measured value either way so a caller can check it against their own budget instead of trusting a boolean. `interpolate_curve3` passes through its points rather than approximating them, solving the global interpolation system with partial pivoting because the matrix is banded but not diagonally dominant for arbitrary point spacing. `loft_surface` interpolates down each column of control points, so the surface passes through interior sections and not merely the first and last.
- Added exact revolution (`revolve_profile_exact`) and exact straight-path fixed-reference sweep (`sweep_profile_exact`) to `axiolid-construct`, both wired into `ReferenceExactCompiler`. A rectangle revolved a full turn about a parallel, non-crossing axis yields an annular tube of two cylinders and two annular planes — analytic surfaces, no tessellation. The straight sweep delegates to `extrude_profile_exact` rather than reimplementing it, so the two cannot drift apart.
- Added a chamfer contract and single-edge reference implementation (`chamfer_extruded_edge`) to `axiolid-construct`: Axiolid previously had no fillet or chamfer capability at any tier. A chamfered prism is a prism over a chamfered cross-section, so the cut is applied in 2D and extruded through the existing exact path. Fillets, variable radii, edge loops, and curved edges each refuse by name.
- Widened exact boolean coverage (`boolean_prisms_exact`) from half-space-bounded difference and intersection to arbitrary coaxial prisms with equal z-spans, including non-convex and multi-ring operands. Two prisms sharing a z-span reduce exactly to a 2D cross-section overlay, which `axiolid-overlay` already decides with certified predicates, so no new numerical machinery was introduced. Verified differentially against the `boolmesh` mesh oracle, which shares no code with this path.

### Changed
- Coplanar triangle pairs in `self_intersections` are now decided by exact 2D region logic rather than reported conservatively. Both triangles are projected onto their dominant plane and overlap is decided by `orient2d`: a proper edge crossing, or a strictly-interior vertex. Triangles that only touch -- shared edge, shared vertex, or a vertex resting exactly on an edge -- share boundary and no area, and are no longer reported. This makes the query usable on exact boolean output, where every split face is coplanar with its siblings.
- Compile refusals now name the missing *family* rather than only the operation. The mesh path previously answered every unresolvable solid with one `Unsupported { operation: Sweep }`, making revolution, swept disk, fixed-reference sweep, and sectioned spine indistinguishable — its own comment claimed naming the capability lets a caller register a provider, but naming only `Sweep` tells a caller a sweep failed, not which provider is missing.

### Fixed
- Fixed the graph compile path ignoring the caller's `memory_budget_bytes`. Budget admission existed in `dispatch/section.rs` and `dispatch/boolean.rs` but not on the compile path, whose evaluation cache accumulates one mesh per node with no bound. The compiler honestly reports `ScratchRequirement::Unbounded`, which by contract fits no declared budget, so both `compile_mesh` and `compile_mesh_batch_into` now admit before compiling; gating only the single-root path would have left the batch path an unguarded way in. Callers who set no budget are unaffected.
- Fixed `native/CMakeLists.txt`'s `axiolid-cargo-build` custom target having no `BYPRODUCTS`: a Ninja-generated build (`cmake -G Ninja`) could not link `Axiolid::axiolid_shared`/`Axiolid::axiolid_static` to the artifacts `cargo build` actually produces, failing with `needed by '...', missing and no known rule to make it` even though every artifact was on disk. Unix Makefiles (the implicit default generator, and the only one Axiolid's own `test-native-cmake.py`/`test-native-fetch.py` harnesses ever exercised) tolerates a custom target lacking `BYPRODUCTS`; Ninja's stricter dependency graph does not. `BYPRODUCTS` now lists every artifact path (shared library, static library, and the Windows import library) so both generators resolve the same real build graph.
- Fixed `scripts/test-downstream-consumers.py` merging a probed process's stderr into stdout before `json.loads`-parsing its output (e.g. `cargo metadata`): a first-use rustup/toolchain notice on a cold Windows or macOS CI runner corrupted the JSON stream and failed the black-box downstream consumer gate with a spurious `Expecting value: line 1 column 1` error. stdout and stderr are now captured on separate pipes; both are still shown together on failure.
- Pinned Python 3.12 via `actions/setup-python` in the `rust-consumers` CI job (`.github/workflows/native.yml`): the job previously relied on the Ubuntu runner's default Python 3.10, which lacks the standard-library `tomllib` module (3.11+ only) that `test-downstream-consumers.py` and `tests/downstream/test_downstream_consumers.py` require, so every push and release failed on `ubuntu-22.04` with `ModuleNotFoundError: No module named 'tomllib'`.
- Fixed `scripts/package-native.py` shipping an empty `INTERFACE_LINK_LIBRARIES` for the Windows STATIC package, and the identical gap in the parallel in-tree `native/CMakeLists.txt` source-build install path and `axiolid_static` target (both `AxiolidTargets.cmake.in`-templated install and direct `add_subdirectory`/`axiolid_fetch(... LINKAGE STATIC)` consumption): `axiolid-capi` statically linked pulls in Rust std's Windows system-library dependencies (`kernel32.lib`, `ntdll.lib`, `userenv.lib`, `ws2_32.lib`, `dbghelp.lib`, confirmed via `rustc --print=native-static-libs` on `x86_64-pc-windows-msvc`), which were never declared to STATIC consumers on any of the three affected paths, so every downstream C/C++ project linking `Axiolid::axiolid_static` on Windows failed at link time with unresolved externals.

## [0.1.1] - 2026-09-03

### Added
- Added `scripts/prepare-release.py`, a version-bump and changelog-rollover tool: validates a strictly-forward semver bump, rejects releasing an empty Unreleased section, dates and rolls the Unreleased section into a versioned heading, and bumps both `[workspace.package].version` and every internal `axiolid-*` path-dependency version requirement in `[workspace.dependencies]` so 0.x caret semantics never strand an internal dependency behind the bumped crate it points to.
- Wired the release pipeline into `scripts/gate.sh`: it now runs the release script unit tests, the publish dependency-order plan, and the lock-free bootstrap package preflight for all 38 publishable archives, so a broken release pipeline fails the same gate as any other regression instead of being discovered only at `workflow_dispatch` time.
- Added an executable downstream integration quickstart guide covering narrow Rust crates, the Rust facade, plain C, and C++/CMake, each pinned to an immutable release and backed by a capability/refusal matrix and CI-executed source examples; see [Downstream integration](./guide/downstream-integration.md).
- Added a machine-checked `2d-curves` downstream closure: points, affine transforms, application-owned unit conversion, and 2D curve vocabulary compile exactly `axiolid-core`, `axiolid-linear`, and `axiolid-curve`, while every solid/CSG, mesh, topology, B-rep, facade, and provider package is forbidden and mutation-tested.
- Added mutation-proven black-box compatibility gates that copy Rust leaf/facade and native C/C++ consumers outside the workspace, pin Rust dependencies to one immutable Git artifact, consume only verified native archives through exported CMake targets, execute semantic success and typed-refusal paths, and run on Linux, macOS, and Windows.
- Added `Axiolid::axiolid` CMake integration for immutable source builds and verified native archives, with shared/static selection, deterministic manifests/checksums, Linux/macOS/Windows Debug/Release CI, and an AArch64 Linux cross-build gate.
- Added the versioned `axiolid-capi` native boundary with generated C11 declarations, opaque globally unique handles, explicit ownership transfer, structured context-owned errors, bounded mesh import/export, Boolean and batch operations, audit/bounds/measurement/transform queries, exact-result classification, typed exact refusal, panic containment, and compiled C/Rust smoke coverage.
- Added the supported `axiolid::application` boundary with explicit portable-provider selection, capability inspection, typed operation/provider/tolerance errors, shared Boolean and section conformance checks, and an isolated downstream consumer/closure fixture. The facade exposes validation, measurement, Boolean and batched subtraction, mesh sectioning, ray queries, and strict exact-profile extrusion without leaking concrete providers or silently falling back from exact to mesh.
- Added constructed intersection curves for certified surface/surface and curve/surface queries: `construct_surface_surface_curves` returns exact degree-1 curves for the affine patch family with a deviation bound valid over the whole curve (distance to a plane is affine, so endpoint residuals bound the segment), and `construct_curve_surface_points` returns isolated crossing points rather than manufacturing an extent. Every other case refuses by name via `IntersectionCurveRefusal`, which distinguishes proven `Disjoint` from `Unresolved`; see [ADR 0038](./adr/0038-constructed-intersection-curves.md).
- Added globally certified closest-point inversion for B-spline surfaces: `invert_surface_certified` and `invert_periodic_surface_certified` answer "which parameters name this point" rather than merely "how far away is it". Uniqueness is proven, not assumed: the retained global-minimizer cover must form a single connected region that is also localized to the requested parameter tolerance, so a pole or a self-touching patch is reported as `SurfaceInversionRefusal::Ambiguous` with its rival candidate boxes instead of an arbitrarily chosen representative. On-surface membership is a separate obligation judged against the certified distance LOWER bound, and a structurally sound refusal is an `Ok` verdict rather than an error.
- Planar certified curve/curve intersection now resolves ownership instead of reporting raw candidate boxes: `CertifiedCurveIntersection2::Degenerate` carries per-box `ClassifiedCurveContact2` values, boundary-shared roots are deduplicated and fused so a shared-endpoint crossing is counted once, tangency is proven from interval derivative hulls over the whole box rather than a float sample, and a transverse root sitting exactly on a shared cell edge is reported as the new `CurveIntersectionDegeneracy::BoundaryCrossing` instead of an unexplained unresolved box.
- Added `axiolid-ray-mesh`: narrow-phase ray/triangle-mesh nearest-hit intersection returning parametric distance, triangle index, barycentric coordinates, and a certified front/back/coplanar side. Double-sided by default, deterministic lowest-index tie-breaking for coincident hits, typed refusals for degenerate triangles and invalid indices, and `nearest_hit_among` for composing with the `axiolid-spatial` BVH broad phase. Exposed through the facade as the `ray-mesh` feature.
- Added `axiolid-oracle`, an independent mapped-3D verification oracle for intersection and inversion results. It maps claimed parameter boxes back into model space through the portable scalar evaluator, shares no subdivision or interval machinery with `axiolid-nurbs`, reports the measured 3D deviation on failure, and soundly refutes overstated global minimum distances; see [ADR 0037](./adr/0037-mapped-3d-verification-oracle.md).
- Added focused exact extrusion construction and graph compilation: sharp filled/hollow rectangles retain planar supports, positive-axis filled circles retain a cylindrical support, every exact result carries closed topology/pcurves/native spans, and unsupported exact families fail closed. `ReferenceExactCompiler` memoises successful roots in a per-batch `NodeId -> ExactBRep` cache that cannot contain meshes.
- Added typed `GeomError::UnsupportedInput` diagnostics for capability-supported but input-family-unsupported requests; graph compilation owns and relabels nested construction refusals.
- Extracted `axiolid-linear` (line/segment/ray/polyline values), `axiolid-predicates` (certified exact-arithmetic predicates), and `axiolid-linear-intersection` (certified line/line and segment/segment classification) so a line-query application compiles five internal packages instead of the kernel. Existing paths such as `axiolid_curve::Line2` and `axiolid_reference::orient2d` are preserved by re-export; see [ADR 0036](./adr/0036-use-case-specific-compilation-closures.md).
- Added `cargo xtask architecture closure check|explain`, declared closure profiles in `architecture/closure-profiles.toml`, and an isolated consumer fixture under `tests/consumers/`, making a minimal dependency closure a machine-checked compatibility promise rather than a claim.
- Added facade features `linear`, `predicates`, and `linear-intersection` as convenience routes; direct leaf dependencies remain the smallest closure.
- Extracted `axiolid-evaluate` (analytic and spline curve/surface evaluation, jets, elementary inversion) from the `axiolid-reference` umbrella and rewired `axiolid-nurbs` onto it. A CAD closure drops from 18 to 11 internal packages, no longer compiling mesh, spatial, measure, primitive, or the mesh contracts. `axiolid_reference::curve::*` and `::surface::*` are preserved by re-export; see [ADR 0036](./adr/0036-use-case-specific-compilation-closures.md).
- Added verified closure profiles `mesh-rule-checker`, `parametric-curves`, and `cad-exact` alongside `linear-intersection-minimal`, each with an isolated downstream fixture, plus generated [closure documentation](./architecture/closure-profiles.md) and `scripts/probe_closure_gate.sh` proving every profile gate can fail.
- Added a source-backed geometry concepts guide with accessible Mermaid architecture diagrams, native ASCII STL models that render interactively on GitHub and Pages, contract equations, dark/mobile support, and a mutation-proven diagram-source gate.
- Stable, typed capability IDs for tessellation, mesh Boolean, mesh section, and graph-to-mesh contracts, plus an application- and vendor-neutral `openbim.geometry` claim/evidence boundary.
- Added the Axiolid favicon and a canonical glossary with automatic first-use links and hover/focus definitions.
- Added L1 `axiolid-brep`: strict owned exact B-rep results with separately typed 3D curve, 2D pcurve, and surface catalogs plus explicit native trim intervals. The facade exposes it through the new `brep` feature; see [ADR 0024](./adr/0024-exact-brep-result-contracts.md).
- Adaptive analytic `Curve3` directrix sampling and validated `parameter_range` trimming for sweeps, with dimension-generic chord subdivision shared by the 2D and 3D flatteners.
- Test suites for `axiolid-curve`, `axiolid-surface`, `axiolid-primitive`, `axiolid-profile`, `axiolid-tessellation-contract`, and `axiolid-backend-cpu`, pinning vocabulary contracts, validation refusals, and CPU feature selection.
- Added analytic rational B-spline surface partials and normals, plus bounded conforming support-surface refinement for pcurve-trimmed curved B-rep faces with holes, periodic charts, guarded structured-grid/Earcut seeds, and shared seam vertices.
- Added the format-neutral `axiolid-nurbs` algorithm crate and `axiolid/nurbs` facade feature with analytic second-order differential geometry, explicitly budgeted curve/surface projection, verified closed-curve seam wrapping, exact curve knot insertion/reversal/split/Bézier decomposition, and exact surface U/V insertion/reversal.
- Added outward-rounded global certificates for clamped NURBS point-to-curve projection and curve-pair minimum distance in 2D/3D, with interval-aware homogeneous knot refinement, unresolved minimizer cells, deterministic work budgets, and pre-allocation Cartesian guards; see [ADR 0025](./adr/0025-certified-nurbs-subdivision-oracle.md).
- Added bounded globally certified closest-point projection for open clamped polynomial and positive-rational NURBS surfaces, retaining every possible global-minimizer parameter box and requiring both distance-gap and parameter resolution; see [ADR 0030](./adr/0030-globally-certified-surface-projection.md).
- Added opt-in verified periodic curve views for 2D/3D wrapped evaluation and canonicalized insertion/split parameters without changing neutral evaluator or control-net topology semantics; see [ADR 0031](./adr/0031-verified-periodic-curve-views.md).
- Added explicit cyclic `PeriodicBSplineSurface` U/V/UV schemas, wrapped jets, alias-safe fixed-topology edits, seam continuity orders, and globally certified periodic-domain projection; see [ADR 0032](./adr/0032-explicit-periodic-bspline-surfaces.md).
- Added bounded planar clamped curve/curve root isolation with exact-sign line/point classification and distinct zero-length point contacts, outward-rounded rational derivative bounds, strict-interior Krawczyk proofs for transverse roots, explicit native-parameter resolution contracts, localized structural overlap/endpoint-tangency outcomes, compact parameter-only DFS work items, hard allocation-safe work ceilings, and explicit unresolved singular or boundary boxes; see [ADR 0026](./adr/0026-certified-planar-nurbs-root-isolation.md).
- Added bounded clamped, internally continuous 3D NURBS curve/surface root isolation (internal knot multiplicity `1..=degree`; valid full-multiplicity internal knots are unsupported by this certified query) with outward tensor rational-Bézier refinement, native-span surface partial enclosures, strict-interior 3×3 Krawczyk proofs for isolated transverse roots, three-parameter resolution certificates, shared refinement/search budgets, fallible per-node allocations, retained partial certificates, and explicit unresolved tangential, singular, or boundary boxes; see [ADR 0027](./adr/0027-certified-nurbs-curve-surface-root-isolation.md).
- Added bounded clamped NURBS surface/surface candidate exclusion and complete transverse trace certificates for single-span polynomial affine patches. The path proves affine control-net identities over exact binary64 values, proves normal transversality with outward intervals, certifies both boundary endpoints through the curve/surface oracle, preserves both native parameterizations, uses fixed bounded boundary work, and leaves curved, coincident, tangential, boundary-owned, and multispan cases unresolved; see [ADR 0028](./adr/0028-certified-affine-surface-surface-tracing.md).
- Added topology-aware integration for one-owner certified affine traces: the boundary-owned rectangle becomes two closed analytic trimmed faces sharing the intersection edge, while the containing rectangle records the same edge as an explicit embedded pcurve. The result retains native endpoint boxes and a conservative residual bound, reserves certified construction storage fallibly, and refuses corners, mixed/dual ownership, curved traces, and incomplete queries; see [ADR 0029](./adr/0029-certified-trace-topology-integration.md).
- Added a backend-neutral mesh plane-section contract and portable scalar oracle with exact binary64 plane-side classification, source-topology contour stitching, explicit mesh-approximation evidence, bounded output and scratch, cancellation, and fail-closed coplanar/non-manifold handling; see [ADR 0033](./adr/0033-mesh-plane-section-contract.md).
- Added a format-neutral authored `OpenProfile` graph declaration for conservative bounded-open exact 2D curve paths, with finite/structurally valid curve, 2D-offset, instance, and trim-selector validation, shared-DAG-linear traversal, exact same-endpoint refusal, explicit no-area/no-width semantics, solid-operation exclusion, and curve-evaluation classification; see [ADR 0034](./adr/0034-authored-open-profile-contract.md).

### Changed
- Refreshed root and package agent instructions, research snapshots, downstream
  repository identity, and ADR amendment markers for the nested ownership tree,
  current package names, maintained test pointers, and architecture gates.
- Format-neutral production-source checks now reject both Protobuf vocabulary and `prost` imports, with independent mutation probes; tests remain free to name transports when verifying rejection.
- Documentation now builds with VitePress 2/Vite 8 and an advisory-free locked dependency graph.
- Replaced the mixed `axiolid-kernel` package with `axiolid-guarantees`, `axiolid-contracts`, `axiolid-mesh-contracts`, operation-specific contract packages, and execution-owned `axiolid-dispatch`.
- Split `axiolid-field` values from `axiolid-field-ops`; the facade now exposes additive `field`, `field-ops`, and `field-navigation` tiers.
- Renamed the generic mesh-valued `GeometryCompiler` API to explicit `MeshCompiler::{compile_mesh, compile_mesh_batch, compile_mesh_batch_into}` and the reference package/type to `axiolid-mesh-compile::ReferenceMeshCompiler`.
- Nested Cargo packages by architectural ownership and renamed implementation packages to role-specific `axiolid-reference`, `axiolid-construct`, and `axiolid-mesh-boolean-boolmesh`; downstream manifests and Rust imports must migrate atomically.
- Renamed or removed every pre-reorganisation package that no longer resolves, so a rev-pinned consumer can map each name to its destination: `axiolid-scalar` -> `axiolid-reference` (with exact predicates split into `axiolid-predicates`), `axiolid-boolmesh` -> `axiolid-mesh-boolean-boolmesh`, `axiolid-compile` -> `axiolid-mesh-compile`, `axiolid-tessellate` -> `axiolid-tessellation-contract`, `axiolid-kernel` -> `axiolid-contracts`/`axiolid-guarantees`/`axiolid-dispatch`, and `axiolid-sweep` -> removed, its construction code now in `axiolid-construct`. See the [crate migration guide](./contributing/crate-migration.md).
- `axiolid-mesh-compile` now converts already-triangular `PolygonMesh` faces directly to `TriMesh` without retriangulation, while continuing to refuse n-gons and faces with holes until an explicit tessellation provider is selected.
- `axiolid-construct` defines explicit `GenerationRequest` and `GeneratedGeometry` contracts. Exact B-rep and tolerance-bearing `TessellatedMesh` outputs are separate variants; raw `TriMesh` cannot cross the tessellation result seam, and unsupported exact construction refuses instead of returning a mesh fallback.
- Documented the kernel's direction: Axiolid is striving to be a multipurpose **exact B-rep kernel**, with tessellation as a requested output rather than the model. Surface/surface intersection and geometric inversion are now in scope; see [ADR 0020](./adr/0020-exact-brep-kernel-model.md). Performance work is explicitly parked behind capability work on the roadmap.
- Topology audit and planar B-rep compilation now reject empty loops or outer shells, invalid outer-bound cardinality, undersized bounds, and zero/non-finite-area bounds instead of silently emitting empty or filled geometry; `BRepHealth` exposes dedicated empty-loop and multiple-outer counters.
- Reject malformed compact knot encodings, non-finite controls/frames/derived evaluations, and non-positive rational weights before or during spline evaluation.
- Curve flattening and curved-face boundary/interior tessellation now preserve explicit outer/bound orientation and fail closed on non-finite error metrics or unmet tolerance, depth, segment, per-face, input, and aggregate work limits.
- Extracted the format-agnostic geometry kernel from the Nehirde workspace and renamed its public crate prefix from `geom-` to `axiolid-`.
- Extracted scalar solid generation — profiles, lofts, sweeps, revolutions, extrusion, and bounded half-space clipping — from the L3 DAG compiler into the new L2 `axiolid-construct` crate. `axiolid-mesh-compile` now owns graph traversal, caching, model-driven directrices, and B-rep tessellation only; see [ADR 0023](./adr/0023-solid-generation-is-an-l2-crate.md).

### Fixed
- Removed the pinned `version = "0.1.0"` requirement from the workspace-level `axiolid-oracle` path dependency (test-only, `publish = false`). `axiolid-nurbs` depends on it only as a `dev-dependency`; the stale version requirement made `cargo package`/`scripts/verify-packages.py` fail with `no matching package named axiolid-oracle` because Cargo cannot strip a versioned path dev-dependency the way it strips a versionless one.
- Normalized generated C-header line endings before freshness comparison so the same committed ABI header verifies on Windows and Unix hosts.
- Hardened source-neutrality checks against dependency aliases while allowing comments, and mutation-verified both behaviors.
- Pinned every documentation and release workflow action to immutable commits; repository-wide regression coverage rejects mutable refs.
- Updated the field gate, its 10/10 mutation probe, and nested ownership documentation for the `axiolid-field` value / `axiolid-field-ops` algorithm split, including all facade feature tiers.
- Covered transient crates.io lookup failures in the release tooling.
- Updated package metadata, generated crate links, documentation navigation, and the GitHub Pages base to the canonical `axiolid/kernel` repository.
- Replaced unstable `cargo publish --workspace` with guarded stable child-first publication.
- Added lock-free bootstrap preflight for all 31 publishable source archives and staged exact unpatched upload verification; `xtask` remains excluded.
- Corrected the README MSRV badge from Rust 1.85 to the workspace-required Rust 1.88.
- Solid admission and boolmesh result validation now reject finite-coordinate meshes when signed-volume accumulation overflows or otherwise becomes non-finite, instead of accepting non-finite volume as outward orientation.
- `orient3d` exact escalation now preserves error-free coordinate-difference tails before evaluating cofactors and uses a bounded exact dyadic fallback when finite inputs would overflow or underflow expansion intermediates. Previously it could certify false signs for exactly coplanar inputs and false zero for extreme finite coordinates.

### Removed
- Removed the `axiolid-sweep` crate and the facade's misleading `sweeps` feature and `axiolid::sweep` module. The crate held a single `Sweeper` trait with no implementors, no tests, and no references. Its former construction code is now properly extracted into `axiolid-construct`; [ADR 0021](./adr/0021-capability-seams-live-in-the-kernel.md) is superseded by [ADR 0023](./adr/0023-solid-generation-is-an-l2-crate.md).
