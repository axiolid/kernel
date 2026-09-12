# Changelog

All notable changes to Axiolid are documented in this file.

## [Unreleased]

### Changed

- `SurfacePairSplitUnresolvedReason` now distinguishes a proven refusal
  from an unimplemented one. `NoPartitionExists` is returned when the
  certified trace is shown to slit BOTH patches -- one endpoint strictly
  interior, the other on a boundary side -- which leaves each face simply
  connected, so no pair of closed trimmed faces exists at any tolerance.
  `UnsupportedEndpointOwnership` keeps its former meaning: an arrangement
  this crate has not built yet. Previously both collapsed to the latter,
  so a caller could not tell "stop, the answer is no" from "retry may
  help". `classify` now returns the reason instead of a bare `None`, and
  endpoint/domain resolution failures report `DegenerateRepresentative`
  rather than being mislabelled as an ownership problem. Breaking: the
  enum gained a variant.

### Added
- Exact intersection curves for elementary surface pairs
  (`exact_surface_intersection`). Plane/plane yields a `Line3`,
  cylinder/plane a `Circle3` or `Ellipse3` with semi-axes `r` and
  `r / cos(theta)`, and sphere/plane a `Circle3` of radius
  `sqrt(r^2 - d^2)`. Each curve is derived symbolically from the
  operands, so it is exact rather than fitted, and is usable as
  B-rep edge geometry directly. Degenerate configurations refuse
  with a typed reason instead of returning a degenerate curve:
  parallel planes, a tangent sphere/plane touch, and an
  axis-parallel cylinder section are all explicit refusals.
  Unsupported pairs say so rather than approximating.
- Certified curved surface analysis (`certify_surface_arcs`). A whole-patch
  transversality bound is structurally zero on curved patches, because the
  swept normal's interval hull straddles the other normal, so curved pairs
  were refused wholesale. Transversality is now certified per cell and every
  cell is returned with what was proven about it: `Empty`, `Transversal`,
  `Tangential` (geometry obstructs, more budget cannot help), or
  `BudgetExhausted` (policy stopped, more budget can help).
- `audit_coverage` proves a region set accounts for the whole parameter
  domain using exact `u128` measure arithmetic, so a dropped or
  double-counted region is reported as `Gap` or `Overlap` instead of
  silently shipping an incomplete result. See ADR 0049.
- A certified chord that partitions BOTH patches now splits both, as
  `CertifiedSurfacePairSplit3::DualSplit`: four closed trimmed faces,
  two per input surface, sharing one intersection edge. Previously this
  arrangement was refused. The single-split path and its
  `unsplit_face`/`embedded_curve` fields are unchanged; the dual type
  omits them because no unsplit face exists.
- Certified surface/surface endpoints are de-duplicated. A chord ending
  on a boundary of both patches is discovered once per surface scanned,
  and the two reports denote the same point.
- Certified curve/surface intersection now certifies roots lying exactly
  ON a patch domain edge. Krawczyk proves a root by mapping a parameter box
  strictly inside itself, which a root on a box face can never satisfy, so
  such roots previously subdivided forever and returned `Unresolved` at
  every tolerance. The pinned parameter is now fixed and the reduced
  two-unknown system certified instead, where the root is interior.
  Boundary-restricted certificates only fire on the surface's own domain
  edges, never on faces introduced by subdivision.
- Exact boolean results now carry face provenance. Every output face of
  `boolean_prisms_exact` reports a `FaceName::Fragment` naming which
  operand and which original profile wall it is part of, so a material
  assignment or an edge selection made before the boolean still resolves
  after it. Recovery is geometric rather than bookkept: an exact overlay
  never invents an edge, so every result edge lies on an input edge, and
  `orient2d`'s filtered-then-exact cascade decides which. Both endpoints
  must lie on the supporting line -- matching one accepts an edge that
  merely touches a corner and attributes it to the wrong wall. A face
  with no supporting input edge is left unnamed rather than guessed.

- Persistent structural names for exact B-rep faces and edges
  (`axiolid-brep::name`). `FaceId`/`SurfaceId` are arena positions, valid
  only inside one assembled value, so nothing could refer to a face across
  an operation -- which is why `EdgeSelector`'s only variant was
  `NearestCorner(Point2)` and the fillet had to re-locate its target by
  nearest-point search. `FaceName` records provenance instead of position
  (`Swept`, `Blend(EdgeName)`, `Fragment{operand, source}`, explicit
  `Anonymous`), and `EdgeName` names an edge by the canonically-ordered
  pair of faces meeting there, so naming it from either side gives one
  name. Names compose: `origin()` strips boolean layers, `is_anonymous()`
  checks the whole chain. The exact extrusion names its caps and walls, a
  blend is named after the edge it replaced rather than its own position,
  and `EdgeSelector::Named` resolves through the same path as the
  positional variant -- refusing an unresolvable name instead of snapping
  to the nearest corner. Naming is opt-in per producer: an operation that
  cannot say where a face came from reports `None` rather than a
  fabricated name. 9 tests, 3/3 mutants killed.
- `parallel-batch` feature on the `boolmesh` provider: independent nodes of
  `union_many`'s reduction tree run concurrently. Deliberately separate from
  the existing `parallel` feature, which threads inside a single solve and is
  a measured regression on realistic IFC. Measured 1.5-2.6x on disjoint grids
  of 64-216 solids, saturating by 8 threads; 16 threads buys nothing. The
  ceiling is structural -- per-level cost is flat while parallel width
  collapses 62 to 1, so the last three levels hold 44.7% of the runtime and
  cannot use more than 4, 2 and 1 threads. Amdahl over the measured level
  costs caps this at 3.28x with infinite threads. Off by default; see
  docs/architecture/threading.md.
- `MeshBoolean::union_many` — batch union with a provider-chosen reduction
  order, plumbed through dispatch and the facade. The trait default folds
  left; `boolmesh` overrides it with a balanced pairwise tree, which issues
  the same n-1 booleans but keeps operands small until the final levels.
  Measured 1.6x at 8 solids to 7.5-9.5x at 125 on a disjoint box grid across
  three runs (the ratio grows with n, so it is a complexity difference);
  1.1x to 1.9x on overlapping grids, where operands merge into one growing
  solid. Union is associative and commutative, so unlike `subtract_many`
  there is no disjointness precondition and no correctness cliff.
- `BoolmeshBoolean::boolean_fast` -- an opt-in alternative to `boolean()`
  using a component-wise winding-number classification instead of one
  query per vertex, matching Manifold's own approach. Wins ~6.6% instructions,
  ~7.1% cycles on a sphere union; not the default because a bug in edge-break
  detection would mislabel a whole connected component instead of one vertex.
  See ADR 0047 addendum.
- `axiolid-reference` now answers booleans of INTERPENETRATING solids
  exactly, where it previously refused with `Unsupported`. Three stages:
  the intersection curve (nodes named by source topology, never by
  position), retriangulation of every cut face against that curve, and
  classification of each resulting piece by exact ray parity. Verified
  against hand-computed volumes -- half-offset unit cubes give union
  1.875, intersection 0.125, difference 0.875. `ScalarBoolean` routes to
  it through a new `Arrangement::Interpenetrating` rather than erroring,
  so the decision sits beside the disjoint and nested cases. This also
  gave conformance a real cross-check: the exact and epsilon-tolerant
  implementations must agree on geometry that could not be compared while
  the oracle refused it. `boolmesh` is unchanged and remains the
  production path. Coplanar faces sharing an AREA are still refused --
  see ADR 0046 for why continuing would produce a plausible wrong answer.
- `axiolid-reference` no longer refuses coplanar faces that share a plane
  but no area. Two solids flush in one plane and metres apart along it
  produced sixteen coplanar face pairs and a blanket refusal; buildings
  are full of that shape. A new `coplanar` module clips triangle against
  triangle on exact `orient2d` signs, and only a genuine shared area --
  which a curve cannot describe -- still refuses.
- `docs/architecture/threading.md: the CPU thread-pool model, why
  `boolmesh`'s `rayon` feature is deliberately off (determinism, plus
  measured net regression on realistic IFC with the crossover point),
  and how a caller sizes worker count. Documents a shipped decision
  whose evidence previously existed only in a stale worktree.

<!-- Versions 0.1.1-0.1.8 were renumbered on 2026-09-09. They were
originally tagged 0.4.0, 0.9.0, 0.9.1, 0.10.0, 0.11.0, 0.12.0, 0.13.0 and
0.14.0, with gaps where no release was cut. None of them ever reached
crates.io -- only 0.1.0 did -- so the published numbering never matched the
git numbering. They are renumbered contiguously here so tags, GitHub
releases and this changelog agree, and so the next crates.io release
follows 0.1.0 honestly. GitHub MILESTONES are a separate axis (capability
themes v0.2-v1.0) and are deliberately unchanged. -->

## [0.1.8] - 2026-09-07

### Added
- Added `CurvatureLaw::Piecewise` to `axiolid-curve`: several curvature laws over one arc-length domain, tiled by interior seams. This is what lets a straight/transition/arc alignment live in a single `Intrinsic2` under one absolute start frame -- decomposing it into separate curves would require an interior start frame whose origin is the position at the seam, a Fresnel-type integral the representation must not compute. `total_turning` sums each piece over its own subinterval in closed form and refuses (`None`) on a malformed law or a seam outside the curve length rather than clamping; `derivative` is per piece and genuinely discontinuous at seams; `is_straight`/`is_constant` stay structural, with constancy requiring pieces that are constant AND mutually equal. Pieces may nest, so a `Composite` transition can sit inside a `Piecewise` alignment.

## [0.1.7] - 2026-09-06

### Added
- Added `CurvatureLaw::Composite` and `Harmonic` to `axiolid-curve`: a polynomial part plus any number of additive sinusoidal terms in one law, so a transition spiral with both a linear ramp and a sine correction -- `k(s) = k0 + (d/L)s - (d/2pi) sin(2 pi s/L)` -- is stored exactly instead of being refused or approximated. `sine_corrected_transition` derives it from the endpoint curvatures and length. The shape is flat and additive rather than a recursive sum, so a given function has one representation, the family stays closed under differentiation and integration, and `is_straight`/`is_constant` stay structural. Existing variants are unchanged

## [0.1.6] - 2026-09-06

### Fixed
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
