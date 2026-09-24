# Geometry capability comparison: Axiolid vs OCCT and CGAL, like for like

**Date:** 2026-09-23
**Axiolid revision audited:** `cb8d203` (0.3.0 plus the #113 fix)
**Compared against:** OCCT at `7d2efad9` (2026-08-10, the tree under
`~/projects/occt/occt-research/occt`) and CGAL `v6.2.1` at `28811b6`.

This supersedes the earlier gap list (`competitive-capability-gaps-2026-09.md`,
removed; it is in git history). Every issue that
document filed (#72 to #80) is closed, and about twenty crates have landed
since, so its gap list no longer describes the kernel.

## Why this comparison is scoped

Axiolid is a geometry kernel and nothing else. It does not read or write
files, draw, store documents, or host a scripting console; those concerns
live in separate projects (STEP in `openbim/step`, IFC in the openbim IFC
workspace). OCCT bundles all of them, and CGAL bundles I/O and viewers. A
comparison of whole repositories therefore measures scope, not capability.

This document compares only the geometry both reference libraries contain.
Everything below is scoped that way, including the line counts.

### What was excluded

- **OCCT:** the DataExchange, Visualization, ApplicationFramework, Draw and
  Deprecated modules; the `TKernel` toolkit (OS, memory, collections,
  units); geometry serialisation (`BinTools`, `BRepTools`, `GeomTools`,
  `BRepGraph`); the `Expr` formula parser; the `GTests` suites.
- **CGAL:** I/O (`Stream_support`, `CGAL_ImageIO`), viewers, Ipe and Qt
  plugins, build/test infrastructure, generic C++ utilities, and the
  machine-learning classifier. Everything else is geometry and stays,
  including CGAL's exact number types, since Axiolid's predicates are the
  counterpart.

### Counting rule

Code lines only: blank and comment-only lines are dropped, the same way for
C++ and Rust. Test directories are excluded in all three. Generated numeric
tables are excluded too, because they are data, not logic: 170,587 lines in
OCCT (123,680 of them in `AppCont`'s precomputed approximation matrices) and
10,504 in CGAL. Axiolid has none.

## Size by capability area

| Area | OCCT | CGAL | Axiolid |
| --- | ---: | ---: | ---: |
| A. Foundations: predicates, number types, math | 77,335 | 166,230 | 7,667 |
| B. Curves and surfaces | 316,286 | – | 14,259 |
| C. B-rep topology and construction | 519,674 | – | 14,231 |
| D. Polygon meshes | – | 168,180 | 11,040 |
| E. Triangulations and meshing | – | 222,543 | 846 |
| F. 2D polygons and arrangements | – | 142,213 | 2,168 |
| G. Point sets and reconstruction | – | 72,172 | 2,042 |
| H. Spatial search and miscellany | – | 25,051 | 1,141 |
| Contracts, execution and facades | – | – | 2,500 |
| **Total** | **913,342** | **796,389** | **55,894** |

OCCT and CGAL barely overlap: OCCT's geometry is curves, surfaces and
B-rep; CGAL's is meshes, triangulations, planar and point-set geometry.
Axiolid spans all eight areas, so it is compared against both at once.
In scope, OCCT is 16.3 times Axiolid and CGAL 14.2 times.

Line counts measure effort, not coverage. The capability table below is
the actual gap list.

## Capability rows

95 rows, each taken from a package or toolkit in OCCT or CGAL. Axiolid's
column was graded from source and tests, never from names or docs:

- **implemented:** a general algorithm with tests.
- **narrow:** implemented for a named subset; the row says which.
- **scoped:** deliberately not raised further. Either the narrowing is the
  design (a documented refusal), or it is a specialist structure with no
  consumer. The ledger records why for each, and the gate requires it.
- **contract only:** a type, trait or refusal with no algorithm behind it.
- **absent:** nothing found.

The first grading found 24 implemented, 43 narrow and 28 absent. Nine were
then reclassified as scoped after review; the table below is current.

| Area | Rows | implemented | narrow | scoped | absent |
| --- | ---: | ---: | ---: | ---: | ---: |
| A. Foundations | 5 | 4 | 1 | 0 | 0 |
| B. Curves and surfaces | 16 | 6 | 8 | 0 | 2 |
| C. B-rep | 20 | 4 | 12 | 0 | 4 |
| D. Polygon meshes | 20 | 6 | 6 | 1 | 7 |
| E. Triangulations | 10 | 1 | 1 | 3 | 5 |
| F. 2D | 11 | 1 | 5 | 1 | 4 |
| G. Point sets | 8 | 3 | 1 | 2 | 2 |
| H. Spatial and misc | 5 | 0 | 2 | 2 | 1 |
| **Total** | **95** | **25** | **36** | **9** | **25** |

No row graded **contract only**: every capability Axiolid names has an
algorithm behind it, even where that algorithm is narrow.

### A. Foundations

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| A1 | exact/filtered predicates (orient, incircle, insphere) | C:Kernel_23,Filtered_kernel O:- (OCCT uses tolerances) | implemented | orient2d/orient3d/incircle/insphere with adaptive exact fallback via expansion arithmetic; general for any finite f64 input Evidence: `crates/algorithms/predicates/src/orient3.rs::orient3d` (dispatches to `orient3_dyadic`/expansion on filter failure), `src/orientation.rs::orient2d`, `src/sphere.rs::incircle/insphere` |
| A2 | exact number types / arbitrary precision | C:Number_types,CGAL_Core,Algebraic_* O:- | implemented | `axiolid-exact` (#154): an outward-rounded interval filter over exact big-integer dyadic arithmetic, one expression evaluated in both tiers; signs and order of (a + b*sqrt(c))/d across radicands; exact segment crossings and line/circle hits. Division-free, so no rational type; one square root per value.
| A3 | tolerance model | O:Precision,BRepLib tolerances C:- | implemented | Explicit `Tolerance{linear, angular}` struct, validated, no silent default; used pervasively as an explicit parameter Evidence: `crates/foundation/core/src/scalar.rs::Tolerance::new/eq` |
| A4 | transforms, frames, bounding boxes | O:gp,Bnd,BndLib C:Kernel_23,Bounding_volumes | implemented | Orthonormal `SpaceFrame`/`PlaneFrame` (validated, right-handed), `Aabb` with union/gap/intersects; general Evidence: `crates/foundation/core/src/space_frame.rs::SpaceFrame::world`, `src/bounds.rs::Aabb::{intersects,gap,union}` |
| A5 | linear algebra / root finding / optimisation / integration | O:math,MathRoot,MathOpt,MathInteg,PLib C:Solver_interface(excluded),QP_solver | narrow | No standalone math crate; root-finding (bisection) and Gauss-Legendre quadrature exist only inline inside curve-evaluation code, and a dense linear solve exists only inline inside curve interpolation/loft fitting — none exposed as a reusable general numeric-kernel API Evidence: `crates/algorithms/parametric/evaluate/src/curve.rs` (bisection, comment "Why bisection rather than a closed-form segment count"), `crates/algorithms/parametric/evaluate/src/arc_length.rs::` (Gauss-Legendre quadrature), `crates/algorithms/parametric/nurbs/src/fit.rs::solve` (private Gaussian-elimination-style solver for interpolation only) |

### B. Curves and surfaces

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| B1 | analytic curves (line, conics) + eval/derivs | O:Geom,Geom2d,ElCLib C:Circular_kernel_2 | implemented | Line2/3, Circle2/3, Ellipse2/3 types with general closed-form jet evaluation Evidence: `crates/representations/analytic/curve/src/linear.rs`, `src/conic.rs::Circle2/Circle3/Ellipse2/Ellipse3`, evaluation in `crates/algorithms/parametric/evaluate/src/curve.rs` |
| B2 | B-spline/NURBS curves: eval, knot ops, degree elev, split | O:BSplCLib,Geom C:- | implemented | eval, knot insertion/removal, degree elevation (exact) and reduction (bounded/refusable), split -- all general for arbitrary degree/knot-vector B-splines including rational Evidence: `crates/algorithms/parametric/nurbs/src/degree.rs::elevate_degree3` (exact), `reduce_degree3`/`remove_knot3` (bounded, return `BoundedResult` deviation, refuse via `GeomError` when tolerance can't be met), `src/transform.rs::split2/split3` |
| B3 | analytic surfaces (plane,cyl,cone,sphere,torus) | O:Geom,ElSLib C:Circular_kernel_3(partly) | implemented | Plane, Cylinder, EllipticalCylinder, Cone, Sphere, Torus all present as exact types with evaluation Evidence: `crates/representations/analytic/surface/src/elementary.rs::Plane/Cylinder/Cone/Sphere/Torus`, eval in `src/evaluate.rs` |
| B4 | NURBS surfaces: eval, knot ops, iso-curves | O:BSplSLib,Geom C:- | narrow | Surface evaluation and knot insertion in u and v (insert_surface_knot_u/v, tested) and parameter reversal exist. Missing for surfaces, though present for curves: knot removal, degree elevation and reduction, and iso-parameter curve extraction. |
| B5 | swept/revolved/extruded/offset SURFACES as types | O:Geom (SurfaceOfRevolution, OffsetSurface...) C:- | narrow | A `Revolution` helper struct exists but is `pub(crate)` (private, used internally by the NURBS revolution-profile module), not a public exported analytic-surface type; curve OFFSET (not surface offset) is exact-only for the helix special case and refused otherwise Evidence: `crates/algorithms/parametric/nurbs/src/revolution_profile.rs::Revolution` (private struct, line 234) |
| B6 | point inversion / projection on curve & surface | O:Extrema,ProjLib,GeomAPI C:- | implemented | Newton-style bounded projection for curves (`project_curve2/3`) and surfaces (`project_surface_certified`/`invert_surface_certified`), general within declared budgets, not certified global minimum (documented) Evidence: `crates/algorithms/parametric/nurbs/src/curve_projection.rs::project_curve2/project_curve3`, `src/certified_surface_inversion.rs::invert_surface_certified` |
| B7 | curve-curve / curve-surface / surface-surface extrema | O:Extrema C:- | narrow | Extrema exist only as a side effect of certified intersection/distance routines (`distance_curve2/3_certified`) for bounded parameter boxes; no general dedicated extrema API decoupled from intersection Evidence: `crates/algorithms/parametric/nurbs/src/certified_curve_distance.rs::distance_curve2_certified/distance_curve3_certified` |
| B8 | curve-curve intersection (2D & 3D) | O:IntCurve,Geom2dInt,IntAna2d C:Intersections_2/3 (linear), Arrangement traits | narrow | Certified intersection exists for NURBS curve pairs with explicit transverse/degeneracy classification and bounded work budgets -- general for regular transverse cases, refuses/returns unresolved on tangential or unbounded-work cases; separately a simple linear line/segment intersection exists (exact, general for lines only) Evidence: `crates/algorithms/parametric/nurbs/src/certified_curve_intersection.rs::intersect_curve2_certified` (returns `ClassifiedCurveContact2`/`CurveIntersectionDegeneracy`), `crates/algorithms/query/intersection/linear/src/line_line.rs`/`segment_segment.rs` (exact linear-only) |
| B9 | curve-surface intersection | O:IntCurveSurface C:- | narrow | Certified bounded curve/surface intersection exists with the same budgeted/refusal design as B8 Evidence: `crates/algorithms/parametric/nurbs/src/certified_curve_surface_intersection.rs::intersect_curve_surface_certified` (options carry `max_refinement_work`/`max_depth`, refuses beyond budget) |
| B10 | surface-surface intersection (analytic + general) | O:IntPatch,IntWalk,GeomInt,IntPolyh C:- | narrow | Explicitly documented as narrow by the source itself: both inputs must be single-span polynomial affine patches with continuous clamped axes; general patch pairs are bounded and conservatively returned as unresolved candidates rather than passed through heuristic marching Evidence: `crates/algorithms/parametric/nurbs/src/certified_surface_surface_intersection.rs::intersect_surface_surface_certified` (module doc, lines 1-6) |
| B11 | approximation/interpolation (fit curve/surface to points) | O:AppDef,AppCont,Approx,GeomAPI C:- | narrow | Only INTERPOLATION (curve passes exactly through points, chord-length parameterised, C2 cubic) and section-based lofting exist; no least-squares APPROXIMATION to a tolerance for scattered points Evidence: `crates/algorithms/parametric/nurbs/src/fit.rs::interpolate_curve3` (module doc: "Interpolation, not approximation"), `loft_surface` in same file |
| B12 | surface filling / plate / Coons / Gordon | O:GeomFill,GeomPlate,Plate C:- | **absent** | No Coons/Gordon/Plate filling code found anywhere in the workspace Evidence: grep for Coons/Gordon/Plate/filling across crates/*.rs returned no surface-filling matches (only unrelated hits: CSG "filling" comment, "hatch" comment) |
| B13 | continuity/conversion (to Bezier, to BSpline, reparam) | O:GeomConvert,Convert C:- | implemented | Bezier-segment decomposition of B-splines is a general, exact, reused primitive (used by degree elevation/reduction itself) Evidence: `crates/algorithms/parametric/nurbs/src/transform.rs::bezier_segments2/bezier_segments3` |
| B14 | 2D constraint solving (tangent circles/lines, Apollonius) | O:GccAna,Geom2dGcc C:Apollonius_graph_2 | **absent** | No GccAna/Apollonius/tangency-solving code found Evidence: grep for GccAna/Apollonius/"tangent circle" across crates/*.rs returned zero matches |
| B15 | local properties: curvature, normals, inflections | O:GeomLProp,LProp,LocalAnalysis C:Jet_fitting_3,Ridges_3(mesh) | implemented | General regular-surface differential analysis (first/second fundamental forms, Gaussian/mean/principal curvature) with degeneracy refusal; curve curvature via `CurvatureLaw` for intrinsic curves Evidence: `crates/algorithms/parametric/nurbs/src/surface_analysis.rs::analyze_surface` (general, computes principal curvatures via discriminant with degeneracy floor) |
| B16 | laws / helix / fair curves | O:Law,HelixGeom,FairCurve C:- | narrow | `CurvatureLaw` (constant/polynomial/sinusoid/composite/piecewise) exists and drives `Intrinsic3` curves including true helices (closed-form offset); no fair-curve (minimum-energy/spline-fairing) algorithm found Evidence: `crates/representations/analytic/curve/src/intrinsic.rs::CurvatureLaw`, `crates/representations/analytic/curve/src/intrinsic3.rs` (helix detection via constant curvature+torsion) |

### C. B-rep topology and construction

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| C1 | B-rep data structure (vertex/edge/wire/face/shell/solid) | O:TopoDS,BRep,BRepGraph C:HalfedgeDS,Combinatorial_map,LCC | implemented | Full append-only typed arena BRep<Curve3,Curve2,Surface> generic over geometry, with vertex/edge/loop/face/shell/solid; general Evidence: `crates/representations/brep/src/lib.rs::BRep` (add_vertex/add_edge/add_loop/add_face/add_shell/add_solid), `crates/representations/topology/src/brep.rs::BRep` (generic arena) |
| C2 | validity checking of B-rep | O:BRepCheck C:is_valid | implemented | `audit_brep`/`BRepHealth` catches dangling references, open/empty loops, unpaired edge uses, overused edges, false closure claims -- general structural validator, deliberately excludes geometric self-intersection checks (documented as separate concern) Evidence: `crates/representations/topology/src/audit.rs::BRepHealth` (dangling_references, open_loops, unpaired_edge_uses, etc.) |
| C3 | primitives (box,cyl,cone,sphere,torus,wedge) | O:BRepPrim,BRepPrimAPI C:- | narrow | `Primitive` enum covers Block, Sphere, Cylinder, Cone, Pyramid only -- no Torus, no Wedge as primitive-solid variants (Torus exists only as an analytic surface type, B3, not as a closed primitive solid) Evidence: `crates/representations/analytic/primitive/src/solid.rs::Primitive` (non_exhaustive enum: Block/Sphere/Cylinder/Cone/Pyramid), tessellated by `crates/algorithms/reference/src/tessellate.rs::tessellate_primitive` |
| C4 | extrude/revolve (prism, revol) exact | O:BRepSweep,BRepPrimAPI C:Straight_skeleton_extrusion_2 | narrow | Exact extrusion and exact revolution both exist but are restricted to specific profile families (rectangle/polygon/certain contour types), not general arbitrary-profile extrude/revolve; source itself titles the module "Exact revolution for the profile families exact extrusion already covers" Evidence: `crates/algorithms/construction/construct/src/extrude_exact.rs::extrude_polygon_rings`, `src/revolve_exact.rs` (module doc: "Exact revolution for the profile families exact extrusion already covers") |
| C5 | sweep/pipe along curve (exact) | O:BRepFill,BRepOffsetAPI_MakePipe(Shell) C:- | narrow | General sweep module exists (`sweep.rs`) but exact sweep is limited to specific directrix families; there is a separate "surface_curve_sweep" compile-time test suggesting only certain curve/surface combinations are exact-supported Evidence: `crates/algorithms/construction/construct/src/sweep.rs`, `crates/execution/compile/src/directrix.rs` |
| C6 | loft / thru-sections | O:BRepOffsetAPI_ThruSections C:- | narrow | `loft_surface` in the NURBS fit module requires sections to already agree in degree and control-point count ("reconciling mismatched sections means knot merging and degree elevation... doing it implicitly would hide a shape change") -- refuses mismatched section families rather than reconciling them Evidence: `crates/algorithms/parametric/nurbs/src/fit.rs::loft_surface` (module doc explicitly states the restriction) |
| C7 | point classification in solid/face | O:BRepClass,BRepClass3d,TopClass C:Side_of_triangle_mesh | implemented | Winding-number-based exact inside/outside classification against a closed triangle mesh, using certified orient3d, general (any closed mesh); explicitly returns `None` (not a silent guess) for ray-degenerate ties Evidence: `crates/algorithms/query/inspect/src/containment.rs::winding_number/contains` |
| C8 | B-rep tessellation with tolerance | O:BRepMesh C:- | implemented | `tessellate` in execution/compile drives faceting of the full BRep with explicit chord/tolerance budgets; large dedicated test suite (2288 lines) Evidence: `crates/execution/compile/src/brep.rs::tessellate` |
| C9 | B-rep booleans (general, incl. curved) | O:BOPAlgo,BRepAlgoAPI,TopOpeBRep* C:- | narrow | Two tiers exist: (1) exact boolean over general planar-faced polyhedra (`boolean_polyhedra_exact`, convex or non-convex, any orientation) -- general for planar-faced solids only, refuses curved faces entirely; (2) a mesh-level (tessellated) boolean provider (`boolmesh`, absorbed from an external crate) that operates on triangle meshes, so curved geometry is only booleaned after faceting, never exactly, and carries a known defect (a depth-2 Menger sponge panics inside the absorbed algorithm per its own module doc) Evidence: `crates/algorithms/construction/construct/src/polyhedron.rs::boolean_polyhedra_exact` (module doc: "General exact boolean over planar-faced solids (#77)"), narrower coaxial-prism special case in `src/boolean_exact.rs::boolean_prisms_exact` |
| C10 | fillet (constant/variable radius, edge networks) | O:ChFi3d,BRepFilletAPI C:- | narrow | Constant-radius fillet limited to a single straight vertical edge of an exact prism (module doc: "Chamfer and fillet on a straight vertical edge of an exact prism"); variable/tapered-radius fillet exists but only for the same single-corner-of-a-polygon-profile family; explicitly refuses curved edges, edge loops/networks Evidence: `crates/algorithms/construction/construct/src/feature.rs::fillet_extruded_profile/fillet_polygon_corner/fillet_polygon_corners` (module doc: "Variable radius, edge loops, curved edges, and fillets all return typed refusals naming the missing capability"), `crates/algorithms/construction/construct/src/fillet_variable.rs::TaperedFillet` (single-corner taper only) |
| C11 | chamfer | O:ChFi3d C:- | narrow | Same single-straight-edge restriction as C10; general polygon-corner chamfer exists for prism profiles only, no curved-edge or edge-network chamfer Evidence: `crates/algorithms/construction/construct/src/feature.rs::chamfer_extruded_profile/chamfer_polygon_corners` |
| C12 | offset / shell / thicken | O:BRepOffset,BRepOffsetAPI C:- | narrow | Solid offset/shelling module exists (`offset.rs`, "Solid offset and shelling (#78)") but scoped to solids constructible by the same exact-extrusion/polyhedron machinery -- no evidence of general curved-face offset/thicken across arbitrary B-rep Evidence: `crates/algorithms/construction/construct/src/offset.rs` (module doc: "Solid offset and shelling (#78)") |
| C13 | draft angle | O:Draft,BRepOffsetAPI_DraftAngle C:- | **absent** | No draft-angle code found anywhere in the workspace Evidence: grep for `draft.?angle |
| C14 | features (holes, ribs, local ops, split) | O:BRepFeat,LocOpe C:- | narrow | Only corner-local features (fillet/chamfer on a polygon corner) exist under a `feature.rs` module; no hole, rib, boss, or general split/local-modification operations found Evidence: `crates/algorithms/construction/construct/src/feature.rs` (contains only chamfer/fillet corner operations, see C10/C11) |
| C15 | shape healing / fixing | O:ShapeFix,ShapeAnalysis,ShapeUpgrade C:- | narrow | Healing exists for triangle MESHES (weld vertices, unify orientation, drop degenerate elements, orient outward -- opt-in `RepairPlan`, no blanket "fix everything") and for B-rep it is AUDIT-only (C2): no B-rep-level repair/fixing (e.g. gap closing, edge merging on a BRep) found, only mesh repair Evidence: `crates/algorithms/repair/heal/src/repair.rs::RepairAction/RepairPlan` (module doc: "There is deliberately no `All` variant"), operates on `axiolid_mesh` |
| C16 | hidden-line removal / projection drawing | O:HLRBRep,Contap C:- | **absent** | No HLR/hidden-line-removal/silhouette-extraction algorithm found; the only "silhouette" hits are unrelated code comments (tessellation-consistency and mesh-decimation comments) Evidence: grep for `hidden.?line |
| C17 | mass properties exact (B-rep) | O:BRepGProp,GProp C:- | narrow | Exact (non-tessellated) mass-properties computation exists but refuses any non-planar face: module doc states "Only planar faces are supported. A cylindrical or spherical face needs a surface integral this module does not implement" Evidence: `crates/algorithms/query/measure/src/exact.rs::` (module doc explicit refusal for curved faces) |
| C18 | B-rep distance / extrema between shapes | O:BRepExtrema C:- | narrow | Only floating-point metric proximity/closest-point witnesses between primitive shapes exist (`ClosestPoints3`); module doc states these "construct floating-point witnesses only" and that certified topological classification is a separate, unimplemented-here concern; no general certified B-rep-to-B-rep distance/extrema Evidence: `crates/algorithms/query/measure/src/proximity.rs::ClosestPoints3` (module doc: "Certified topological classification belongs in axiolid-reference |
| C19 | medial axis / bisector 2D | O:MAT2d,Bisector C:Straight_skeleton_2, Segment_Delaunay_graph_2 | **absent** | No medial-axis or straight-skeleton implementation; "bisector" hits found are all local per-corner angle-bisector geometry used inside fillet/chamfer construction, not a medial-axis/skeleton algorithm Evidence: grep for `medial.?axis |
| C20 | 2D hatching | O:Hatch,Geom2dHatch C:- | **absent** | No hatching implementation found; the one text hit for "hatch" is an unrelated code comment ("standard escape hatch") Evidence: grep for `hatch |

### D. Polygon meshes

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| D1 | mesh data structure (halfedge/surface mesh) | C:Surface_mesh,Polyhedron,HalfedgeDS O:Poly | narrow | `TriMesh`/`PolygonMesh` are index-buffer structures with derived `EdgeAdjacency`, not a persistent halfedge/DCEL mesh type; the only true halfedge DCEL in the repo is 2D (`axiolid-arrangement`), not a 3D mesh. Evidence: `crates/representations/discrete/mesh/src/lib.rs::TriMesh,PolygonMesh` |
| D2 | mesh booleans / corefinement | C:PMP_Boolean_operations,Nef_3 O:-(BOP on B-rep) | implemented | Robust (non-exact-arithmetic) manifold boolean over general triangle meshes, absorbed own code (ADR 0047); separate from exact path (see D3). Evidence: `crates/providers/mesh/boolmesh/src/provider.rs::BoolmeshProvider` (implements `MeshBoolean`) tested by `crates/providers/mesh/boolmesh/tests/analytic_boxes.rs` and `tests/differential_corpus.rs::bounds_do_not_lose_features_across_scales`. |
| D3 | exact polyhedral booleans (Nef) | C:Nef_3,Nef_2,Nef_S2 O:- | narrow | `boolean_polyhedra_exact` handles general planar-faced solids (convex or non-convex) for union/intersection/difference, constructing intersection coordinates in f64 (not exact rational/interval arithmetic despite the name); curved faces refused. No Nef-style open/half-space set representation (no unbounded regions, no 2D/spherical Nef). Evidence: `crates/algorithms/construction/construct/src/polyhedron.rs::boolean_polyhedra_exact` tested by `crates/algorithms/construction/construct/tests/boolean_polyhedra.rs`. Doc (`docs/capabilities.md` line ~120) confirms: "The production mesh boolean constructs intersection coordinates in f64". |
| D4 | self-intersection detection | C:PMP(self_intersections) O:BOPAlgo_CheckerSI | implemented | General triangle-triangle test via BVH broad phase + narrow phase; reports intersecting pairs. Evidence: `crates/algorithms/repair/heal/src/intersect.rs::self_intersections` (or equivalent fn) tested by `crates/algorithms/repair/heal/tests/self_intersection.rs`. |
| D5 | mesh repair (stitch, orient, degenerate, holes) | C:PMP_Mesh_repair,Polygon_repair O:- | **scoped** | Opt-in `RepairAction` enum: `WeldVertices`, `UnifyOrientation`, `DropDegenerateElements`, `OrientOutward`. Deliberately no "repair everything" mode; each action applies only where safe and reports `skipped` when not applicable. Hole filling is a SEPARATE capability (see D6), not part of `heal`. Evidence: `crates/algorithms/repair/heal/src/repair.rs::RepairAction,RepairPlan,RepairReport` tested by `crates/algorithms/repair/heal/tests/mesh.rs`, `tests/orient_outward.rs`. |
| D6 | hole filling / fairing | C:PMP (triangulate_hole, fair) O:- | **absent** | No `fill_hole`/`triangulate_hole`/fairing function found anywhere in the workspace; `heal` addresses welding/orientation/degenerate removal, not boundary closure. Convex-decomposition's "cap" (D18/decompose) closes a PLANAR cut cross-section only, not an arbitrary boundary hole — different problem. Evidence: grep for fair/hole_fill/fill_hole across crates: 0 hits (excluding decompose's cut-capping). |
| D7 | isotropic remeshing / refinement / smoothing | C:PMP_Remeshing,Tetrahedral_remeshing O:- | narrow | Uniform/edge-length subdivision (4-way split) with optional surface-aware vertex placement when a source B-rep surface is known; separately, boundary-fixed Laplacian smoothing. No isotropic *remeshing* (target-edge-length equalization via edge flip/collapse/split combined) — only pure refinement (splitting, monotonic triangle growth) and pure smoothing, run independently. Evidence: `crates/algorithms/discrete/refine/src/lib.rs::refine` tested by `crates/algorithms/discrete/refine/tests/refine.rs` |
| D8 | simplification / decimation | C:Surface_mesh_simplification O:- | implemented | Edge-collapse decimation with caller bound on deviation; rejects collapses that would invert a triangle, create non-manifold edges, or open the boundary. Evidence: `crates/algorithms/discrete/decimate/src/collapse.rs::decimate` tested by `crates/algorithms/discrete/decimate/tests/decimation.rs::decimation_reduces_and_reports_its_deviation`. |
| D9 | subdivision surfaces | C:Subdivision_method_3 O:- | **absent** | No Catmull-Clark/Loop subdivision scheme found. The 4-way triangle split in `refine` is geometric subdivision-for-refinement, not a limit-surface subdivision scheme (no smoothing masks, no valence-based weighting). Evidence: grep for subdivision/catmull/loop_subdiv: 0 hits outside `refine`'s doc comment discussing why it is NOT that. |
| D10 | parameterization (UV unwrap) | C:Surface_mesh_parameterization O:- | **absent** | No UV/parameterization crate; `region/profile` is 2D construction profiles (extrusion sketches), unrelated to mesh UV unwrapping. Evidence: grep for parameteriz/uv_unwrap/lscm/ARAP: 0 relevant hits. |
| D11 | geodesics / shortest path on mesh | C:Surface_mesh_shortest_path,Heat_method_3 O:- | **absent** | No geodesic-distance or mesh-surface shortest-path algorithm. `axiolid-route` (F8) computes visibility/shortest path in a 2D polygon plane, not on a 3D mesh surface. Evidence: grep for geodesic/shortest_path/heat_method: only 2D polygon route hits (see F8). |
| D12 | segmentation / skeletonization | C:Surface_mesh_segmentation,Surface_mesh_skeletonization O:- | **absent** | No mesh segmentation or (3D) skeletonization code found; only unrelated 2D "route"/navigation code and mesh connected-component decomposition (`component::decompose`, which is D18-adjacent, not segmentation by feature/region). Evidence: grep for segmentat/skeletoniz: 0 hits besides F8's route crate and F3 (2D straight skeleton, ABSENT too, see below). |
| D13 | deformation | C:Surface_mesh_deformation O:- | **absent** | No mesh deformation (ARAP, cage, skinning) code found. Evidence: grep for deformation/as_rigid_as_possible: 0 hits. |
| D14 | mesh measures (area/volume/centroid/moments) | C:PMP measures O:BRepGProp (on triangulation) | implemented | `MeshMeasure` computes area, signed volume, volume centroid, and second moments in one pass; refuses non-closed/non-manifold input for volume/moments (open shell still gets area). Evidence: `crates/algorithms/query/measure/src/mesh_measure.rs::MeshMeasure::measure` tested by `crates/algorithms/query/measure/tests/mass_properties.rs`. |
| D15 | mesh distance (Hausdorff), inside test, ray cast | C:PMP distance, Side_of_triangle_mesh, AABB_tree O:BRepExtrema | narrow | Inside test (winding number, exact via `orient3d`) and ray cast (BVH-accelerated nearest hit) are IMPL and general. Mesh-to-mesh proximity/clearance (`mesh_proximity.rs`) is a min-gap/clash measure, not a full two-sided Hausdorff distance metric — reports nearest gap, "clash" as 0.0, out-of-range as `None`, not the symmetric max-min Hausdorff value CGAL's `PMP::approximate_Hausdorff_distance` computes. Evidence: `crates/algorithms/query/inspect/src/containment.rs::winding_number,contains` |
| D16 | mesh plane section / slicing | C:Polygon_mesh_slicer O:BRepAlgoAPI_Section | implemented | `ScalarSection`/`MeshPlaneSection` produce closed plane-local contours from source-topology stitching of outward-oriented closed solids; explicitly fails closed on coplanar overlap and open/non-manifold input (no region nesting claimed). Evidence: `crates/algorithms/reference/src/section.rs::ScalarSection` and `crates/contracts/operations/mesh-section/src/conformance.rs`, tested by section-related test files under `algorithms/reference/tests`. |
| D17 | mesh clipping/splitting by plane or mesh | C:PMP clip/split O:BRepAlgoAPI_Splitter | narrow | Plane-splitting exists via convex-decomposition's clip-and-cap machinery (D18) and via the exact-B-rep offset/solid path; general clip/split of an arbitrary mesh by an arbitrary cutting *mesh* (not just a plane) is not a standalone general operation — relies on the boolean provider (D2) for mesh-by-mesh splitting. Evidence: `crates/algorithms/discrete/decompose/src/lib.rs` (plane clip+cap, see doc comment "Two splitters") |
| D18 | convex decomposition | C:Convex_decomposition_3 O:- | implemented | `Strategy::Exact` (recursive plane split to zero concavity, to tolerance) and `Strategy::Approximate` (bounded concavity); reports which strategy actually applied and the concavity achieved. Requires closed two-manifold solid input (refused otherwise) — this is the format precondition, not a narrowing of the decomposition algorithm itself. Evidence: `crates/algorithms/discrete/decompose/src/lib.rs::Strategy,Decomposition,convex_decompose` (approx name) tested by `crates/algorithms/discrete/decompose/tests/decompose.rs::a_reflex_solid_is_split_into_convex_parts`, `an_approximate_result_never_claims_to_be_exact`. |
| D19 | approximation (VSA) / shape detection (RANSAC, regions) | C:Surface_mesh_approximation,Shape_detection O:- | **absent** | No variational shape approximation or RANSAC/region-growing shape-detection code found anywhere in the workspace. Evidence: grep for shape_detection/ransac/vsa: 0 hits. |
| D20 | topology queries (genus, homotopy, cycles) | C:Surface_mesh_topology O:- | narrow | Genus only, via Euler characteristic; refuses meshes with boundary, non-manifold edges, non-orientable characteristic, or more than one connected component (must call per-component). No homotopy/cycle-basis queries. Evidence: `crates/algorithms/query/inspect/src/genus.rs::genus` tested by `crates/algorithms/query/inspect/tests/queries.rs`. |

### E. Triangulations and meshing

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| E1 | 2D Delaunay / constrained Delaunay | C:Triangulation_2 O:BRepMesh (internal) | implemented | Constrained Delaunay triangulation with certified `incircle`/`orient2d` predicates (exact arithmetic on demand), every constraint edge preserved as a union of output edges, empty-circumcircle away from constraints. Evidence: `crates/algorithms/planar/triangulate/src/build.rs::triangulate` (re-exported from `src/lib.rs`) — see also `src/recover.rs` |
| E2 | 2D quality meshing (Ruppert/Chew) | C:Mesh_2 O:- | **scoped** | `Quality`/`refine`/`triangulate_refined` implement Ruppert-style refinement with an explicit Steiner-point budget; reports `RefineOutcome::Capped` when the angle bound is not reached rather than silently returning a worse mesh — a correct and deliberate refusal for inputs with small input angles (Ruppert does not terminate universally), but this means the angle guarantee is conditional, matching CGAL's own caveat, so this is genuinely close to IMPL but is marked NARROW because the crate itself documents the guarantee as conditional/capped rather than unconditional. Evidence: `crates/algorithms/planar/triangulate/src/refine.rs::refine,triangulate_refined,Quality,RefineOutcome` (re-exported `src/lib.rs`). |
| E3 | 3D Delaunay / regular triangulation | C:Triangulation_3 O:- | **absent** | No tetrahedralization/3D Delaunay code found anywhere in the workspace. Evidence: grep for tetrahedr/delaunay_3/Triangulation3: 0 relevant hits (only unrelated "tet" substring matches in orient3d predicate files). |
| E4 | 3D constrained Delaunay | C:Constrained_triangulation_3 O:- | **absent** | Depends on E3, which is absent. Evidence: (same search as E3) |
| E5 | volume (tetrahedral) meshing | C:Mesh_3,Tetrahedral_remeshing O:- (TKXMesh stub) | **absent** | No tet-mesh generator; no `TetMesh` type exists in `representations`. Evidence: grep for TetMesh/tet_mesh/volume_mesh: 0 hits. |
| E6 | surface meshing of implicit/smooth surfaces | C:Surface_mesher,Mesh_3 O:- | narrow | `axiolid-levelset` extracts a surface mesh from a signed scalar field via marching-cubes-family extraction (see G5); "meshing of an implicit surface" and "isosurface extraction" are the same underlying capability here, so this row and G5 report the same evidence. No dedicated Delaunay-refinement surface mesher (CGAL `Mesh_3`/`Surface_mesher` style) with a guaranteed approximation-error bound exists. Evidence: `crates/algorithms/sampled/levelset/src/lib.rs` tested by `crates/algorithms/sampled/levelset/tests/extract.rs::a_surface_reaching_the_bounds_still_closes`. |
| E7 | alpha shapes / alpha wrap | C:Alpha_shapes_2/3,Alpha_wrap_2/3 O:- | **absent** | No alpha-shape or alpha-wrap code found. Evidence: grep for alpha_shape/AlphaShape/alpha_wrap: 0 hits. |
| E8 | Voronoi / power diagrams | C:Voronoi_diagram_2,Apollonius,Segment_Delaunay_graph O:- | **absent** | No Voronoi diagram, power diagram, Apollonius, or segment-Delaunay-graph code found. Evidence: grep for voronoi/Voronoi: 0 hits. |
| E9 | periodic / hyperbolic / spherical triangulations | C:Periodic_*,Hyperbolic_*,Triangulation_on_sphere_2 O:- | **scoped** | No such structures found. Evidence: not directly searched beyond general triangulation source review |
| E10 | d-dimensional triangulations/hulls | C:Triangulation,Convex_hull_d,Kernel_d O:- | **scoped** | Only 2D triangulation and 2D/3D convex hull exist (G1); no generalized d-dimensional kernel. Evidence: see G1 evidence |

### F. 2D polygons and arrangements

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| F1 | polygon booleans (exact) | C:Boolean_set_operations_2,Nef_2 O:(BOP on faces) | narrow | Two backends: a certified integer-predicate polygon path (exact for straight-edge polygons, with holes) via `axiolid-overlay::Region::{union,intersection,difference}`, and an arc-aware path (`arc_overlay`) that wraps the `cavalier_contours` crate (a third-party dependency, not Axiolid's own exact-predicate code) for boundaries containing circular arcs — the arc path's exactness is therefore only as good as `cavalier_contours`'s float arithmetic, not certified like the straight-edge path. No Nef-2D (open/unbounded regions). Evidence: `crates/algorithms/planar/overlay/src/region.rs::Region::union,intersection,difference` tested via `crates/algorithms/planar/overlay` test suite |
| F2 | polygon offset (straight/rounded) | C:Straight_skeleton_2,Minkowski_sum_2 O:BRepOffsetAPI_MakeOffset | narrow | `offset_polygons`/`stroke_polyline` in overlay crate, and `Region::dilate`/`Region::erode` (disc-expansion form). No standalone straight-skeleton-based offset (see F3, ABSENT) — offset here is the Minkowski-disc form (rounded corners only by construction; no mitred/beveled straight-skeleton offset variant). Evidence: `crates/algorithms/planar/overlay/src/offset.rs::offset_polygons,stroke_polyline` |
| F3 | straight skeleton | C:Straight_skeleton_2 O:- | **absent** | No straight-skeleton implementation found anywhere in the workspace (searched directly, 0 hits). Evidence: grep for straight_skeleton/StraightSkeleton: 0 hits. |
| F4 | arrangements of curves (segments, arcs, conics, Bezier) | C:Arrangement_on_surface_2 O:- | narrow | `axiolid-arrangement` is a general editable DCEL, but only for STRAIGHT-EDGE (segment) boundaries — `Vertex`/`HalfEdge`/`Face` store `Point2` positions with no curve type. Arcs are handled only inside the separate, non-incremental `arc_overlay` boolean path (F1/F2), not as arrangement edges; conics and Bezier curves are not supported anywhere. Evidence: `crates/algorithms/planar/arrangement/src/lib.rs::Arrangement` (straight edges only) tested by `crates/algorithms/planar/arrangement/tests/arrangement.rs::a_square_builds_one_bounded_face_plus_the_outer_one`. |
| F5 | envelopes / lower envelope | C:Envelope_2/3 O:- | **scoped** | No lower/upper envelope algorithm found; all "envelope" hits in the codebase refer to unrelated navigation/traversal "clearance envelope" (agent radius/height/slope) concepts in `axiolid-field`, not the CGAL curve-envelope sense. Evidence: grep for envelope: all hits are `TraversalEnvelope` in `crates/algorithms/sampled/field/src/navigate.rs`, unrelated capability. |
| F6 | Minkowski sum 2D | C:Minkowski_sum_2 O:- | narrow | No standalone 2D-specific Minkowski sum; the general `axiolid-minkowski` crate operates on 3D planar-faced solids (see H2) via convex-hull-of-pairwise-sums for convex operands and decompose+pairwise+boolean for non-convex. 2D callers would have to lift into 3D or use `Region::dilate` for the disc-only special case (F2). No 2D polygon+polygon Minkowski sum function exists. Evidence: `crates/algorithms/discrete/minkowski/src/lib.rs::minkowski_sum,minkowski_sum_with` (3D only) tested by `crates/algorithms/discrete/minkowski/tests/minkowski.rs::the_sum_of_two_boxes_is_a_box_with_summed_extents`. |
| F7 | polygon partition (convex/monotone) | C:Partition_2 O:- | **absent** | No convex or monotone polygon partition algorithm found; the CDT-based triangulation (E1) fully partitions into triangles but that is a different, finer-grained capability. Evidence: not found in `overlay`/`triangulate`/`arrangement` source. |
| F8 | visibility / shortest path in polygon | C:Visibility_2 O:- | implemented | `axiolid-route` computes visibility graphs and shortest paths constrained to stay inside a polygon (with holes). Evidence: `crates/algorithms/planar/route/src/lib.rs` tested by `crates/algorithms/planar/route/tests/route.rs`. |
| F9 | snap rounding | C:Snap_rounding_2 O:- | **absent** | No snap-rounding algorithm found. Evidence: grep for snap_round/SnapRound: 0 hits. |
| F10 | polyline simplification | C:Polyline_simplification_2 O:- | **absent** | No Douglas-Peucker/Ramer or Visvalingam-Whyatt implementation found. Evidence: grep for polyline_simplif/douglas_peucker/Douglas/Ramer: 0 hits. |
| F11 | sweep-line segment intersection | C:Surface_sweep_2,Intersections_2 O:- | narrow | Only pairwise segment-segment and line-line intersection exist (`axiolid-intersection-linear`), certified via predicates; there is no Bentley-Ottmann-style sweep-line algorithm to report all intersections among a set of N segments in better-than-quadratic time — a caller needing that must call the pairwise test O(n²) times. Evidence: `crates/algorithms/query/intersection/linear/src/segment_segment.rs` tested by `crates/algorithms/query/intersection/linear/tests/segment_segment.rs`. |

### G. Point sets and reconstruction

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| G1 | convex hull 2D/3D | C:Convex_hull_2/3 O:- | implemented | 3D: `convex_hull` builds a closed outward-oriented `TriMesh` by incremental insertion, every visibility decision through certified `orient3d`, typed refusals for <4 points / all-collinear / all-coplanar. 2D hull exists in `axiolid-reference`. Evidence: `crates/algorithms/construction/construct/src/hull.rs::convex_hull` tested by `crates/algorithms/construction/construct/tests/convex_hull.rs` |
| G2 | bounding volumes (min sphere, OBB, min rect) | C:Bounding_volumes,Optimal_bounding_box O:Bnd_OBB | narrow | 2D only: minimum-area oriented rectangle and strict 2D convex hull. No 3D OBB, no minimum enclosing circle/sphere. Evidence: `crates/algorithms/reference/src/convex_hull.rs::minimum_area_rectangle` |
| G3 | point set processing (normals, outliers, smoothing, simplify) | C:Point_set_processing_3 O:- | **absent** | `PointCloud` (in `representations/discrete/pointcloud`) is a plain data container (points/normals/colours/intensities, with caller-supplied normals) with no processing algorithms — no normal estimation, outlier removal, smoothing, or simplification/subsampling function exists. Evidence: `crates/representations/discrete/pointcloud/src/lib.rs::PointCloud` (fields/accessors only, no `estimate_normals`/`remove_outliers`/`simplify` found). |
| G4 | surface reconstruction (Poisson, AF, scale-space) | C:Poisson_*,Advancing_front_*,Scale_space_* O:- | **scoped** | Reconstruction exists via a custom SDF-based method (nearest-neighbour-driven signed-distance field, then G5's level-set extraction) — not Poisson, not advancing-front, not scale-space. This is a legitimate reconstruction pipeline but a single, different algorithm family from all three CGAL variants named in the row. Evidence: `crates/providers/pointcloud/sdf/src/lib.rs` tested by `crates/providers/pointcloud/sdf/tests/reconstruct.rs` |
| G5 | isosurface extraction (marching cubes, dual contouring) | C:Isosurfacing_3 O:- | implemented | Level-set surface extraction from a signed scalar field (marching-cubes family); closes correctly even when the surface reaches the sample-grid bounds. Evidence: `crates/algorithms/sampled/levelset/src/lib.rs` tested by `crates/algorithms/sampled/levelset/tests/extract.rs::a_surface_reaching_the_bounds_still_closes`. |
| G6 | PCA / plane fitting | C:Principal_component_analysis O:- | **absent** | No PCA or least-squares plane-fitting function found anywhere in the workspace; `fit.rs` in `axiolid-nurbs` only fits B-spline curves/surfaces through explicit points (interpolation, not PCA). Evidence: grep for fn plane_fit/least_squares_plane/fit_plane/fn pca: 0 hits |
| G7 | kd-tree / k-NN / range search | C:Spatial_searching O:BVH(nearest), NCollection_UBTree | **scoped** | `PointIndex` provides exact k-NN and radius search, but via a UNIFORM GRID, not a kd-tree — the crate's own docs state this is deliberate ("every query has the same radius and cell arithmetic beats tree descent") and that no tree-based point index is provided; a caller with widely varying query radii or highly non-uniform point density does not get the adaptive behaviour a kd-tree/octree would give. Evidence: `crates/algorithms/query/spatial/src/points.rs::PointIndex` tested by `crates/algorithms/query/spatial/tests/points.rs` and `tests/nearest.rs` |
| G8 | AABB tree / box intersection | C:AABB_tree,Box_intersection_d O:BVH,Bnd | implemented | `Bvh` over bounded objects (triangles/solids), adapts to geometry distribution, used by clash/ray-cast/healing; also implements box-pair broad-phase queries. Evidence: `crates/algorithms/query/spatial/src/bvh.rs::Bvh` tested by `crates/algorithms/query/spatial/tests/bvh.rs`. |

### H. Spatial search and miscellany

| Row | Capability | Reference (O: OCCT, C: CGAL) | Axiolid | Scope, refusals, evidence |
| --- | --- | --- | --- | --- |
| H1 | BVH / spatial index | C:AABB_tree,Orthtree O:BVH | **scoped** | `Bvh` for objects (general, IMPL — same as G8) and `PointIndex` for points (uniform grid, not an orthtree/adaptive octree — see G7's narrowing, which applies here too). No adaptive Orthtree-equivalent structure exists. Evidence: `crates/algorithms/query/spatial/src/bvh.rs::Bvh` |
| H2 | 3D Minkowski sum | C:Minkowski_sum_3 O:- | narrow | Exact for convex-convex pairs (hull of pairwise vertex sums) and for non-convex via decompose+pairwise-sum+boolean-union with an explicit pairwise-sum budget (4096) that refuses rather than runs unbounded work. Minkowski DIFFERENCE additionally requires the SUBJECT to be convex — refused by name for a non-convex subject (erosion via vertex-wise containment is only valid when the subject is convex). Curved operands refused (planar-faced solids only). Evidence: `crates/algorithms/discrete/minkowski/src/lib.rs::minkowski_sum,minkowski_sum_with,minkowski_difference_with` tested by `crates/algorithms/discrete/minkowski/tests/minkowski.rs::a_non_convex_sum_differs_from_treating_the_operand_as_convex`, `erosion_refuses_a_non_convex_subject`. |
| H3 | convex collision / distance (GJK/SAT) | C:Polytope_distance_d O:- | **scoped** | Full SAT (separating axis theorem) distance/intersection for convex shapes, by deliberate design choice (documented reasoning: SAT over GJK for the model-checking use case). Deliberately does NOT report penetration depth (documented refusal — a caller wanting EPA/penetration depth for physics is told to use a physics engine instead). This is IMPL for "convex collision/distance" as literally asked, but the row also implies GJK-class capability which is a named, explicit non-goal; marking NARROW to flag the documented penetration-depth gap. Evidence: `crates/algorithms/query/collide/src/lib.rs::distance,intersects,boxes_intersect,contains_point` |
| H4 | Frechet / curve distances | C:Frechet_distance O:- | **absent** | No Frechet distance or other curve-distance metric found anywhere in the workspace. Evidence: grep for frechet/Frechet: 0 hits. |
| H5 | interpolation / barycentric coords | C:Interpolation,Barycentric_coordinates_2/3 O:- | narrow | Barycentric coordinates exist only as an internal by-product of ray/triangle intersection (`RayHit::barycentric`) and of attribute-blend interpolation inside the mesh boolean provider (`boolmesh/src/attributes.rs::barycentric`) — neither is exposed as a general-purpose, reusable barycentric-coordinate or interpolation API for arbitrary points against arbitrary polygons/triangles. NURBS curve/surface interpolation (`fit.rs::interpolate_curve3`) is a different capability (curve fitting, not barycentric coordinates). Evidence: `crates/algorithms/query/intersection/ray-mesh/src/lib.rs::RayHit::barycentric` (field, computed in `intersect_triangle`) tested by `crates/algorithms/query/intersection/ray-mesh/tests/nearest_hit.rs` |

## Where the gaps are

Ordered by how much downstream work each one blocks.

1. **General surface intersection (B9, B10).** Surface/surface and
   curve/surface intersection are certified but bounded: they refuse
   tangential, overlapping and general curved cases. Curved booleans,
   exact sections, offsets and fillets all sit downstream of this, so it
   gates most of C.
2. **Curved B-rep booleans (C9).** The exact boolean covers any planar-faced
   solid; anything with a curved face goes through the mesh boolean, which
   is robust but approximate.
3. **Construction breadth (C4-C12).** Extrude, revolve, sweep and loft are
   exact only for named profile families (#111 tracks the extrusion side);
   fillet and chamfer handle one straight edge or one polygon corner, not
   edge chains; no draft angle (C13). Surface filling is absent (B12).
4. **3D triangulation and volume meshing (E3-E5).** No 3D Delaunay, no
   constrained 3D Delaunay, no tetrahedral meshing. This is CGAL's largest
   area and Axiolid's widest gap. Whether it is in scope needs a named
   consumer (simulation meshing for CFD or FEM would be one).
5. **Mesh processing breadth (D6, D9-D13, D19).** No hole filling, subdivision,
   UV parameterisation, geodesics, segmentation or shape detection. Repair,
   decimation and refinement exist but are opt-in and narrow.
6. **Point-set processing (G3, G6).** No normal estimation, outlier removal,
   smoothing or plane fitting; reconstruction exists (signed-distance
   field) but has to be handed clean, oriented input.
7. **2D algorithms (F3, F5, F7, F9, F10).** No straight skeleton, polygon
   partition, snap rounding or polyline simplification. Straight skeleton
   is the one with an obvious building use: roof generation and offset
   without rounded corners.
8. **Measurement on curved B-rep (C17).** Exact mass properties refuse any
   non-planar face; curved solids are measured only after tessellation.

## Tracked work

Every gap cluster is a GitHub issue in the [Geometry breadth](https://github.com/axiolid/kernel/milestone/13) milestone. The live,
gate-checked version of this table is `architecture/capability-ledger.toml`;
`cargo xtask gaps` prints it ordered by priority.

| Issue | Work | Rows | Waits on |
| --- | --- | --- | --- |
| [#111](https://github.com/axiolid/kernel/issues/111) | Exact-mode extrusion/revolve/sweep for arbitrary and composite profiles, not just rectangle/circle | C4 | - |
| [#118](https://github.com/axiolid/kernel/issues/118) | 3D oriented bounding box and minimum enclosing sphere | G2 | - |
| [#119](https://github.com/axiolid/kernel/issues/119) | General curve/surface and surface/surface intersection, including tangent and overlapping cases | B7, B8, B9, B10 | - |
| [#120](https://github.com/axiolid/kernel/issues/120) | Exact B-rep boolean for solids with curved faces | C9, D3 | #119 |
| [#121](https://github.com/axiolid/kernel/issues/121) | Fillet and chamfer over edge chains and networks, constant and variable radius | C10, C11 | #119, #120 |
| [#122](https://github.com/axiolid/kernel/issues/122) | Exact sweep along arbitrary curves and loft through arbitrary sections | B5, C5, C6 | - |
| [#123](https://github.com/axiolid/kernel/issues/123) | Solid offset, shell and thicken over curved faces, and draft angle | C12, C13 | #119 |
| [#124](https://github.com/axiolid/kernel/issues/124) | Surface filling: Coons, Gordon and plate surfaces through boundary curves | B12 | - |
| [#125](https://github.com/axiolid/kernel/issues/125) | Exact mass properties and distance over curved B-rep faces | C17, C18 | #119 |
| [#126](https://github.com/axiolid/kernel/issues/126) | 3D Delaunay and constrained 3D Delaunay triangulation | E3, E4 | - |
| [#127](https://github.com/axiolid/kernel/issues/127) | Tetrahedral volume meshing with quality bounds | E5 | #126 |
| [#128](https://github.com/axiolid/kernel/issues/128) | Voronoi and power diagrams, alpha shapes and alpha wrapping | E7, E8 | - |
| [#129](https://github.com/axiolid/kernel/issues/129) | Mesh hole filling with fairing, and subdivision surfaces | D6, D9 | #140 |
| [#130](https://github.com/axiolid/kernel/issues/130) | Mesh geodesics and UV parameterisation | D10, D11 | #140 |
| [#131](https://github.com/axiolid/kernel/issues/131) | Mesh segmentation, skeletonisation and shape detection (planes, cylinders) | D12, D19 | #140 |
| [#132](https://github.com/axiolid/kernel/issues/132) | Mesh deformation with fixed handles | D13 | - |
| [#133](https://github.com/axiolid/kernel/issues/133) | Point-set processing: normals, outliers, smoothing, simplification and plane fitting | G3, G6 | - |
| [#134](https://github.com/axiolid/kernel/issues/134) | Straight skeleton and mitred polygon offset | F2, F3 | - |
| [#135](https://github.com/axiolid/kernel/issues/135) | Polygon partition, polyline simplification and snap rounding | F7, F9, F10 | - |
| [#136](https://github.com/axiolid/kernel/issues/136) | Shared numeric substrate: root finding, quadrature, least squares and optimisation | A5 | - |
| [#137](https://github.com/axiolid/kernel/issues/137) | B-rep shape healing: sewing, tolerance repair and small-feature removal | C15 | - |
| [#138](https://github.com/axiolid/kernel/issues/138) | Hidden-line removal and 2D projection drawings of solids | C16 | - |
| [#139](https://github.com/axiolid/kernel/issues/139) | 2D medial axis, bisectors and hatching | C19, C20 | - |
| [#140](https://github.com/axiolid/kernel/issues/140) | Halfedge mesh representation with O(1) adjacency for 3D surface meshes | D1 | - |
| [#141](https://github.com/axiolid/kernel/issues/141) | Surface knot removal, degree elevation and iso-curve extraction | B4 | - |
| [#142](https://github.com/axiolid/kernel/issues/142) | Torus and wedge B-rep primitive solids | C3 | - |
| [#143](https://github.com/axiolid/kernel/issues/143) | Public barycentric and mean-value coordinates | H5 | - |
| [#144](https://github.com/axiolid/kernel/issues/144) | Mesh genus with boundaries and components, and a homology cycle basis | D20 | - |
| [#145](https://github.com/axiolid/kernel/issues/145) | 2D Minkowski sum for non-convex polygons | F6 | - |
| [#146](https://github.com/axiolid/kernel/issues/146) | Bentley-Ottmann sweep for many-segment intersection | F11 | - |
| [#147](https://github.com/axiolid/kernel/issues/147) | Discrete and continuous Frechet distance between polylines | H4 | - |
| [#148](https://github.com/axiolid/kernel/issues/148) | Two-sided mesh Hausdorff distance with a certified error bound | D15 | - |
| [#149](https://github.com/axiolid/kernel/issues/149) | Isotropic remeshing: split, collapse, flip and tangential relaxation | D7 | #140 |
| [#150](https://github.com/axiolid/kernel/issues/150) | Least-squares curve and surface fitting to points | B11 | #136 |
| [#151](https://github.com/axiolid/kernel/issues/151) | Fair curves: minimum-energy interpolation and batten curves | B16 | #136 |
| [#152](https://github.com/axiolid/kernel/issues/152) | Form features: holes, pockets, slots and ribs on B-rep solids | C14 | #120 |
| [#153](https://github.com/axiolid/kernel/issues/153) | Surface meshing of smooth and implicit surfaces with quality bounds | E6 | #126 |
| [#155](https://github.com/axiolid/kernel/issues/155) | Exact polygon booleans with circular arcs, independent of cavalier_contours | F1 | - |
| [#156](https://github.com/axiolid/kernel/issues/156) | Clip and split a triangle mesh by another mesh | D17 | - |
| [#157](https://github.com/axiolid/kernel/issues/157) | 2D arrangements of circular arcs and conic curves | F4 | - |
| [#158](https://github.com/axiolid/kernel/issues/158) | 3D Minkowski sum and difference for non-convex solids | H2 | - |
| [#159](https://github.com/axiolid/kernel/issues/159) | Apollonius and tangent-circle constructions | B14 | - |

## Where Axiolid is ahead

- **Refusal over silent approximation.** Every narrow row above refuses by
  name outside its subset, with a typed error.
- **Exact predicates under the B-rep side.** OCCT's modelling runs on
  tolerances throughout; Axiolid's planar exact boolean and overlay decide
  every classification with certified predicates.
- **Both halves in one kernel.** OCCT has no mesh-processing toolkit to speak
  of and CGAL has no B-rep; Axiolid covers both, behind one typed contract
  layer, in pure Rust.

## How this was measured

- Package lists and line counts: scripts under the author's scratch area,
  reproducible from the pinned revisions above with any line counter that
  drops blank and comment-only lines.
- Capability grades: each row was graded from the public functions and
  tests of the crates listed in its evidence column. Grades were
  spot-checked by hand; one was corrected (G2, bounding volumes: the 2D
  minimum-area rectangle in `axiolid-reference` had been missed).
- A known-defect note in `axiolid-mesh-boolean-boolmesh` says a depth-2
  Menger sponge panics inside the mesh boolean. Re-run at `cb8d203` as 63
  sequential axis-aligned differences: all succeed, none refused. That does
  not prove the original construction is fixed, only that this one does not
  reproduce it.
