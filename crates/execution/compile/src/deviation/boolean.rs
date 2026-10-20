//! The deviation of a boolean's mesh, measured against the exact result
//! (#235).
//!
//! # Measured, not derived
//!
//! The mesh boolean cuts operand meshes that each lie within their own
//! bound of the exact operands, but those bounds are one-sided (exact
//! surface to mesh): a mesh may carry extra area, and an operand mesh can
//! reach outside the exact operand by any amount without breaking its
//! bound. Where it does, the mesh boolean can remove a strip of the other
//! operand's surface that the exact boolean keeps. Even with two-sided
//! operand bounds the cut curve moves by `(e_a + e_b) / sin(angle)`
//! between the operands' surfaces, unbounded as they meet tangentially.
//! So no bound is derived from the operands' bounds here.
//!
//! Instead, where [`ReferenceExactCompiler`] builds the exact result of the
//! same boolean node (differences of placed extrusions, #228, and whatever
//! it compiles exactly later), the distance from that exact B-rep's surface
//! to the boolean's mesh is certified directly, by the branch and bound of
//! [`crate::certify`] over every face of the exact result:
//!
//! - A face's patch is its exact surface over its parameter box, with
//!   `A, B, C` the surface's certified second-derivative bounds
//!   ([`SurfaceBoundOracle`]).
//! - Its domain is the region its loops' pcurves enclose (even-odd). Each
//!   pcurve span is flattened in parameters until its certified chord
//!   bound ([`chord_bound2`]) is within `delta`, and the loops must close
//!   within a gap the result explains (see "Loops closed by a reading"
//!   below). Every point of the exact boundary is then within
//!   `slack = delta + gap` of the polygon, and each polygon point within
//!   `slack` of the boundary, so a cell further than `slack` from every
//!   polygon segment does not meet the boundary. Such a cell is wholly
//!   inside or wholly outside, and the polygon's even-odd count at its
//!   centre says which (the straight-line homotopy from each pcurve span to
//!   its chord stays within `slack` of the chord, so it never crosses the
//!   centre and the winding parity is the exact loop's).
//! - A cell certainly outside is dropped; one within `slack` of the
//!   boundary is kept and bounded whole, which bounds its part inside the
//!   face. Only a cell certainly inside raises the lower bound.
//!
//! The faces' closures cover the exact surface, so the worst cell bounds
//! the distance from every exact surface point to the mesh: the report's
//! one-sided quantity, with nothing derived from the operands. The bound is
//! [`DeviationBound::Certified`].
//!
//! # Relative to the exact compiler's result
//!
//! The bound is to the B-rep [`ReferenceExactCompiler`] returns, and to
//! nothing else. Since #228 that result is the exact boolean of operands
//! moved by at most the caller's tolerance (ADR 0080; the guarantees in
//! `axiolid_brep_boolean`'s crate docs): where faces agree only up to
//! rounding, the general boolean decides them within the tolerance. Near
//! tangency the result of such a decision can differ from the boolean of
//! the unperturbed operands by more than the tolerance, so the bound is
//! not claimed for that boolean, and the contribution's detail says what
//! it is measured against. No perturbation term is added: none is proved.
//!
//! The exact compiler reports whether any such decision fired (#236,
//! [`ReferenceExactCompiler::compile_exact_with_report`]). When none did,
//! its result is the exact boolean of the operands as given, so the bound
//! holds for that boolean too, and the detail says so (operands placed by
//! exact axis matrices, or crossing only transversally). Otherwise the
//! detail keeps "operands within tolerance".
//!
//! One such decision is named on its own (#324): faces of the two operands
//! read as one surface within tolerance (a niche's inner cylinder a unit
//! in the last place off a curved wall's, a cap a rounding residue off a
//! slab's face). The furthest that reading moved an operand is a
//! [`DeviationBound::Certified`] contribution of [`DeviationPath::Boolean`]
//! under the detail [`COINCIDENT_OPERANDS`], as a snap's move is: it says
//! how far the operands the exact result is the boolean of lie from the
//! given ones, not how far the mesh lies from that boolean.
//!
//! # Warped mesh operands (#266)
//!
//! The exact result is compiled with
//! `ReferenceExactCompiler::compile_flattened_with_report`, which also
//! takes mesh operands: an authored polygon or triangle mesh, and a B-rep
//! whose faces are all planes with straight edges or warped faces without
//! a surface, each as the exact polyhedron of the triangles the mesh
//! compiler makes of it (`crate::exact`'s `polyhedral` module), and a
//! union or intersection when an operand is such a mesh. A warped face's
//! triangles lie within its slab width `w` of every reading of the face,
//! both ways (`crate::planar`), so that exact result is the exact boolean
//! of operands perturbed by at most the largest such `w`: ADR 0080's
//! contract, with the perturbation reported. The contribution is the
//! measured distance `d` to it plus `w` (rounded up), under the detail
//! [`FLATTENED`]; the search settles within the budget less `w`, so a
//! boolean whose slab widths fit the budget reports within it. An
//! unwarped mesh (`w = 0`) is the exact polyhedron it describes and is
//! measured like any exact operand. A flattened boolean the exact
//! compiler refuses is unbounded by the refusal's name, like any other.
//!
//! # Which booleans are measured
//!
//! Only those whose result is emitted ([`emitted_booleans`]): a boolean
//! that is an operand of another (a wall's first opening, under its
//! second) is replaced in the report by the outer boolean's measurement,
//! so measuring it too would only cost an exact compilation and a search.
//!
//! # When to stop
//!
//! The search stops as soon as its worst cell is within the requested
//! budget, rather than tightening to within 10% of the sampled maximum: a
//! boolean's report answers whether the mesh is within the budget, and
//! the bound is sound wherever the search stops. Above the budget it
//! tightens as before, so a miss is reported as close to the truth as the
//! search gets.
//!
//! An exact result whose faces are all planes is the exception (#314): its
//! search runs to within 10% of the sampled maximum or to the rounding of
//! the mesh ([`crate::certify`]'s `certify_affine`), so a mesh that is the
//! exact result reports a bound near rounding instead of just under the
//! request. Once a real deviation is shown (far above rounding), by the
//! samples or by the slab width a flattened warped operand adds to the
//! report (#266), the budget settles it again, so such a boolean costs
//! what it did. That needs the boundary cells to converge too: on a plane, a
//! cell the trim crosses is bounded over its pieces in the face
//! ([`Domain::clip`]), and the slack is the chord bound the flattening
//! reached (`0` for a straight trim) rather than the `delta` it allowed.
//!
//! # Loops closed by a reading
//!
//! A loop's consecutive pcurve spans repeat the point they share, and the
//! measured gap between the repeats goes into the slack, so the bound is
//! sound for any gap. Which gaps are admitted only says which loops the
//! result explains (#334). A loop that closes exactly is admitted within
//! the rounding of its parameters (`1e-9` of their extent), as before, and
//! keeps its slack, so #314's tight planar bound is unchanged. Where the
//! exact boolean read features as one within tolerance (merged points,
//! coincident faces), the pcurves of the result's faces close only within
//! those readings, so a gap is also admitted up to the largest linear
//! reading the boolean reported
//! ([`BooleanReport::linear`],
//! over the whole boolean: the report is not kept per face), never beyond
//! the caller's linear tolerance, and turned into parameters by the
//! surface's largest first derivative over the face, so the admitted gap
//! spans no more than that reading on the surface. A wider gap, or one with no
//! reading to explain it, stays refused by name.
//!
//! # Unbounded, by name
//!
//! A boolean the exact compiler refuses (a union or intersection of placed
//! extrusions, an operand that is neither a placed extrusion nor a mesh it
//! flattens, a configuration the general boolean cannot build) carries the
//! refusal's name. A face whose
//! surface family has no derivative bounds, whose pcurve family has no
//! certified chord bound, or whose loops do not close in its parameters
//! within the gap admitted above (a periodic face trimmed by two separate
//! circles) is unbounded by name.
//!
//! The exact result is compiled at the node's own options, so under an
//! instance it is the local solid and the mesh the local mesh; the
//! instance scales the bound like any other.

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::Curve2;
use axiolid_mesh::TriMesh;
use std::collections::HashSet;

use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};
use axiolid_reference::bound::{chord_bound2, continuity_breaks2, SurfaceBoundOracle};
use axiolid_surface::Surface;

use axiolid_brep_boolean::{BooleanReport, ToleranceDecisionKind};

use super::{Deviation, DeviationBound, DeviationPath, COINCIDENT_OPERANDS};
use crate::certify::{certify_affine, certify_within, Cell, Coverage, Patch};
use crate::ReferenceExactCompiler;

/// The detail of a boolean measured against the exact compiler's result:
/// the exact boolean of the operands moved by at most the tolerance (#228).
const MEASURED: &str = "measured against the exact compiler's result, operands within tolerance";
/// The detail of a boolean measured against the exact compiler's result
/// when no decision within tolerance fired: the exact boolean of the
/// operands as given (#236).
const MEASURED_EXACT: &str =
    "measured against the exact compiler's result, the exact boolean of the given operands";
/// The detail of a boolean with a warped mesh operand (#266), measured
/// against the exact boolean of its operands flattened to their triangle
/// meshes, the bound widened by the largest slab width they were
/// flattened by.
const FLATTENED: &str =
    "measured against the exact boolean of operands flattened within their slab width, plus it";
/// The detail of a boolean with no exact result to measure against.
const NO_EXACT: &str = "no exact result";

/// Pcurve segments one face's domain may be flattened into.
const MAX_SEGMENTS: usize = 1 << 16;
/// Halvings of one pcurve span before its chord bound must fit.
const MAX_DEPTH: u32 = 40;
/// Initial cells per face along one parameter direction.
const MAX_GRID: usize = 32;

/// The deviation of boolean node `id`, compiled to `mesh` at `options`.
pub(crate) fn of_boolean(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    mesh: &TriMesh,
) -> Deviation {
    let (exact, report, flattened) =
        match ReferenceExactCompiler::new().compile_flattened_with_report(graph, id, options) {
            Ok(compiled) => compiled,
            Err(error) => {
                return Deviation::one(
                    DeviationPath::Boolean,
                    NO_EXACT,
                    DeviationBound::Unbounded(refusal(&error)),
                )
            }
        };
    // A warped operand was flattened to its triangles (#266): the result
    // is the exact boolean of operands moved by up to `flattened`. No
    // decision within tolerance otherwise: the exact compiler's result is
    // the boolean of the operands as given, and the bound holds for it
    // (#236).
    let detail = if flattened > 0.0 {
        FLATTENED
    } else if report.is_exact() {
        MEASURED_EXACT
    } else {
        MEASURED
    };
    let budget = crate::compiler::chord_error(options);
    // The flattening takes its share of the budget first (#266): the
    // search settles once within what is left, so a boolean whose
    // operands' slab widths fit the budget still reports within it. One
    // that does not fit cannot meet it whatever the search finds, so it
    // settles within the budget too rather than tightening further.
    // The slab width is a deviation the report already carries: like a
    // sampled one, it lets the budget settle an all-planar search (#314).
    let reading = closing_reading(&report, options);
    let measured = if flattened < budget {
        let left = budget - flattened;
        measure(&exact, mesh, left, left, flattened, reading)
    } else {
        measure(&exact, mesh, budget, budget, flattened, reading)
    };
    let bound = match measured {
        Ok(value) => DeviationBound::Certified(widened(value, flattened)),
        Err(reason) => DeviationBound::Unbounded(reason),
    };
    let mut deviation = Deviation::one(DeviationPath::Boolean, detail, bound);
    // Faces of the operands read as one surface within tolerance (#324):
    // the exact result is the boolean of operands moved by up to that much
    // onto each other, and the move is named, as a snap's is.
    if let Some(moved) = report
        .decisions()
        .iter()
        .find(|d| d.kind == ToleranceDecisionKind::CoincidentSupports)
        .map(|d| d.linear)
    {
        deviation.add(
            DeviationPath::Boolean,
            COINCIDENT_OPERANDS,
            DeviationBound::Certified(moved),
        );
    }
    deviation
}

/// The gap within which an exact face's loops may close beyond rounding
/// (#334): the largest linear reading the boolean took, never beyond the
/// caller's linear tolerance (see the module's "Loops closed by a
/// reading").
fn closing_reading(report: &BooleanReport, options: &ExecutionOptions) -> Scalar {
    report.linear().min(options.tolerance().linear())
}

/// The measured distance to the exact boolean of the flattened operands,
/// widened by how far they were flattened (#266): their sum, rounded up
/// when the `f64` sum rounded down (Knuth's two-sum gives the error
/// exactly), so never below the real sum.
fn widened(measured: Scalar, flattened: Scalar) -> Scalar {
    let sum = measured + flattened;
    let back = sum - flattened;
    let error = (flattened - (sum - back)) + (measured - back);
    if error > 0.0 {
        sum.next_up()
    } else {
        sum
    }
}

/// The boolean nodes whose meshes `root` emits: reached from `root`
/// through instances and collections, not through another boolean's
/// operands.
pub(crate) fn emitted_booleans(graph: &GeometryGraph, root: NodeId) -> HashSet<NodeId> {
    let mut emitted = HashSet::new();
    let mut seen = HashSet::new();
    let mut stack = vec![root];
    while let Some(id) = stack.pop() {
        if !seen.insert(id) {
            continue;
        }
        match graph.get(id) {
            Some(GeometryNode::Instance(instance)) => stack.push(instance.source),
            Some(GeometryNode::InstanceAtStation(placed)) => stack.push(placed.source),
            Some(GeometryNode::Collection(members)) => stack.extend(members.iter().copied()),
            Some(GeometryNode::SolidOperation(SolidOperation::Boolean { .. })) => {
                emitted.insert(id);
            }
            _ => {}
        }
    }
    emitted
}

/// The exact compiler's refusal, named.
fn refusal(error: &GeomError) -> &'static str {
    match error {
        GeomError::UnsupportedInput { input, .. } => input,
        GeomError::Degenerate(_) => "boolean: the exact result is degenerate",
        GeomError::BudgetExceeded { resource } => resource,
        _ => "boolean: the exact compiler cannot build the result",
    }
}

/// A certified bound on the distance from every point of `exact`'s faces
/// to `mesh`, or why there is none; the search settles once it is within
/// `enough` (see the module's "When to stop"), an all-planar one only once
/// a deviation is shown, sampled or the `shown` one the report adds. A
/// face's loops may close within `reading`, the largest linear reading the
/// boolean took (see the module's "Loops closed by a reading").
fn measure(
    exact: &ExactBRep,
    mesh: &TriMesh,
    target: Scalar,
    enough: Scalar,
    shown: Scalar,
    reading: Scalar,
) -> Result<Scalar, &'static str> {
    let topology = exact.topology();
    let mut faces = Vec::with_capacity(topology.faces().len());
    for face in topology.faces() {
        let surface = face
            .surface
            .and_then(|id| exact.surfaces().get(id.index()))
            .ok_or("boolean: an exact face has no surface")?;
        let oracle = SurfaceBoundOracle::new(surface)
            .ok_or("boolean: an exact face's surface family has no derivative bounds")?;
        let domain = Domain::new(exact, face, surface, target, reading)?;
        faces.push(ExactFace {
            surface,
            oracle,
            domain,
        });
    }
    let mut cells = Vec::new();
    for (k, face) in faces.iter().enumerate() {
        face.seed(k, &mut cells)?;
    }
    let patches: Vec<&dyn Patch> = faces.iter().map(|f| f as &dyn Patch).collect();
    // All faces planar (#314): the search runs to the measured distance.
    let bound = if faces.iter().all(|f| f.domain.affine) {
        certify_affine(&patches, cells, mesh, enough, shown)
    } else {
        certify_within(&patches, cells, mesh, target, enough)
    };
    bound.ok_or("boolean: an exact face could not be bounded against the mesh")
}

/// One face of the exact result: its surface over its trimmed domain.
struct ExactFace<'a> {
    surface: &'a Surface,
    oracle: SurfaceBoundOracle<'a>,
    domain: Domain,
}

impl ExactFace<'_> {
    /// Cover the domain's box with cells of roughly square image.
    fn seed(&self, patch: usize, cells: &mut Vec<Cell>) -> Result<(), &'static str> {
        let (lo, hi) = (self.domain.lo, self.domain.hi);
        let bounds = self
            .oracle
            .bounds((lo.x, hi.x), (lo.y, hi.y))
            .ok_or("boolean: an exact face's surface cannot be bounded over its domain")?;
        let eu = bounds.du * (hi.x - lo.x);
        let ev = bounds.dv * (hi.y - lo.y);
        let side = eu.min(ev).max(eu.max(ev) / MAX_GRID as Scalar);
        let count = |extent: Scalar| -> usize {
            if side > 0.0 && extent.is_finite() {
                ((extent / side).ceil() as usize).clamp(1, MAX_GRID)
            } else {
                1
            }
        };
        let (nu, nv) = (count(eu), count(ev));
        for i in 0..nu {
            for j in 0..nv {
                let at = |a: Scalar, b: Scalar, k: usize, n: usize| {
                    a + (b - a) * k as Scalar / n as Scalar
                };
                cells.push(Cell {
                    patch,
                    x: (at(lo.x, hi.x, i, nu), at(lo.x, hi.x, i + 1, nu)),
                    y: (at(lo.y, hi.y, j, nv), at(lo.y, hi.y, j + 1, nv)),
                });
            }
        }
        Ok(())
    }
}

impl Patch for ExactFace<'_> {
    fn jet(&self, x: Scalar, y: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        let point = axiolid_reference::surface::evaluate(self.surface, x, y).ok()?;
        let (du, dv) = axiolid_reference::surface::partials(self.surface, x, y).ok()?;
        (du.is_finite() && dv.is_finite()).then_some((point, du, dv))
    }

    fn second(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        let bounds = self.oracle.bounds(x, y)?;
        // Across a knot line of full multiplicity the second partials do
        // not bound the Taylor remainder.
        bounds
            .smooth
            .then_some([bounds.duu, bounds.duv, bounds.dvv])
    }

    fn coverage(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Coverage {
        self.domain.coverage(x, y)
    }

    fn clip(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Option<Vec<Vec<Point2>>> {
        self.domain.clip(x, y)
    }
}

/// A face's domain in its surface parameters: the polygon of its loops'
/// flattened pcurves, every exact boundary point within `slack` of it.
struct Domain {
    segments: Vec<[Point2; 2]>,
    slack: Scalar,
    lo: Point2,
    hi: Point2,
    /// The surface is affine over the domain's box (a plane, #314): its
    /// boundary cells are clipped ([`Self::clip`]), with each segment's
    /// unit normal into the domain in `inward`.
    affine: bool,
    inward: Vec<Vec2>,
    /// Rounding of the parameters' coordinates.
    rounding: Scalar,
}

impl Domain {
    fn new(
        exact: &ExactBRep,
        face: &axiolid_topology::Face<axiolid_brep::SurfaceId>,
        surface: &Surface,
        target: Scalar,
        reading: Scalar,
    ) -> Result<Self, &'static str> {
        let topology = exact.topology();
        let mut rings: Vec<Vec<Point2>> = Vec::new();
        // First pass: a box of the pcurves' ends and middles, only to size
        // the flattening budget.
        let mut uses = Vec::new();
        for bound in &face.bounds {
            let wire = topology
                .loops()
                .get(bound.loop_id.index())
                .ok_or("boolean: an exact face has a dangling loop")?;
            let mut spans = Vec::with_capacity(wire.edges.len());
            for (k, edge_use) in wire.edges.iter().enumerate() {
                let curve = edge_use
                    .pcurve
                    .and_then(|id| exact.curves2().get(id.index()))
                    .ok_or("boolean: an exact face's trim has no pcurve")?;
                let span = exact
                    .pcurve_interval(bound.loop_id, k)
                    .ok_or("boolean: an exact face's pcurve has no interval")?;
                spans.push((curve, span.start, span.end));
            }
            uses.push(spans);
        }
        let mut lo = Point2::splat(Scalar::INFINITY);
        let mut hi = Point2::splat(Scalar::NEG_INFINITY);
        for spans in &uses {
            for &(curve, a, b) in spans {
                for t in [a, b, 0.5 * (a + b)] {
                    let p = axiolid_reference::evaluate2(curve, t)
                        .map_err(|_| "boolean: an exact face's pcurve cannot be evaluated")?;
                    lo = lo.min(p);
                    hi = hi.max(p);
                }
            }
        }
        if !(lo.is_finite() && hi.is_finite()) {
            return Err("boolean: an exact face has no finite domain");
        }
        // A parameter step moves the surface by at most `g` per unit, so a
        // chord bound of `delta` in parameters is `g delta` on the surface;
        // a hundredth of the budget keeps the boundary cells thin.
        let bounds = SurfaceBoundOracle::new(surface)
            .and_then(|o| o.bounds((lo.x, hi.x), (lo.y, hi.y)))
            .ok_or("boolean: an exact face's surface cannot be bounded over its domain")?;
        let g = bounds.du.max(bounds.dv).max(Scalar::MIN_POSITIVE);
        let delta = 0.01 * target / g;
        let affine = bounds.smooth && [bounds.duu, bounds.duv, bounds.dvv] == [0.0; 3];
        let mut count = 0_usize;
        let mut chord: Scalar = 0.0;
        for spans in &uses {
            let mut ring = Vec::new();
            for &(curve, a, b) in spans {
                flatten(curve, a, b, delta, &mut ring, &mut count, &mut chord)?;
            }
            rings.push(ring);
        }
        // A plane's slack is the chord bounds the flattening reached, not
        // what it allowed (#314): zero for straight trims.
        let delta = if affine { chord } else { delta };
        let scale = (hi - lo).length().max(lo.length()).max(hi.length());
        // A loop closes within the rounding of its parameters, or within
        // the boolean's largest reading (#334): a gap of `p` in parameters
        // spans at most `g p` on the surface, so one within `reading / g`
        // spans no more than the reading that explains it.
        // The slack takes the gap measured, not the one admitted, so a
        // loop that closes exactly keeps its tight slack (#314).
        let closing = (1e-9 * scale.max(1.0)).max(reading / g);
        let mut gap: Scalar = 0.0;
        // Consecutive spans repeat (within rounding) the point they share,
        // and the ring closes from its last point to its first: measure
        // those gaps to refuse a loop that does not close and to fold the
        // rounding into the slack.
        for spans in &uses {
            let ends: Vec<(Point2, Point2)> = spans
                .iter()
                .map(|&(curve, a, b)| {
                    Ok((
                        axiolid_reference::evaluate2(curve, a)
                            .map_err(|_| "boolean: an exact face's pcurve cannot be evaluated")?,
                        axiolid_reference::evaluate2(curve, b)
                            .map_err(|_| "boolean: an exact face's pcurve cannot be evaluated")?,
                    ))
                })
                .collect::<Result<_, &'static str>>()?;
            for k in 0..ends.len() {
                let next = ends[(k + 1) % ends.len()].0;
                gap = gap.max((ends[k].1 - next).length());
            }
        }
        if gap > closing {
            return Err("boolean: an exact face's loops do not close in its parameters");
        }
        Ok(Self::from_rings(&rings, delta, gap, scale, affine))
    }

    /// The domain of `rings`, each exact boundary point within `delta` of
    /// a chord and every gap between consecutive chords within `gap`.
    fn from_rings(
        rings: &[Vec<Point2>],
        delta: Scalar,
        gap: Scalar,
        scale: Scalar,
        affine: bool,
    ) -> Self {
        let mut segments = Vec::new();
        for ring in rings {
            // A plane's boundary cells are clipped, which needs each corner
            // as one vertex (#314): the repeats of consecutive spans'
            // shared points (and the close) are dropped. Each dropped point
            // is within `gap` of the one kept, so the polygon moves by at
            // most `gap`, which the slack holds as it held the gap
            // segments. Otherwise the gaps stay segments of their own.
            let mut kept: Vec<Point2> = Vec::with_capacity(ring.len());
            for &p in ring {
                if !affine || kept.last().is_none_or(|&q| (p - q).length() > gap) {
                    kept.push(p);
                }
            }
            while affine && kept.len() > 1 && (kept[kept.len() - 1] - kept[0]).length() <= gap {
                kept.pop();
            }
            for pair in kept.windows(2) {
                segments.push([pair[0], pair[1]]);
            }
            if let (Some(&first), Some(&last)) = (kept.first(), kept.last()) {
                segments.push([last, first]);
            }
        }
        let slack = (delta + gap) * (1.0 + 1e-9) + 4.0 * Scalar::EPSILON * scale;
        // The domain lies inside its boundary's box, and the boundary
        // within `slack` of the polygon: the polygon's box, grown, holds
        // it.
        let (lo, hi) = rings.iter().flatten().fold(
            (
                Point2::splat(Scalar::INFINITY),
                Point2::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        let inward = if affine {
            (0..segments.len())
                .map(|k| inward_normal(&segments, k))
                .collect()
        } else {
            Vec::new()
        };
        Self {
            segments,
            slack,
            lo: lo - Point2::splat(slack),
            hi: hi + Point2::splat(slack),
            affine,
            inward,
            rounding: 16.0 * Scalar::EPSILON * scale,
        }
    }

    /// The segments within `slack` of the box `[lo, hi]`.
    fn near(&self, lo: Point2, hi: Point2) -> impl Iterator<Item = usize> + '_ {
        self.segments
            .iter()
            .enumerate()
            .filter(move |(_, &[a, b])| {
                let (slo, shi) = (a.min(b), a.max(b));
                let apart = (lo - shi).max(slo - hi);
                if apart.x > self.slack || apart.y > self.slack {
                    return false;
                }
                segment_box_distance(a, b, lo, hi) <= self.slack
            })
            .map(|(k, _)| k)
    }

    /// The cell's part within `slack` of the domain, as convex pieces
    /// (#314), when one segment or two at a shared corner come within
    /// `slack` of it; `None` (bound the cell whole) otherwise.
    ///
    /// No other segment comes within `slack` of the cell, so a near
    /// segment's ends lie further than `slack` from it, and its line meets
    /// the cell (grown by `slack`) only along the segment: the distance to
    /// the cell is convex along the line, and above `slack` at the ends.
    /// So in the cell the domain is one half-plane `H` (one segment), or
    /// near a corner `H1 ∩ H2` (the domain turns inward there) or
    /// `H1 ∪ H2` (it turns outward). The points within `slack` of it lie in
    /// the same of the half-planes moved out by `slack`: a mitred offset,
    /// which holds the rounded one. Each piece is the cell clipped by
    /// them, and a union gives one piece per half-plane.
    fn clip(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Option<Vec<Vec<Point2>>> {
        if !self.affine {
            return None;
        }
        let (lo, hi) = (Point2::new(x.0, y.0), Point2::new(x.1, y.1));
        let mut near = self.near(lo, hi);
        let first = near.next()?;
        let second = near.next();
        if near.next().is_some() {
            return None;
        }
        let cell = vec![lo, Point2::new(hi.x, lo.y), hi, Point2::new(lo.x, hi.y)];
        // Moved out by the slack, and by the rounding of the clip itself.
        let shift = self.slack + self.rounding;
        let half = |polygon: &[Point2], k: usize| {
            clip_half(polygon, self.segments[k][0], self.inward[k], shift)
        };
        let pieces = match second {
            None => vec![half(&cell, first)],
            Some(second) => {
                // The corner, from `a` through `v` to `b`.
                let (i, j) = if self.segments[first][1] == self.segments[second][0] {
                    (first, second)
                } else if self.segments[second][1] == self.segments[first][0] {
                    (second, first)
                } else {
                    return None;
                };
                let (v, b) = (self.segments[j][0], self.segments[j][1]);
                // Inward only when clearly so: the union holds the
                // intersection, so a doubtful corner takes it.
                let turn = (b - v).dot(self.inward[i]);
                if turn > 1e-9 * (b - v).length() {
                    vec![half(&half(&cell, i), j)]
                } else {
                    vec![half(&cell, i), half(&cell, j)]
                }
            }
        };
        Some(pieces.into_iter().filter(|p| !p.is_empty()).collect())
    }

    fn coverage(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Coverage {
        let (lo, hi) = (Point2::new(x.0, y.0), Point2::new(x.1, y.1));
        if self.near(lo, hi).next().is_some() {
            return Coverage::Boundary;
        }
        let centre = Point2::new(0.5 * (x.0 + x.1), 0.5 * (y.0 + y.1));
        if self.inside(centre) {
            Coverage::Inside
        } else {
            Coverage::Outside
        }
    }

    /// Even-odd count of the polygon's crossings of the ray to `+x`.
    fn inside(&self, p: Point2) -> bool {
        let mut inside = false;
        for &[a, b] in &self.segments {
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > p.x {
                    inside = !inside;
                }
            }
        }
        inside
    }
}

/// Points of `curve` from `a` to `b` whose chords are certified within
/// `delta` of it, appended to `out`; `chord` is raised to the largest
/// chord bound used.
fn flatten(
    curve: &Curve2,
    a: Scalar,
    b: Scalar,
    delta: Scalar,
    out: &mut Vec<Point2>,
    count: &mut usize,
    chord: &mut Scalar,
) -> Result<(), &'static str> {
    let (lo, hi) = (a.min(b), a.max(b));
    let mut cuts: Vec<Scalar> = continuity_breaks2(curve, 1)
        .into_iter()
        .filter(|&t| t > lo && t < hi)
        .collect();
    cuts.sort_by(Scalar::total_cmp);
    if a > b {
        cuts.reverse();
    }
    let mut params = vec![a];
    let mut from = a;
    for to in cuts.into_iter().chain(core::iter::once(b)) {
        split(curve, from, to, delta, 0, &mut params, chord)?;
        from = to;
    }
    *count += params.len();
    if *count > MAX_SEGMENTS {
        return Err("boolean: an exact face's trim needs too many segments to certify");
    }
    for t in params {
        out.push(
            axiolid_reference::evaluate2(curve, t)
                .map_err(|_| "boolean: an exact face's pcurve cannot be evaluated")?,
        );
    }
    Ok(())
}

/// Halve `[a, b]` until each piece's chord bound fits, pushing each piece's
/// end parameter in order and raising `chord` to the bounds that fit.
fn split(
    curve: &Curve2,
    a: Scalar,
    b: Scalar,
    delta: Scalar,
    depth: u32,
    params: &mut Vec<Scalar>,
    chord: &mut Scalar,
) -> Result<(), &'static str> {
    let (lo, hi) = (a.min(b), a.max(b));
    if let Some(d) = chord_bound2(curve, lo, hi).filter(|&d| d <= delta) {
        *chord = chord.max(d);
        params.push(b);
        return Ok(());
    }
    if depth >= MAX_DEPTH {
        return Err("boolean: an exact face's pcurve family has no certified chord bound");
    }
    let m = 0.5 * (a + b);
    split(curve, a, m, delta, depth + 1, params, chord)?;
    split(curve, m, b, delta, depth + 1, params, chord)
}

/// The unit normal of segment `k` into the domain: its left normal when
/// the ray along it from the segment's middle crosses the other segments
/// an odd number of times (even-odd, as [`Domain::inside`]), else the
/// right. A segment with no length has none (zero), and clips nothing.
fn inward_normal(segments: &[[Point2; 2]], k: usize) -> Vec2 {
    let [a, b] = segments[k];
    let along = b - a;
    let length = along.length();
    if length.is_nan() || length <= 0.0 {
        return Vec2::ZERO;
    }
    let n = Vec2::new(-along.y, along.x) / length;
    let m = (a + b) * 0.5;
    let mut odd = false;
    for (j, &[p, q]) in segments.iter().enumerate() {
        if j == k {
            continue;
        }
        // `s` across the ray, `t` along it.
        let (sp, sq) = (n.perp_dot(p - m), n.perp_dot(q - m));
        if (sp > 0.0) != (sq > 0.0) {
            let (tp, tq) = (n.dot(p - m), n.dot(q - m));
            if tp + (tq - tp) * (sp / (sp - sq)) > 0.0 {
                odd = !odd;
            }
        }
    }
    if odd {
        n
    } else {
        -n
    }
}

/// The convex `polygon` clipped to `{p : (p - a) . n >= -shift}`
/// (Sutherland-Hodgman); empty when nothing is left.
fn clip_half(polygon: &[Point2], a: Point2, n: Vec2, shift: Scalar) -> Vec<Point2> {
    let level = |p: Point2| (p - a).dot(n) + shift;
    let mut out = Vec::with_capacity(polygon.len() + 1);
    for k in 0..polygon.len() {
        let (p, q) = (polygon[k], polygon[(k + 1) % polygon.len()]);
        let (lp, lq) = (level(p), level(q));
        if lp >= 0.0 {
            out.push(p);
        }
        if (lp >= 0.0) != (lq >= 0.0) {
            out.push(p + (q - p) * (lp / (lp - lq)));
        }
    }
    out
}

/// Distance from segment `ab` to the box `[lo, hi]`: zero when they meet,
/// else the least of the endpoints' distances to the box and the box
/// corners' distances to the segment (for a segment and a convex polygon
/// that miss each other, the closest pair has a vertex of one).
fn segment_box_distance(a: Point2, b: Point2, lo: Point2, hi: Point2) -> Scalar {
    let in_box = |p: Point2| p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y;
    if in_box(a) || in_box(b) {
        return 0.0;
    }
    let corners = [lo, Point2::new(hi.x, lo.y), hi, Point2::new(lo.x, hi.y)];
    for k in 0..4 {
        if segments_cross(a, b, corners[k], corners[(k + 1) % 4]) {
            return 0.0;
        }
    }
    let to_box = |p: Point2| (lo - p).max(p - hi).max(Point2::ZERO).length();
    let to_segment = |p: Point2| {
        let ab = b - a;
        let len = ab.length_squared();
        let t = if len > 0.0 {
            ((p - a).dot(ab) / len).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (p - (a + ab * t)).length()
    };
    corners
        .iter()
        .map(|&c| to_segment(c))
        .fold(to_box(a).min(to_box(b)), Scalar::min)
}

/// Whether closed segments `pq` and `rs` meet (a touch counts).
fn segments_cross(p: Point2, q: Point2, r: Point2, s: Point2) -> bool {
    let orient = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a);
    let (d1, d2) = (orient(r, s, p), orient(r, s, q));
    let (d3, d4) = (orient(p, q, r), orient(p, q, s));
    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
    {
        return true;
    }
    let on = |a: Point2, b: Point2, c: Point2, d: Scalar| {
        d == 0.0
            && c.x >= a.x.min(b.x)
            && c.x <= a.x.max(b.x)
            && c.y >= a.y.min(b.y)
            && c.y <= a.y.max(b.y)
    };
    on(r, s, p, d1) || on(r, s, q, d2) || on(p, q, r, d3) || on(p, q, s, d4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{BooleanOperator, Tolerance, Transform3, Vec3};
    use axiolid_exact_compile_contract::ExactCompiler;
    use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
    use axiolid_mesh_compile_contract::MeshCompiler;
    use axiolid_model::{GeometryGraphBuilder, Instance};
    use axiolid_profile::{CircleProfile, Profile, RectangleProfile, SectionProfile};

    /// Tightened to the sampled maximum (no budget stop), the bound still
    /// covers the exact result's trim curves, where the cut lies and an
    /// I-beam's hole edge is furthest from the mesh.
    #[test]
    fn a_tight_bound_covers_the_cut() {
        let mut g = GeometryGraphBuilder::new();
        let section = g
            .push(GeometryNode::Profile(Profile::Section(SectionProfile::I {
                depth: 0.3,
                width: 0.15,
                web_thickness: 0.0071,
                flange_thickness: 0.0107,
                fillet_radius: Some(0.015),
                flange_edge_radius: None,
                flange_slope: None,
            })))
            .unwrap();
        let beam = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile: section,
                direction: Vec3::Z,
                depth: 1.0,
            }))
            .unwrap();
        let disk = g
            .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
                radius: 0.05,
                thickness: None,
            })))
            .unwrap();
        let hole = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile: disk,
                direction: Vec3::Z,
                depth: 0.2,
            }))
            .unwrap();
        let placed = g
            .push(GeometryNode::Instance(Instance {
                source: hole,
                transform: Transform3::from_translation(Vec3::new(-0.1, 0.02, 0.5))
                    * Transform3::from_rotation_y(core::f64::consts::FRAC_PI_2),
            }))
            .unwrap();
        let cut = g
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: beam,
                right: placed,
                operator: BooleanOperator::Difference,
            }))
            .unwrap();
        let graph = g.finish(vec![cut]).unwrap();
        let budget = 1e-3;
        let options = ExecutionOptions::new(Tolerance::new(budget, 1e-9).unwrap())
            .with_chord_error(budget)
            .unwrap();
        let mesh = crate::ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_reported(&graph, cut, &options)
            .unwrap()
            .mesh;
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, cut, &options)
            .unwrap();
        let bound = measure(&exact, &mesh, budget, 0.0, 0.0, 0.0).unwrap();
        let index = crate::certify::TriangleIndex::new(&mesh).unwrap();
        let topology = exact.topology();
        let mut worst: Scalar = 0.0;
        for face in topology.faces() {
            let surface = &exact.surfaces()[face.surface.unwrap().index()];
            for bound in &face.bounds {
                let wire = &topology.loops()[bound.loop_id.index()];
                for (k, edge_use) in wire.edges.iter().enumerate() {
                    let curve = &exact.curves2()[edge_use.pcurve.unwrap().index()];
                    let span = exact.pcurve_interval(bound.loop_id, k).unwrap();
                    for s in 0..=2000 {
                        let t = span.start + (span.end - span.start) * Scalar::from(s) / 2000.0;
                        let uv = axiolid_reference::evaluate2(curve, t).unwrap();
                        let p = axiolid_reference::surface::evaluate(surface, uv.x, uv.y).unwrap();
                        worst = worst.max(index.nearest(p));
                    }
                }
            }
        }
        assert!(worst <= bound, "{worst} above the bound {bound}");
        assert!(bound <= 1.2 * worst, "{bound} for a sampled {worst}");
    }

    /// A chain's inner booleans are not emitted unless something else
    /// emits them too: here the first cut is also a collection member.
    #[test]
    fn only_emitted_booleans_are_measured() {
        let mut g = GeometryGraphBuilder::new();
        let profile = g
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x: 1.0,
                    y: 1.0,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap();
        let block = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth: 1.0,
            }))
            .unwrap();
        let cut = |g: &mut GeometryGraphBuilder, left| {
            g.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right: block,
                operator: BooleanOperator::Difference,
            }))
            .unwrap()
        };
        let first = cut(&mut g, block);
        let second = cut(&mut g, first);
        let third = cut(&mut g, second);
        let placed = g
            .push(GeometryNode::Instance(Instance {
                source: third,
                transform: Transform3::IDENTITY,
            }))
            .unwrap();
        let both = g
            .push(GeometryNode::Collection(vec![placed, first]))
            .unwrap();
        let graph = g.finish(vec![placed, both]).unwrap();
        assert_eq!(emitted_booleans(&graph, placed), HashSet::from([third]));
        assert_eq!(
            emitted_booleans(&graph, both),
            HashSet::from([third, first])
        );
    }

    fn square(lo: (Scalar, Scalar), hi: (Scalar, Scalar)) -> Vec<Point2> {
        vec![
            Point2::new(lo.0, lo.1),
            Point2::new(hi.0, lo.1),
            Point2::new(hi.0, hi.1),
            Point2::new(lo.0, hi.1),
        ]
    }

    fn to_polygon(p: Point2, domain: &Domain) -> Scalar {
        domain
            .segments
            .iter()
            .map(|&[a, b]| {
                let ab = b - a;
                let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
                (p - (a + ab * t)).length()
            })
            .fold(Scalar::INFINITY, Scalar::min)
    }

    /// In a convex polygon of either winding, up to `tolerance`.
    fn in_piece(p: Point2, piece: &[Point2], tolerance: Scalar) -> bool {
        let n = piece.len();
        let sides: Vec<Scalar> = (0..n)
            .map(|k| {
                let edge = piece[(k + 1) % n] - piece[k];
                edge.perp_dot(p - piece[k]) / edge.length().max(Scalar::MIN_POSITIVE)
            })
            .collect();
        sides.iter().all(|&s| s >= -tolerance) || sides.iter().all(|&s| s <= tolerance)
    }

    /// The clipped pieces of a plane's boundary cells (#314) hold every
    /// point the face may have in the cell (in the polygon, or within the
    /// slack of it), and reach no further than the slack: at straight
    /// edges, at the outward corners of an outer loop and at the inward
    /// corners of a hole and of an L.
    #[test]
    fn clipped_pieces_hold_the_domain_near_its_boundary_and_no_more() {
        let delta = 0.02;
        let l_shape = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 1.5),
            Point2::new(1.5, 1.5),
            Point2::new(1.5, 4.0),
            Point2::new(0.0, 4.0),
        ];
        let cases = [
            vec![
                square((0.0, 0.0), (4.0, 4.0)),
                square((1.0, 1.0), (2.5, 2.0)).into_iter().rev().collect(),
            ],
            vec![l_shape],
        ];
        for rings in cases {
            let domain = Domain::from_rings(&rings, delta, 0.0, 4.0, true);
            let reach = domain.slack + domain.rounding;
            let (mut clipped, mut corners) = (0, 0);
            for offset in [0.0, 0.07, 0.13] {
                for i in 0..20 {
                    for j in 0..20 {
                        let x0 = -0.5 + offset + 0.25 * Scalar::from(i);
                        let y0 = -0.5 + offset * 0.7 + 0.25 * Scalar::from(j);
                        let (x, y) = ((x0, x0 + 0.25), (y0, y0 + 0.25));
                        if domain.coverage(x, y) != Coverage::Boundary {
                            continue;
                        }
                        let Some(pieces) = domain.clip(x, y) else {
                            continue;
                        };
                        clipped += 1;
                        corners += usize::from(
                            domain
                                .near(Point2::new(x.0, y.0), Point2::new(x.1, y.1))
                                .count()
                                == 2,
                        );
                        for s in 0..=24 {
                            for t in 0..=24 {
                                let p = Point2::new(
                                    x.0 + (x.1 - x.0) * Scalar::from(s) / 24.0,
                                    y.0 + (y.1 - y.0) * Scalar::from(t) / 24.0,
                                );
                                let may = domain.inside(p) || to_polygon(p, &domain) <= delta;
                                if may {
                                    assert!(
                                        pieces.iter().any(|piece| in_piece(p, piece, 1e-12)),
                                        "{p:?} of the cell {x:?} x {y:?} is in no piece: {pieces:?}"
                                    );
                                }
                            }
                        }
                        for &v in pieces.iter().flatten() {
                            // A mitred corner reaches `sqrt 2` of the
                            // slack out at a right angle.
                            assert!(
                                domain.inside(v)
                                    || to_polygon(v, &domain) <= 2f64.sqrt() * reach + 1e-12,
                                "{v:?} lies {} outside the domain",
                                to_polygon(v, &domain)
                            );
                        }
                    }
                }
            }
            assert!(
                clipped > 100 && corners > 4,
                "{clipped} cells, {corners} at corners"
            );
        }
    }

    /// A plane's slack is the chord bound its flattening reached (#314):
    /// every point of a curved trim lies within it of the polygon, and a
    /// face with straight trims only has a slack at rounding.
    #[test]
    fn a_plane_slack_holds_its_trims() {
        let mut g = GeometryGraphBuilder::new();
        let rectangle = g
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x: 6.0,
                    y: 0.3,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap();
        let wall = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile: rectangle,
                direction: Vec3::Z,
                depth: 3.0,
            }))
            .unwrap();
        let disk = g
            .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
                radius: 0.4,
                thickness: None,
            })))
            .unwrap();
        let window = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile: disk,
                direction: Vec3::Z,
                depth: 1.0,
            }))
            .unwrap();
        let across = g
            .push(GeometryNode::Instance(Instance {
                source: window,
                transform: Transform3::from_translation(Vec3::new(0.8, 0.5, 1.5))
                    * Transform3::from_rotation_x(core::f64::consts::FRAC_PI_2),
            }))
            .unwrap();
        let cut = g
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: wall,
                right: across,
                operator: BooleanOperator::Difference,
            }))
            .unwrap();
        let graph = g.finish(vec![cut]).unwrap();
        let target = 1e-3;
        let options = ExecutionOptions::new(Tolerance::new(target, 1e-9).unwrap())
            .with_chord_error(target)
            .unwrap();
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, cut, &options)
            .unwrap();
        let topology = exact.topology();
        let (mut curved, mut straight) = (0, 0);
        for face in topology.faces() {
            let surface = &exact.surfaces()[face.surface.unwrap().index()];
            let domain = Domain::new(&exact, face, surface, target, 0.0).unwrap();
            if !domain.affine {
                continue;
            }
            let mut worst: Scalar = 0.0;
            let mut trims = 0;
            for bound in &face.bounds {
                let wire = &topology.loops()[bound.loop_id.index()];
                for (k, edge_use) in wire.edges.iter().enumerate() {
                    let curve = &exact.curves2()[edge_use.pcurve.unwrap().index()];
                    trims += usize::from(!matches!(curve, Curve2::Line(_)));
                    let span = exact.pcurve_interval(bound.loop_id, k).unwrap();
                    for s in 0..=4000 {
                        let t = span.start + (span.end - span.start) * Scalar::from(s) / 4000.0;
                        let uv = axiolid_reference::evaluate2(curve, t).unwrap();
                        worst = worst.max(to_polygon(uv, &domain));
                    }
                }
            }
            assert!(worst <= domain.slack, "{worst:e} beyond {:e}", domain.slack);
            if trims == 0 {
                straight += 1;
                assert!(domain.slack <= 1e-12, "{:e}", domain.slack);
            } else {
                curved += 1;
                assert!(worst > 1e-9, "a flattened circle: {worst:e}");
            }
        }
        assert!(
            curved >= 2 && straight >= 2,
            "{curved} curved, {straight} straight"
        );
    }

    /// A wall `5 x 3 x 0.25` less a door on its bottom edge and a window,
    /// their caps `1e-8` outside and inside its faces (#334): read as one surface within tolerance, the
    /// exact result's loops close only within the readings. They are
    /// admitted within the boolean's largest reading, capped by the
    /// caller's tolerance, and refused without it or with a reading too
    /// small to explain them.
    #[test]
    fn loops_close_within_the_booleans_readings() {
        let e = 1e-8;
        let mut g = GeometryGraphBuilder::new();
        let rectangle = |g: &mut GeometryGraphBuilder, x, y| {
            g.push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x,
                    y,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap()
        };
        let solid = |g: &mut GeometryGraphBuilder, profile, depth, at: Vec3| {
            let solid = g
                .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                    profile,
                    direction: Vec3::Z,
                    depth,
                }))
                .unwrap();
            g.push(GeometryNode::Instance(Instance {
                source: solid,
                transform: Transform3::from_translation(at),
            }))
            .unwrap()
        };
        let wall = rectangle(&mut g, 5.0, 3.0);
        let door = rectangle(&mut g, 1.0, 2.2);
        let window = rectangle(&mut g, 1.0, 1.0);
        let mut cut = solid(&mut g, wall, 0.25, Vec3::ZERO);
        for (profile, at, depth) in [
            (door, Vec3::new(-1.5, -0.4, -e), 0.25 + 2.0 * e),
            (window, Vec3::new(1.0, 0.0, e), 0.25 - 2.0 * e),
        ] {
            let tool = solid(&mut g, profile, depth, at);
            cut = g
                .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                    left: cut,
                    right: tool,
                    operator: BooleanOperator::Difference,
                }))
                .unwrap();
        }
        let graph = g.finish(vec![cut]).unwrap();
        let target = 1e-3;
        let options = ExecutionOptions::new(Tolerance::new(target, 1e-9).unwrap())
            .with_chord_error(target)
            .unwrap();
        let (exact, report) = ReferenceExactCompiler::new()
            .compile_exact_with_report(&graph, cut, &options)
            .unwrap();
        let reading = report.linear();
        assert!(reading >= e && reading <= 10.0 * e, "{report:?}");
        // The reading, capped by the caller's tolerance; none when exact.
        assert_eq!(closing_reading(&report, &options), reading);
        let fine = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).unwrap());
        assert_eq!(closing_reading(&report, &fine), 1e-9);
        assert_eq!(closing_reading(&BooleanReport::default(), &options), 0.0);
        let open = "boolean: an exact face's loops do not close in its parameters";
        let domains = |reading: Scalar| {
            exact
                .topology()
                .faces()
                .iter()
                .map(|face| {
                    let surface = &exact.surfaces()[face.surface.unwrap().index()];
                    Domain::new(&exact, face, surface, target, reading).map(|d| d.slack)
                })
                .collect::<Vec<_>>()
        };
        // Every face closes within the reading, the slack holding its gap.
        let admitted = domains(reading);
        assert!(admitted.iter().all(Result::is_ok), "{admitted:?}");
        assert!(
            admitted.iter().flatten().any(|&slack| slack >= e),
            "{admitted:?}"
        );
        // Without it, or with a thousandth of it, some face stays open.
        for short in [0.0, 1e-3 * reading] {
            assert!(
                domains(short).contains(&Err(open)),
                "{short:e}: {:?}",
                domains(short)
            );
        }
    }

    /// The bound of an all-planar boolean still holds its exact result
    /// when the mesh is off it (#314): one corner of the wall's mesh is
    /// moved, and every exact point, its corners included, is within the
    /// bound. The samples show the deviation, so the search settles at
    /// the budget as for a curved result; without a budget it tightens to
    /// near the sampled worst.
    #[test]
    fn a_planar_bound_holds_a_moved_corner() {
        let mut g = GeometryGraphBuilder::new();
        let rectangle = |g: &mut GeometryGraphBuilder, x, y| {
            g.push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x,
                    y,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap()
        };
        let extrude = |g: &mut GeometryGraphBuilder, profile, depth| {
            g.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth,
            }))
            .unwrap()
        };
        let profile = rectangle(&mut g, 3.0, 0.25);
        let wall = extrude(&mut g, profile, 3.0);
        let profile = rectangle(&mut g, 1.0, 2.0);
        let door = extrude(&mut g, profile, 1.0);
        let door = g
            .push(GeometryNode::Instance(Instance {
                source: door,
                transform: Transform3::from_translation(Vec3::new(0.3, 0.5, 0.0))
                    * Transform3::from_rotation_x(core::f64::consts::FRAC_PI_2),
            }))
            .unwrap();
        let cut = g
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: wall,
                right: door,
                operator: BooleanOperator::Difference,
            }))
            .unwrap();
        let graph = g.finish(vec![cut]).unwrap();
        let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
        let mut mesh = crate::ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_reported(&graph, cut, &options)
            .unwrap()
            .mesh;
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, cut, &options)
            .unwrap();
        // Every face is a plane: the search runs to the floor.
        let untouched = measure(&exact, &mesh, 1e-3, 1e-3, 0.0, 0.0).unwrap();
        assert!(untouched <= 1e-11, "{untouched:e}");
        // A deviation the report adds anyway (a slab width, #266) lets the
        // budget settle it, as a sampled one does.
        let shown = measure(&exact, &mesh, 1e-3, 1e-3, 1e-6, 0.0).unwrap();
        assert!(shown > 1e-6 && shown <= 1e-3, "{shown:e}");
        // Three mesh corners in turn, each moved off every face it is on.
        let corners = mesh.positions.len();
        for k in [0, corners / 3, 2 * corners / 3] {
            let original = mesh.positions[k];
            mesh.positions[k] = original + Vec3::new(1.0, -2.0, 1.5) * 1e-5;
            let budget = 1e-3;
            let bound = measure(&exact, &mesh, budget, budget, 0.0, 0.0).unwrap();
            let index = crate::certify::TriangleIndex::new(&mesh).unwrap();
            let topology = exact.topology();
            let mut worst: Scalar = 0.0;
            for face in topology.faces() {
                let surface = &exact.surfaces()[face.surface.unwrap().index()];
                for bound in &face.bounds {
                    let wire = &topology.loops()[bound.loop_id.index()];
                    for (e, edge_use) in wire.edges.iter().enumerate() {
                        let curve = &exact.curves2()[edge_use.pcurve.unwrap().index()];
                        let span = exact.pcurve_interval(bound.loop_id, e).unwrap();
                        for s in 0..=200 {
                            let t = span.start + (span.end - span.start) * Scalar::from(s) / 200.0;
                            let uv = axiolid_reference::evaluate2(curve, t).unwrap();
                            let p =
                                axiolid_reference::surface::evaluate(surface, uv.x, uv.y).unwrap();
                            worst = worst.max(index.nearest(p));
                        }
                    }
                }
            }
            assert!(
                worst <= bound && bound <= budget,
                "corner {k}: {worst:e} sampled, the bound {bound:e}"
            );
            // The corner's old place is on the exact surface (within
            // rounding): it shows the deviation.
            let moved = index.nearest(original);
            assert!(moved > 1e-6, "corner {k}: moved by {moved:e}");
            // Settled at the budget, not searched on towards the samples.
            assert!(
                bound > 10.0 * worst.max(moved),
                "corner {k}: {bound:e} for {worst:e} sampled, {moved:e} at the corner"
            );
            // Without a budget it tightens towards them.
            if k == 0 {
                let tight = measure(&exact, &mesh, budget, 0.0, 0.0, 0.0).unwrap();
                assert!(
                    worst <= tight && tight <= 3.0 * worst.max(moved),
                    "corner {k}: {tight:e} for {worst:e} sampled, {moved:e} at the corner"
                );
            }
            mesh.positions[k] = original;
        }
    }
}
