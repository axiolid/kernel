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
//!   within a rounding gap. Every point of the exact boundary is then within
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
//! # Unbounded, by name
//!
//! A boolean the exact compiler refuses (a union or intersection of placed
//! operands, an operand that is not a placed extrusion, a configuration the
//! general boolean cannot build) carries the refusal's name. A face whose
//! surface family has no derivative bounds, whose pcurve family has no
//! certified chord bound, or whose loops do not close in its parameters
//! (a periodic face trimmed by two separate circles) is unbounded by name.
//!
//! The exact result is compiled at the node's own options, so under an
//! instance it is the local solid and the mesh the local mesh; the
//! instance scales the bound like any other.

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::Curve2;
use axiolid_mesh::TriMesh;
use std::collections::HashSet;

use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};
use axiolid_reference::bound::{chord_bound2, continuity_breaks2, SurfaceBoundOracle};
use axiolid_surface::Surface;

use super::{Deviation, DeviationBound, DeviationPath};
use crate::certify::{certify_within, Cell, Coverage, Patch};
use crate::ReferenceExactCompiler;

/// The detail of a boolean measured against the exact compiler's result:
/// the exact boolean of the operands moved by at most the tolerance (#228).
const MEASURED: &str = "measured against the exact compiler's result, operands within tolerance";
/// The detail of a boolean measured against the exact compiler's result
/// when no decision within tolerance fired: the exact boolean of the
/// operands as given (#236).
const MEASURED_EXACT: &str =
    "measured against the exact compiler's result, the exact boolean of the given operands";
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
    let (exact, report) =
        match ReferenceExactCompiler::new().compile_exact_with_report(graph, id, options) {
            Ok(compiled) => compiled,
            Err(error) => {
                return Deviation::one(
                    DeviationPath::Boolean,
                    NO_EXACT,
                    DeviationBound::Unbounded(refusal(&error)),
                )
            }
        };
    // No decision within tolerance: the exact compiler's result is the
    // boolean of the operands as given, and the bound holds for it (#236).
    let detail = if report.is_exact() {
        MEASURED_EXACT
    } else {
        MEASURED
    };
    let target = crate::compiler::chord_error(options);
    let bound = match measure(&exact, mesh, target, target) {
        Ok(value) => DeviationBound::Certified(value),
        Err(reason) => DeviationBound::Unbounded(reason),
    };
    Deviation::one(DeviationPath::Boolean, detail, bound)
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
/// `enough` (see the module's "When to stop").
fn measure(
    exact: &ExactBRep,
    mesh: &TriMesh,
    target: Scalar,
    enough: Scalar,
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
        let domain = Domain::new(exact, face, surface, target)?;
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
    certify_within(&patches, cells, mesh, target, enough)
        .ok_or("boolean: an exact face could not be bounded against the mesh")
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
}

/// A face's domain in its surface parameters: the polygon of its loops'
/// flattened pcurves, every exact boundary point within `slack` of it.
struct Domain {
    segments: Vec<[Point2; 2]>,
    slack: Scalar,
    lo: Point2,
    hi: Point2,
}

impl Domain {
    fn new(
        exact: &ExactBRep,
        face: &axiolid_topology::Face<axiolid_brep::SurfaceId>,
        surface: &Surface,
        target: Scalar,
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
        let mut count = 0_usize;
        for spans in &uses {
            let mut ring = Vec::new();
            for &(curve, a, b) in spans {
                flatten(curve, a, b, delta, &mut ring, &mut count)?;
            }
            rings.push(ring);
        }
        let scale = (hi - lo).length().max(lo.length()).max(hi.length());
        let closing = 1e-9 * scale.max(1.0);
        let mut gap: Scalar = 0.0;
        let mut segments = Vec::new();
        for ring in &rings {
            // Each span contributes its own points, so consecutive spans
            // repeat (within rounding) the point they share; the closing
            // pair is the last point and the first.
            for pair in ring.windows(2) {
                segments.push([pair[0], pair[1]]);
            }
            if let (Some(&first), Some(&last)) = (ring.first(), ring.last()) {
                segments.push([last, first]);
            }
        }
        // The gaps between consecutive spans (and the ring's close) are
        // already segments above; measure them to refuse a loop that does
        // not close and to fold the rounding into the slack.
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
        let slack = (delta + gap) * (1.0 + 1e-9) + 4.0 * Scalar::EPSILON * scale;
        // The domain lies inside its boundary's box, and the boundary
        // within `slack` of the polygon: the polygon's box, grown, holds
        // it (the box above only sized `delta`).
        let (lo, hi) = rings.iter().flatten().fold(
            (
                Point2::splat(Scalar::INFINITY),
                Point2::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        Ok(Self {
            segments,
            slack,
            lo: lo - Point2::splat(slack),
            hi: hi + Point2::splat(slack),
        })
    }

    fn coverage(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Coverage {
        let (lo, hi) = (Point2::new(x.0, y.0), Point2::new(x.1, y.1));
        let near = self.segments.iter().any(|&[a, b]| {
            let (slo, shi) = (a.min(b), a.max(b));
            let apart = (lo - shi).max(slo - hi);
            if apart.x > self.slack || apart.y > self.slack {
                return false;
            }
            segment_box_distance(a, b, lo, hi) <= self.slack
        });
        if near {
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
/// `delta` of it, appended to `out`.
fn flatten(
    curve: &Curve2,
    a: Scalar,
    b: Scalar,
    delta: Scalar,
    out: &mut Vec<Point2>,
    count: &mut usize,
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
        split(curve, from, to, delta, 0, &mut params)?;
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
/// end parameter in order.
fn split(
    curve: &Curve2,
    a: Scalar,
    b: Scalar,
    delta: Scalar,
    depth: u32,
    params: &mut Vec<Scalar>,
) -> Result<(), &'static str> {
    let (lo, hi) = (a.min(b), a.max(b));
    let fits = chord_bound2(curve, lo, hi).is_some_and(|d| d <= delta);
    if fits {
        params.push(b);
        return Ok(());
    }
    if depth >= MAX_DEPTH {
        return Err("boolean: an exact face's pcurve family has no certified chord bound");
    }
    let m = 0.5 * (a + b);
    split(curve, a, m, delta, depth + 1, params)?;
    split(curve, m, b, delta, depth + 1, params)
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
        let bound = measure(&exact, &mesh, budget, 0.0).unwrap();
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
}
