//! A boolean's curved operand meshed by the angle it meets a plane at
//! (#380, ADR 0080).
//!
//! The mesh cut of two faceted operands moves the cut by the operands'
//! chord sag over the sine of the angle their surfaces meet at: where a
//! curved operand meets a plane of the other obliquely, the section curve
//! lands further from the exact one than either flattening alone (the box
//! intersected with a bent tube of #378 reports 1.0116e-3 against a budget
//! of 1e-3). So a curved operand is meshed at its chord budget times `s`,
//! the sine of the smallest angle at which it meets a **plane** face of the
//! other operand along their section.
//!
//! - **Closed forms.** `s` comes from a plane against a cylinder, cone,
//!   sphere or torus ([`plane_section_sine`]), over the curved face's
//!   parameter box where it matters (a cylinder's or cone's axial range, a
//!   torus's angle about its axis). Every other family has no form and
//!   contributes nothing.
//! - **Which faces.** The surfaces come from the operands' **leaves**,
//!   compiled exactly in their own frames (never the exact boolean) and
//!   placed by the chain of instances above them, plus the planes of
//!   half-space leaves ([`inventory`]). A pair is read only where the
//!   faces' enclosing balls meet. A boolean whose operands are certainly
//!   planar skips all of it ([`certainly_planar`]).
//! - **Which `s`.** The least over the pairs, floored at [`FLOOR`] (a
//!   grazing contact would refine without limit) and rounded down to a
//!   multiple of 1/256 ([`quantized`]), so constructions that differ by
//!   rounding get the same chord and share cache keys. An operand with no
//!   such pair keeps its chord: `1`. A plane that only touches the face,
//!   or reaches over to its far side by no more than the cut's linear
//!   tolerance ([`shallow_reach`]: the face's greatest distance to the
//!   shallower side of the plane, which is what the rounding of a
//!   placement far from the origin needs), is no pair: the mesh cut
//!   decides a contact within its tolerance by a skin within rounding and
//!   the exact referee (#243, #400), and a finer chord turns that into a
//!   pinch the cut refuses. The reach is read off the face, not off the
//!   section's angle: a slanted plane through a thin pipe meets it at a
//!   sine whose `radius (1 - cos)` is a fraction of the depth it cuts, and
//!   a plane tangent to a tube at one point may cross it at another. A
//!   crossing that holds a tangent point reads `0` there and gets the
//!   floor.
//! - **Nested booleans.** A boolean's factor shrinks what its operand is
//!   meshed at by that much *at the most*: a leaf under several booleans
//!   (a wall less a row of openings) gets the least of their factors, not
//!   their product (`Budget::shrunk` in the compiler).
//! - **Every compile.** The factors depend on the graph alone, never on
//!   whether a deviation report is requested, so the mesh does not.
//!
//! This is an accuracy measure, not a certificate: the boolean's report
//! still measures, and a pair with no closed form (two curved faces, an
//! elliptical cylinder, a B-spline, an authored mesh) or a contact grazing
//! below the floor (about 14.5 degrees) can still miss the budget and
//! reports so. Anything unreadable (a leaf the exact compiler refuses, a
//! chain that is not rigid) contributes nothing and never refuses the
//! compile. Of a bounded half-space only the plane is read, not the walls
//! of its boundary prism. The parameter box of a face is the hull of its
//! pcurves' ends, middles and quarter points, as the deviation search sizes
//! its domains: approximate, which suits a heuristic.

use std::collections::{HashMap, HashSet};
use std::f64::consts::TAU;

use axiolid_brep::ExactBRep;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point2, Point3, Scalar, Transform3, Vec3};
use axiolid_curve::Curve2;
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};
use axiolid_primitive::Primitive;
use axiolid_profile::Profile;
use axiolid_reference::bound::SurfaceBoundOracle;
use axiolid_surface::Surface;

use crate::ReferenceExactCompiler;

/// The smallest factor a curved operand's chord is shrunk by: at most four
/// times the triangles on a doubly curved operand, two on a singly curved.
pub(crate) const FLOOR: Scalar = 0.25;

/// `s` is rounded down to a multiple of `1 / GRID`.
const GRID: Scalar = 256.0;

/// The rounding noise a computed `s` is credited before it is rounded down:
/// a plane square to an axis reads `1 - 1e-16` in a rotated frame, which
/// must not fall to the grid step below.
const NOISE: Scalar = 1e-9;

/// Relative slack on a cosine or a feasibility test, so a tangent contact
/// read with rounding still meets.
const SLACK: Scalar = 1e-9;

/// A plane this near parallel to an axis (`|n x a|`) is read as parallel.
const PARALLEL: Scalar = 1e-12;

/// The factors of a boolean's `left` and `right` operands (in the order the
/// graph names them): `1` for an operand with no curved face meeting a
/// plane of the other, and for a region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Factors {
    pub(crate) left: Scalar,
    pub(crate) right: Scalar,
}

impl Factors {
    pub(crate) const ONE: Self = Self {
        left: 1.0,
        right: 1.0,
    };
}

/// Leaf B-reps compiled exactly in their own frames, per node and the bits
/// of the tolerance it was compiled at (the same node under instances that
/// scale the tolerance differently compiles once each), for one compile:
/// `None` where the exact compiler refused the leaf.
pub(crate) type LeafCache = HashMap<(NodeId, [u64; 2]), Option<ExactBRep>>;

/// The factors of the boolean `id`'s operands, or [`Factors::ONE`] for any
/// other node. Infallible: whatever is unreadable contributes nothing.
pub(crate) fn operand_factors(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    leaves: &mut LeafCache,
) -> Factors {
    let Some(GeometryNode::SolidOperation(SolidOperation::Boolean { left, right, .. })) =
        graph.get(id)
    else {
        return Factors::ONE;
    };
    if certainly_planar(graph, *left) && certainly_planar(graph, *right) {
        return Factors::ONE;
    }
    let left = inventory(graph, *left, options, leaves);
    let right = inventory(graph, *right, options, leaves);
    let reach = options.tolerance().linear();
    Factors {
        left: factor(&left, &right, reach),
        right: factor(&right, &left, reach),
    }
}

/// The factor of the operand whose faces are `curved`, against the planes
/// of the `other`: [`quantized`] of the least sine over the pairs whose
/// balls meet within `reach`, or `1` with none.
fn factor(curved: &[Face], other: &[Face], reach: Scalar) -> Scalar {
    let mut least: Option<Scalar> = None;
    for face in curved {
        let Shape::Curved { surface, lo, hi } = &face.shape else {
            continue;
        };
        for plane in other {
            let Shape::Plane { point, normal } = plane.shape else {
                continue;
            };
            if !face.meets(plane, reach) {
                continue;
            }
            let Some(s) = plane_section_sine(point, normal, surface, *lo, *hi) else {
                continue;
            };
            // A plane that only touches the face is no section: the mesh
            // cut decides a contact within its tolerance by a skin within
            // rounding and the exact referee (#243, #400), which a finer
            // chord turns into a pinch. It touches where the face lies to
            // one side of it, or reaches over to the other by no more than
            // the tolerance; a plane that cuts through the face with a
            // tangent point inside the crossing is a section, and reads the
            // floor. Where the reach is unknown the pair counts.
            if shallow_reach(point, normal, surface, *lo, *hi)
                .is_some_and(|reach_in| reach_in <= reach)
            {
                continue;
            }
            least = Some(least.map_or(s, |t| t.min(s)));
        }
    }
    least.map_or(1.0, quantized)
}

/// `s` as a chord factor: rounded down to a multiple of `1 / GRID`, at
/// least [`FLOOR`].
fn quantized(s: Scalar) -> Scalar {
    ((s * GRID + NOISE).floor() / GRID).max(FLOOR)
}

/// One face of an operand's leaves, in the boolean's frame.
struct Face {
    shape: Shape,
    /// Where the face lies, `None` where unbounded or not bounded.
    ball: Option<Ball>,
}

enum Shape {
    /// A plane through `point` with unit `normal`.
    Plane { point: Point3, normal: Vec3 },
    /// A cylinder, cone, sphere or torus over a parameter box.
    Curved {
        surface: Surface,
        lo: Point2,
        hi: Point2,
    },
}

/// A ball holding a face.
#[derive(Debug, Clone, Copy)]
struct Ball {
    centre: Point3,
    radius: Scalar,
}

impl Face {
    /// Whether the faces' balls meet within `reach`; an unbounded face
    /// meets everything.
    fn meets(&self, other: &Self, reach: Scalar) -> bool {
        match (self.ball, other.ball) {
            (Some(a), Some(b)) => a.centre.distance(b.centre) <= a.radius + b.radius + reach,
            _ => true,
        }
    }
}

/// The faces of everything under `id`, in the frame of the boolean above
/// it: the leaves compiled exactly and placed by the chain of instances
/// above them, and the planes of half-spaces. Empty where the walk is
/// longer than the graph is large (a shared sub-graph read through too many
/// paths), as for [`crate::compiler`]'s regions: unreadable contributes
/// nothing.
fn inventory(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    leaves: &mut LeafCache,
) -> Vec<Face> {
    let mut faces = Vec::new();
    let mut work = vec![(id, Transform3::IDENTITY)];
    let mut visits = 0_usize;
    while let Some((mut id, mut transform)) = work.pop() {
        visits += 1;
        if visits > graph.len() {
            return Vec::new();
        }
        while let Some(node) = graph.get(id) {
            match node {
                GeometryNode::Instance(instance) => {
                    transform *= instance.transform;
                    id = instance.source;
                }
                GeometryNode::InstanceAtStation(placed) => {
                    let Ok(placement) = crate::station::placement(graph, id) else {
                        break;
                    };
                    transform *= placement.transform;
                    id = placed.source;
                }
                GeometryNode::SolidOperation(SolidOperation::Boolean { left, right, .. }) => {
                    work.push((*right, transform));
                    work.push((*left, transform));
                    break;
                }
                GeometryNode::Collection(members) => {
                    work.extend(members.iter().rev().map(|&member| (member, transform)));
                    break;
                }
                GeometryNode::HalfSpace(half_space) => {
                    faces.extend(half_space_plane(half_space, &transform));
                    break;
                }
                GeometryNode::SolidOperation(operation) => {
                    if let Some((half_space, ..)) = crate::half_space_boundary::bounded(operation) {
                        if let Some(GeometryNode::HalfSpace(half_space)) = graph.get(half_space) {
                            faces.extend(half_space_plane(half_space, &transform));
                        }
                    } else {
                        leaf_faces(graph, id, &transform, options, leaves, &mut faces);
                    }
                    break;
                }
                GeometryNode::TriMesh(_) | GeometryNode::PolygonMesh(_) => break,
                _ => {
                    leaf_faces(graph, id, &transform, options, leaves, &mut faces);
                    break;
                }
            }
        }
    }
    faces
}

/// The boundary plane of `half_space` placed by `transform`, or `None` for
/// a placement with no inverse. A normal moves by the inverse transpose.
fn half_space_plane(
    half_space: &axiolid_primitive::HalfSpace,
    transform: &Transform3,
) -> Option<Face> {
    let matrix = transform.matrix3;
    if matrix.determinant() == 0.0 || !matrix.is_finite() {
        return None;
    }
    let normal = (matrix.inverse().transpose() * half_space.boundary.normal).normalize_or_zero();
    let point = transform.transform_point3(half_space.boundary.origin);
    (normal != Vec3::ZERO && point.is_finite()).then_some(Face {
        shape: Shape::Plane { point, normal },
        ball: None,
    })
}

/// The faces of the leaf `id` compiled exactly (once, in its own frame) and
/// placed by `transform`; none where either refuses.
fn leaf_faces(
    graph: &GeometryGraph,
    id: NodeId,
    transform: &Transform3,
    options: &ExecutionOptions,
    leaves: &mut LeafCache,
    faces: &mut Vec<Face>,
) {
    let tolerance = options.tolerance();
    let Some(exact) = leaves
        .entry((
            id,
            [tolerance.linear().to_bits(), tolerance.angular().to_bits()],
        ))
        .or_insert_with(|| {
            ReferenceExactCompiler::new()
                .compile_exact(graph, id, options)
                .ok()
        })
        .as_ref()
    else {
        return;
    };
    let placed;
    let exact = if *transform == Transform3::IDENTITY {
        exact
    } else {
        match exact.transformed(transform) {
            Ok(moved) => {
                placed = moved;
                &placed
            }
            Err(_) => return,
        }
    };
    brep_faces(exact, faces);
}

/// The planes and the curved faces `plane_section_sine` has a form for.
fn brep_faces(exact: &ExactBRep, out: &mut Vec<Face>) {
    for face in exact.topology().faces() {
        let Some(surface) = face.surface.and_then(|id| exact.surfaces().get(id.index())) else {
            continue;
        };
        let bounds = parameter_box(exact, face);
        let shape = match surface {
            Surface::Plane(plane) => {
                let normal = plane.frame.z.normalize_or_zero();
                if normal == Vec3::ZERO || !plane.frame.origin.is_finite() {
                    continue;
                }
                Shape::Plane {
                    point: plane.frame.origin,
                    normal,
                }
            }
            Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_) => {
                let Some((lo, hi)) = bounds else {
                    continue;
                };
                Shape::Curved {
                    surface: surface.clone(),
                    lo,
                    hi,
                }
            }
            _ => continue,
        };
        let ball = bounds.and_then(|(lo, hi)| ball(surface, lo, hi));
        out.push(Face { shape, ball });
    }
}

/// The hull of the face's pcurves' ends, middles and quarter points.
fn parameter_box(
    exact: &ExactBRep,
    face: &axiolid_topology::Face<axiolid_brep::SurfaceId>,
) -> Option<(Point2, Point2)> {
    let topology = exact.topology();
    let mut lo = Point2::splat(Scalar::INFINITY);
    let mut hi = Point2::splat(Scalar::NEG_INFINITY);
    for bound in &face.bounds {
        let wire = topology.loops().get(bound.loop_id.index())?;
        for (k, edge_use) in wire.edges.iter().enumerate() {
            let curve = edge_use
                .pcurve
                .and_then(|id| exact.curves2().get(id.index()))?;
            let span = exact.pcurve_interval(bound.loop_id, k)?;
            for share in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let t = span.start + (span.end - span.start) * share;
                let p = axiolid_reference::evaluate2(curve, t).ok()?;
                lo = lo.min(p);
                hi = hi.max(p);
            }
        }
    }
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

/// A ball around the surface over the box: its middle, and the first and
/// second derivative bounds times the half extents (a Taylor bound). `None`
/// where the surface cannot be bounded there.
fn ball(surface: &Surface, lo: Point2, hi: Point2) -> Option<Ball> {
    let bounds = SurfaceBoundOracle::new(surface)?.bounds((lo.x, hi.x), (lo.y, hi.y))?;
    if !bounds.smooth {
        return None;
    }
    let (hx, hy) = (0.5 * (hi.x - lo.x), 0.5 * (hi.y - lo.y));
    let radius = bounds.du * hx
        + bounds.dv * hy
        + 0.5 * (bounds.duu * hx * hx + 2.0 * bounds.duv * hx * hy + bounds.dvv * hy * hy);
    let centre = axiolid_reference::surface::evaluate(surface, lo.x + hx, lo.y + hy).ok()?;
    (radius.is_finite() && centre.is_finite()).then_some(Ball { centre, radius })
}

/// Whether the operand `id` is made of planes alone, read from the graph
/// without any construction: instances and booleans of such, half-spaces,
/// meshes, blocks, and extrusions of a polygonal profile.
fn certainly_planar(graph: &GeometryGraph, id: NodeId) -> bool {
    let mut seen = HashSet::new();
    let mut work = vec![id];
    while let Some(id) = work.pop() {
        if !seen.insert(id) {
            continue;
        }
        let planar = |profile: NodeId| matches!(graph.get(profile), Some(GeometryNode::Profile(p)) if polygonal(p));
        match graph.get(id) {
            Some(GeometryNode::Instance(instance)) => work.push(instance.source),
            Some(GeometryNode::InstanceAtStation(placed)) => work.push(placed.source),
            Some(GeometryNode::Collection(members)) => work.extend(members),
            Some(GeometryNode::SolidOperation(SolidOperation::Boolean { left, right, .. })) => {
                work.extend([*left, *right]);
            }
            Some(
                GeometryNode::HalfSpace(_)
                | GeometryNode::TriMesh(_)
                | GeometryNode::PolygonMesh(_)
                | GeometryNode::Primitive(Primitive::Block { .. }),
            ) => {}
            Some(GeometryNode::SolidOperation(operation)) => {
                let this = match operation {
                    SolidOperation::Extrusion { profile, .. } => planar(*profile),
                    SolidOperation::TaperedExtrusion {
                        start_profile,
                        end_profile,
                        ..
                    } => planar(*start_profile) && planar(*end_profile),
                    other => crate::half_space_boundary::bounded(other).is_some(),
                };
                if !this {
                    return false;
                }
            }
            _ => return false,
        }
    }
    true
}

/// Whether every edge of `profile` is a straight line.
fn polygonal(profile: &Profile) -> bool {
    match profile {
        Profile::Rectangle(r) => r.outer_radius.is_none() && r.inner_radius.is_none(),
        Profile::Contour(contour) => {
            core::iter::once(&contour.outer)
                .chain(&contour.holes)
                .all(|ring| {
                    ring.segments
                        .iter()
                        .all(|s| matches!(s.curve, Curve2::Line(_) | Curve2::Polyline(_)))
                })
        }
        Profile::Derived { basis, .. } => polygonal(basis),
        Profile::Composite(members) => members.iter().all(polygonal),
        _ => false,
    }
}

/// A plane against the frame of a revolved surface: its unit axis `a`, and
/// with `n` the plane's unit normal and `p` a point on it, `alpha = n . a`,
/// `m` the part of `n` square to the axis, at angle `psi_n` about it, and
/// `delta = (p - o) . n` the plane's offset from the frame's origin.
struct Relation {
    alpha: Scalar,
    m: Scalar,
    psi_n: Scalar,
    delta: Scalar,
}

impl Relation {
    fn new(frame: &Frame3, p: Point3, n: Vec3) -> Option<Self> {
        let (x, y, a) = (
            frame.x.normalize_or_zero(),
            frame.y.normalize_or_zero(),
            frame.z.normalize_or_zero(),
        );
        let n = n.normalize_or_zero();
        if [x, y, a, n].contains(&Vec3::ZERO) || !p.is_finite() || !frame.origin.is_finite() {
            return None;
        }
        let (nx, ny) = (n.dot(x), n.dot(y));
        Some(Self {
            alpha: n.dot(a),
            m: nx.hypot(ny),
            psi_n: ny.atan2(nx),
            delta: (p - frame.origin).dot(n),
        })
    }
}

/// The sine of the smallest angle at which the surface `curved` meets the
/// plane through `p` with normal `n`, along their section inside the
/// surface's parameter box `[lo, hi]`; `None` when they do not meet there,
/// or when the family has no closed form (see the module's closed forms).
///
/// With `N` the surface's unit normal, the angle's cosine is `|n . N|`. On
/// a surface of revolution about `a`, a point at angle `psi` has
/// `n . e(psi) = m cos(psi - psi_n)`, written `m c`.
fn plane_section_sine(
    p: Point3,
    n: Vec3,
    curved: &Surface,
    lo: Point2,
    hi: Point2,
) -> Option<Scalar> {
    let u = (lo.x.min(hi.x), lo.x.max(hi.x));
    let v = (lo.y.min(hi.y), lo.y.max(hi.y));
    let positive = |value: Scalar| (value.is_finite() && value > 0.0).then_some(value);
    match curved {
        Surface::Cylinder(c) => cylinder(
            &Relation::new(&c.frame, p, n)?,
            positive(c.radius.abs())?,
            v,
        ),
        Surface::Cone(c) => cone(&Relation::new(&c.frame, p, n)?, (c.radius, c.semi_angle), v),
        Surface::Sphere(s) => sphere(&Relation::new(&s.frame, p, n)?, positive(s.radius.abs())?),
        Surface::Torus(t) => torus(
            &Relation::new(&t.frame, p, n)?,
            (positive(t.major_radius)?, positive(t.minor_radius)?),
            (u, v),
        ),
        _ => None,
    }
}

/// Whether `angle` lies in the arc `[lo, hi]` of a turn (`whole` where the
/// arc is the full turn).
fn within(angle: Scalar, (lo, hi): (Scalar, Scalar), whole: bool) -> bool {
    let along = (angle - lo).rem_euclid(TAU);
    whole || along <= hi - lo + 1e-12 || along >= TAU - 1e-12
}

/// `sqrt(1 - cos^2)`, the sine of an angle whose cosine is `cos`.
fn sine(cos: Scalar) -> Scalar {
    (1.0 - cos * cos).max(0.0).sqrt()
}

/// The interval of cosines `c = cos(psi - psi_n)` a plane meets a surface at
/// between two ends (where `c` is monotone), within `[-1, 1]`; `None` when
/// it lies outside.
fn cosines(ends: [Scalar; 2]) -> Option<(Scalar, Scalar)> {
    let (lo, hi) = (ends[0].min(ends[1]), ends[0].max(ends[1]));
    let (lo, hi) = (lo.max(-1.0 - SLACK), hi.min(1.0 + SLACK));
    (lo <= hi).then_some((lo.clamp(-1.0, 1.0), hi.clamp(-1.0, 1.0)))
}

/// A cylinder of radius `r` over its axial range `v`: the section is
/// `r m c + v alpha = delta`, and the normal is radial, `n . N = m c`. The
/// angle is least where `|c|` is greatest within the range the axial extent
/// allows, so a plane nearly parallel to the axis is read over the face
/// only.
fn cylinder(rel: &Relation, r: Scalar, v: (Scalar, Scalar)) -> Option<Scalar> {
    if rel.m <= PARALLEL {
        // Square to the axis: a circle at one height, met square.
        let height = rel.delta / rel.alpha;
        return (v.0 <= height && height <= v.1).then_some(1.0);
    }
    let at = |v: Scalar| (rel.delta - v * rel.alpha) / (r * rel.m);
    let (lo, hi) = cosines([at(v.0), at(v.1)])?;
    Some(sine(rel.m * lo.abs().max(hi.abs())))
}

/// The axial range `v` of a cone of base radius `r0` and semi-angle `beta`
/// where its radius `r0 + v tan(beta)` is positive; `None` where none is.
fn cone_span((r0, beta): (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<(Scalar, Scalar)> {
    let slope = beta.tan();
    let (mut v0, mut v1) = v;
    if slope > 0.0 {
        v0 = v0.max(-r0 / slope);
    } else if slope < 0.0 {
        v1 = v1.min(-r0 / slope);
    }
    (v0 <= v1).then_some((v0, v1))
}

/// A cone of base radius `r0` and semi-angle `beta` over its axial range
/// `v` (kept where its radius `r0 + v tan(beta)` is positive): the section
/// is `rho(v) m c + v alpha = delta`, and `n . N = cos(beta) m c -
/// sin(beta) alpha`. `c` is monotone in `v` (no pole where `rho > 0`), so
/// the angle is least at an end of the range it spans.
fn cone(rel: &Relation, (r0, beta): (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<Scalar> {
    if !(r0.is_finite() && beta.is_finite()) {
        return None;
    }
    let (sin, cos) = beta.sin_cos();
    let slope = beta.tan();
    let (v0, v1) = cone_span((r0, beta), v)?;
    if rel.m <= PARALLEL {
        let height = rel.delta / rel.alpha;
        return (v0 <= height && height <= v1).then(|| sine(sin * rel.alpha));
    }
    let at = |v: Scalar| {
        let rho = (r0 + v * slope).max(Scalar::MIN_POSITIVE);
        (rel.delta - v * rel.alpha) / (rel.m * rho)
    };
    let (lo, hi) = cosines([at(v0), at(v1)])?;
    let cosine = |c: Scalar| (cos * rel.m * c - sin * rel.alpha).abs();
    Some(sine(cosine(lo).max(cosine(hi))))
}

/// A sphere of radius `r`: a plane at `delta` from its centre meets it at
/// one angle everywhere, `n . N = delta / r`.
fn sphere(rel: &Relation, r: Scalar) -> Option<Scalar> {
    let cos = rel.delta / r;
    (cos.abs() <= 1.0 + SLACK).then(|| sine(cos))
}

/// A torus of major radius `big` and minor radius `small` over the angles
/// `u` about its axis and `v` about its tube: with `g(psi) = delta - big m
/// c` the section has `n . N = g / small`, and exists where `|g| <= small
/// sqrt(m^2 c^2 + alpha^2)`. That is an interval of `c`, between the roots
/// of `(big^2 - small^2) m^2 c^2 - 2 delta big m c + delta^2 - small^2
/// alpha^2`, and `|g|` is linear in `c`: along the section the angle is
/// least where `c` is extreme, at a root (the section turns in `psi`), at
/// `c = +-1`, at an end of `u`, or at an end of `v`; each counts where it
/// lies in the face (a tangent point on a part the face lacks is no pair).
fn torus(
    rel: &Relation,
    (big, small): (Scalar, Scalar),
    (u, v): ((Scalar, Scalar), (Scalar, Scalar)),
) -> Option<Scalar> {
    if big <= small {
        return None;
    }
    let (whole_u, whole_v) = (u.1 - u.0 >= TAU - 1e-9, v.1 - v.0 >= TAU - 1e-9);
    // The cosine of the angle between the plane's and the surface's normals
    // at minor angle `phi`, where the section passes at cosine `c`.
    let normal = |c: Scalar, phi: Scalar| phi.cos() * rel.m * c + phi.sin() * rel.alpha;
    let mut sines: Vec<Scalar> = Vec::new();
    // The sections at cosine `c` (at angles `psis` about the axis): at each,
    // the minor angles `phi` with `(big + small cos phi) m c + small sin phi
    // alpha = delta` that lie in the face.
    let mut at = |c: Scalar, psis: &[Scalar]| {
        if !psis.iter().any(|&psi| within(psi, u, whole_u)) {
            return;
        }
        let reach = (rel.m * rel.m * c * c + rel.alpha * rel.alpha).sqrt();
        let q = (rel.delta - big * rel.m * c) / small;
        if reach <= PARALLEL {
            // The normal is square to every ring: any minor angle meets it
            // square, if the plane is there at all.
            if q.abs() <= SLACK {
                sines.push(1.0);
            }
            return;
        }
        if q.abs() > reach * (1.0 + SLACK) {
            return;
        }
        let spread = (q / reach).clamp(-1.0, 1.0).acos();
        let phase = rel.alpha.atan2(rel.m * c);
        for phi in [phase + spread, phase - spread] {
            if within(phi, v, whole_v) {
                sines.push(sine(normal(c, phi)));
            }
        }
    };
    if rel.m <= PARALLEL {
        // Square to the axis: the same minor angles at every `psi`.
        let sin_phi = rel.delta / (small * rel.alpha);
        if sin_phi.abs() <= 1.0 + SLACK {
            let phi = sin_phi.clamp(-1.0, 1.0).asin();
            if [phi, core::f64::consts::PI - phi]
                .into_iter()
                .any(|phi| within(phi, v, whole_v))
            {
                sines.push(sine(rel.delta / small));
            }
        }
    } else {
        let a = (big * big - small * small) * rel.m * rel.m;
        let b = -2.0 * rel.delta * big * rel.m;
        let c = rel.delta * rel.delta - small * small * rel.alpha * rel.alpha;
        let psis = |cos: Scalar| {
            let angle = cos.clamp(-1.0, 1.0).acos();
            [rel.psi_n + angle, rel.psi_n - angle]
        };
        let discriminant = b * b - 4.0 * a * c;
        if discriminant >= 0.0 {
            let root = discriminant.sqrt();
            for cos in [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)] {
                if cos.abs() <= 1.0 + SLACK {
                    at(cos.clamp(-1.0, 1.0), &psis(cos));
                }
            }
        }
        at(1.0, &[rel.psi_n]);
        at(-1.0, &[rel.psi_n + core::f64::consts::PI]);
        if !whole_u {
            for psi in [u.0, u.1] {
                at((psi - rel.psi_n).cos(), &[psi]);
            }
        }
        if !whole_v {
            // The section along an end of the minor angle.
            for phi in [v.0, v.1] {
                let ring = big + small * phi.cos();
                let cos = (rel.delta - small * phi.sin() * rel.alpha) / (ring * rel.m);
                if ring > 0.0 && cos.abs() <= 1.0 + SLACK {
                    let cos = cos.clamp(-1.0, 1.0);
                    if psis(cos).iter().any(|&psi| within(psi, u, whole_u)) {
                        sines.push(sine(normal(cos, phi)));
                    }
                }
            }
        }
    }
    sines.into_iter().min_by(Scalar::total_cmp)
}

/// How far the face over the box `[lo, hi]` reaches to the shallower side of
/// the plane through `p` with normal `n`: the least of its greatest distances
/// to either side, `0` where it lies to one side (a plane that only touches
/// it). That is the true penetration of a plane that grazes the surface and
/// the depth of the cap it cuts off, where a radius and an angle would read
/// only a plane parallel to a cylinder's axis, or a sphere. `None` where the
/// family has no form here: the pair then counts.
///
/// The box over-reaches a trimmed face, which only counts more pairs.
fn shallow_reach(p: Point3, n: Vec3, curved: &Surface, lo: Point2, hi: Point2) -> Option<Scalar> {
    let u = (lo.x.min(hi.x), lo.x.max(hi.x));
    let v = (lo.y.min(hi.y), lo.y.max(hi.y));
    let (least, greatest) = match curved {
        Surface::Cylinder(c) => {
            cylinder_offsets(&Relation::new(&c.frame, p, n)?, c.radius.abs(), v)
        }
        Surface::Cone(c) => {
            cone_offsets(&Relation::new(&c.frame, p, n)?, (c.radius, c.semi_angle), v)?
        }
        Surface::Sphere(s) => {
            let rel = Relation::new(&s.frame, p, n)?;
            let r = s.radius.abs();
            (-r - rel.delta, r - rel.delta)
        }
        Surface::Torus(t) => torus_offsets(
            &Relation::new(&t.frame, p, n)?,
            (t.major_radius, t.minor_radius),
            (u, v),
        )?,
        _ => return None,
    };
    let reach = greatest.max(0.0).min((-least).max(0.0));
    reach.is_finite().then_some(reach)
}

/// The range of `n . (x - p)` over a cylinder of radius `r` and axial range
/// `v`: `r m c + v alpha - delta`, with `c` over a turn.
fn cylinder_offsets(rel: &Relation, r: Scalar, v: (Scalar, Scalar)) -> (Scalar, Scalar) {
    let along = [v.0 * rel.alpha, v.1 * rel.alpha];
    let (least, greatest) = (along[0].min(along[1]), along[0].max(along[1]));
    (
        least - r * rel.m - rel.delta,
        greatest + r * rel.m - rel.delta,
    )
}

/// The range of `n . (x - p)` over a cone: `rho(v) m c + v alpha - delta`,
/// linear in `v` for `c = +-1`, so read at the ends of its span.
fn cone_offsets(
    rel: &Relation,
    (r0, beta): (Scalar, Scalar),
    v: (Scalar, Scalar),
) -> Option<(Scalar, Scalar)> {
    if !(r0.is_finite() && beta.is_finite()) {
        return None;
    }
    let slope = beta.tan();
    let (v0, v1) = cone_span((r0, beta), v)?;
    let ring = |v: Scalar| (r0 + v * slope).max(0.0) * rel.m;
    let (mut least, mut greatest) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
    for v in [v0, v1] {
        least = least.min(v * rel.alpha - ring(v) - rel.delta);
        greatest = greatest.max(v * rel.alpha + ring(v) - rel.delta);
    }
    Some((least, greatest))
}

/// The range of `n . (x - p)` over a torus face: `(R + r cos phi) m cos(psi -
/// psi_n) + r sin phi alpha - delta`. With `R > r` the ring is positive, so
/// at each minor angle `phi` the greatest is at the greatest cosine of the
/// angle about the axis the face holds, the least at its least; over `phi` each
/// is then a sinusoid, read at its ends and its peak.
fn torus_offsets(
    rel: &Relation,
    (big, small): (Scalar, Scalar),
    (u, v): ((Scalar, Scalar), (Scalar, Scalar)),
) -> Option<(Scalar, Scalar)> {
    if !(big.is_finite() && small.is_finite() && small > 0.0 && big > small) {
        return None;
    }
    let whole_u = u.1 - u.0 >= TAU - 1e-9;
    let whole_v = v.1 - v.0 >= TAU - 1e-9;
    let (c_least, c_greatest) = if whole_u {
        (-1.0, 1.0)
    } else {
        let at = |psi: Scalar| (psi - rel.psi_n).cos();
        let (a, b) = (at(u.0), at(u.1));
        (
            if within(rel.psi_n + core::f64::consts::PI, u, false) {
                -1.0
            } else {
                a.min(b)
            },
            if within(rel.psi_n, u, false) {
                1.0
            } else {
                a.max(b)
            },
        )
    };
    // The greatest of `a cos(phi) + b sin(phi)` over the minor angles.
    let sinusoid = |a: Scalar, b: Scalar| {
        let end = |phi: Scalar| a * phi.cos() + b * phi.sin();
        if whole_v || within(b.atan2(a), v, false) {
            a.hypot(b)
        } else {
            end(v.0).max(end(v.1))
        }
    };
    let (k_least, k_greatest) = (rel.m * c_least, rel.m * c_greatest);
    Some((
        big * k_least - small * sinusoid(-k_least, -rel.alpha) - rel.delta,
        big * k_greatest + small * sinusoid(k_greatest, rel.alpha) - rel.delta,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{BooleanOperator, Plane3, Tolerance, Vec2};
    use axiolid_model::{GeometryGraphBuilder, Instance};
    use axiolid_profile::{CircleProfile, RectangleProfile};
    use axiolid_surface::{Cone, Cylinder, Sphere, Torus};
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_6};

    fn frame(origin: Point3, x: Vec3, y: Vec3, z: Vec3) -> Frame3 {
        Frame3 { origin, x, y, z }
    }

    fn upright() -> Frame3 {
        frame(Point3::ZERO, Vec3::X, Vec3::Y, Vec3::Z)
    }

    /// The torus of the issue's bent tube: `R = 1`, `r = 0.1`, about an
    /// axis through `(-1, 0, 0)` along `-y`, over a quarter turn.
    fn bent_tube() -> Surface {
        Surface::Torus(Torus {
            frame: frame(Point3::new(-1.0, 0.0, 0.0), Vec3::X, Vec3::Z, Vec3::NEG_Y),
            major_radius: 1.0,
            minor_radius: 0.1,
        })
    }

    fn quarter() -> (Point2, Point2) {
        (Point2::new(0.0, 0.0), Point2::new(FRAC_PI_2, TAU))
    }

    fn column(radius: Scalar) -> Surface {
        Surface::Cylinder(Cylinder {
            frame: upright(),
            radius,
        })
    }

    fn sine_of(
        p: Point3,
        n: Vec3,
        surface: &Surface,
        (lo, hi): (Point2, Point2),
    ) -> Option<Scalar> {
        plane_section_sine(p, n, surface, lo, hi)
    }

    #[test]
    fn a_plane_through_the_bent_tube_meets_it_at_27_degrees() {
        let s = sine_of(
            Point3::new(-0.2, 0.0, 0.0),
            Vec3::X,
            &bent_tube(),
            quarter(),
        )
        .expect("a section");
        assert!((s - 17.0_f64.sqrt() / 9.0).abs() < 1e-12, "{s}");
        assert!((s - 0.458123).abs() < 1e-6, "{s}");
        assert_eq!(quantized(s), 117.0 / 256.0);
        assert_eq!(quantized(s), 0.45703125);
    }

    #[test]
    fn a_plane_along_the_tube_axis_meets_the_base_square() {
        let s = sine_of(Point3::ZERO, Vec3::Z, &bent_tube(), quarter()).expect("a section");
        assert!((s - 1.0).abs() < 1e-12, "{s}");
    }

    #[test]
    fn an_oblique_plane_cuts_a_cylinder_at_its_tilt() {
        let tilt = FRAC_PI_6 * 2.0;
        let n = Vec3::new(tilt.sin(), 0.0, tilt.cos());
        let s = sine_of(
            Point3::new(0.0, 0.0, 2.0),
            n,
            &column(0.2),
            (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0)),
        )
        .expect("a section");
        assert!((s - 0.5).abs() < 1e-12, "{s}");
    }

    #[test]
    fn a_nearly_parallel_plane_is_read_over_the_face_only() {
        // `alpha = 0.01`: met everywhere along an unrestricted axis at an
        // angle of `0.01`, but only a stretch of it lies in the face, where
        // `c` reaches `0.5` at the most.
        let alpha: Scalar = 0.01;
        let n = Vec3::new((1.0 - alpha * alpha).sqrt(), 0.0, alpha);
        let s = sine_of(
            Point3::new(0.1 / n.x, 0.0, 0.0),
            n,
            &column(0.2),
            (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0)),
        )
        .expect("a section");
        assert!((s - 0.866025).abs() < 1e-6, "{s}");
    }

    #[test]
    fn a_sphere_cut_off_centre() {
        let sphere = Surface::Sphere(Sphere {
            frame: upright(),
            radius: 0.5,
        });
        let s = sine_of(
            Point3::new(0.0, 0.0, 0.3),
            Vec3::Z,
            &sphere,
            (Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2)),
        )
        .expect("a section");
        assert!((s - 0.8).abs() < 1e-12, "{s}");
        assert_eq!(quantized(s), 0.796875);
    }

    #[test]
    fn a_cone_cut_square_to_its_axis() {
        let cone = Surface::Cone(Cone {
            frame: upright(),
            radius: 0.5,
            semi_angle: FRAC_PI_6,
        });
        let s = sine_of(
            Point3::new(0.0, 0.0, 0.5),
            Vec3::Z,
            &cone,
            (Point2::new(0.0, 0.0), Point2::new(TAU, 1.0)),
        )
        .expect("a section");
        assert!((s - 0.866025).abs() < 1e-6, "{s}");
    }

    #[test]
    fn a_cone_cut_along_its_axis_reads_the_slant() {
        // A plane parallel to the axis at 0.2 meets the cone where
        // `rho c = 0.2`: `c` falls from 0.4 at the base to 0.186 at the top,
        // and `n . N = cos(beta) c`, greatest at the base.
        let cone = Surface::Cone(Cone {
            frame: upright(),
            radius: 0.5,
            semi_angle: FRAC_PI_6,
        });
        let s = sine_of(
            Point3::new(0.2, 0.0, 0.0),
            Vec3::X,
            &cone,
            (Point2::new(0.0, 0.0), Point2::new(TAU, 1.0)),
        )
        .expect("a section");
        let cos = FRAC_PI_6.cos() * 0.4;
        assert!((s - (1.0 - cos * cos).sqrt()).abs() < 1e-12, "{s}");
    }

    #[test]
    fn a_tangent_plane_reads_zero_and_a_shallow_one_the_floor() {
        let s = sine_of(
            Point3::new(0.1, 0.0, 0.0),
            Vec3::X,
            &column(0.1),
            (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0)),
        );
        assert_eq!(s, Some(0.0));
        assert_eq!(quantized(0.0), FLOOR);
        assert_eq!(FLOOR, 0.25);
        // A plane 1 % inside the radius grazes at 8 degrees: the floor too.
        let s = sine_of(
            Point3::new(0.099, 0.0, 0.0),
            Vec3::X,
            &column(0.1),
            (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0)),
        )
        .expect("a section");
        assert!((s - (1.0 - 0.99_f64 * 0.99).sqrt()).abs() < 1e-12, "{s}");
        assert_eq!(quantized(s), FLOOR);
    }

    #[test]
    fn a_plane_that_misses_returns_none() {
        let whole = (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0));
        assert_eq!(
            sine_of(Point3::new(0.2, 0.0, 0.0), Vec3::X, &column(0.1), whole),
            None
        );
        let sphere = Surface::Sphere(Sphere {
            frame: upright(),
            radius: 0.5,
        });
        assert_eq!(
            sine_of(Point3::new(0.0, 0.0, 0.6), Vec3::Z, &sphere, whole),
            None
        );
        assert_eq!(
            sine_of(Point3::new(0.2, 0.0, 0.0), Vec3::X, &bent_tube(), quarter()),
            None
        );
        // Past the face's axial range, though not the unbounded surface.
        assert_eq!(
            sine_of(Point3::new(0.0, 0.0, 5.0), Vec3::Z, &column(0.1), whole),
            None
        );
    }

    #[test]
    fn a_square_plane_in_a_turned_frame_is_met_square_to_the_grid() {
        // Rounding in a placed frame must not drop `1` to the step below.
        let turn = Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2))
            * Transform3::from_rotation_z(0.6);
        let moved = Surface::Cylinder(Cylinder {
            frame: frame(
                turn.transform_point3(Point3::ZERO),
                turn.transform_vector3(Vec3::X),
                turn.transform_vector3(Vec3::Y),
                turn.transform_vector3(Vec3::Z),
            ),
            radius: 0.4,
        });
        let n = turn.transform_vector3(Vec3::Z);
        let s = sine_of(
            turn.transform_point3(Point3::new(0.0, 0.0, 1.0)),
            n,
            &moved,
            (Point2::new(0.0, 0.0), Point2::new(TAU, 2.0)),
        )
        .expect("a section");
        assert_eq!(quantized(s), 1.0);
    }

    #[test]
    fn a_torus_face_cut_short_keeps_the_roots_it_holds() {
        // A face of the bent tube cut short at 30 degrees: the roots at
        // 27.27 and 43.4 degrees are of the full quarter; of this one only
        // the first lies in it.
        let (lo, hi) = (Point2::new(0.0, 0.0), Point2::new(0.5, TAU));
        let s = sine_of(Point3::new(-0.2, 0.0, 0.0), Vec3::X, &bent_tube(), (lo, hi))
            .expect("a section");
        assert!((s - 17.0_f64.sqrt() / 9.0).abs() < 1e-12, "{s}");
        // Cut short before it: no section.
        let (lo, hi) = (Point2::new(0.0, 0.0), Point2::new(0.4, TAU));
        assert_eq!(
            sine_of(Point3::new(-0.2, 0.0, 0.0), Vec3::X, &bent_tube(), (lo, hi)),
            None
        );
    }

    /// The least sine of a cone or a torus against the plane, found
    /// independently of the closed forms: sweep the angle `psi` about the
    /// axis, solve the plane's equation for the other parameter with the
    /// surface's explicit points and normals, and take the least. The sweep
    /// misses the exact root of a feasibility interval by its step, so it
    /// reads at most a little above the closed form.
    fn swept(
        p: Point3,
        n: Vec3,
        curved: &Surface,
        u: (Scalar, Scalar),
        v: (Scalar, Scalar),
    ) -> Scalar {
        let (frame, o) = match curved {
            Surface::Cone(c) => (c.frame, c.frame.origin),
            Surface::Torus(t) => (t.frame, t.frame.origin),
            _ => unreachable!(),
        };
        let n = n.normalize();
        let mut least = Scalar::INFINITY;
        let steps = 400_000;
        for k in 0..=steps {
            let psi = u.0 + (u.1 - u.0) * k as Scalar / steps as Scalar;
            let e = frame.x * psi.cos() + frame.y * psi.sin();
            let (en, an) = (e.dot(n), frame.z.dot(n));
            let mut found = Vec::new();
            match curved {
                Surface::Cone(c) => {
                    // (o + rho e + v a - p) . n = 0 is linear in v.
                    let slope = c.semi_angle.tan();
                    let a = (o - p).dot(n) + c.radius * en;
                    let b = slope * en + an;
                    if b.abs() > 1e-12 {
                        let height = -a / b;
                        let rho = c.radius + height * slope;
                        if (v.0..=v.1).contains(&height) && rho > 0.0 {
                            let (sin, cos) = c.semi_angle.sin_cos();
                            found.push(e * cos - frame.z * sin);
                        }
                    }
                }
                Surface::Torus(t) => {
                    // (o + (R + r cos w) e + r sin w a - p) . n = 0.
                    let q = ((p - o).dot(n) - t.major_radius * en) / t.minor_radius;
                    let amplitude = en.hypot(an);
                    if q.abs() <= amplitude {
                        let phase = an.atan2(en);
                        for w in [
                            phase + (q / amplitude).acos(),
                            phase - (q / amplitude).acos(),
                        ] {
                            if (w - v.0).rem_euclid(TAU) <= v.1 - v.0 + 1e-12 {
                                found.push(e * w.cos() + frame.z * w.sin());
                            }
                        }
                    }
                }
                _ => unreachable!(),
            }
            for normal in found {
                least = least.min((1.0 - n.dot(normal).powi(2)).max(0.0).sqrt());
            }
        }
        least
    }

    #[test]
    fn the_closed_forms_agree_with_a_sweep_of_the_surface() {
        let whole = (0.0, TAU);
        // A cone cut by a tilted plane: `alpha` and `m` both live.
        let cone = Surface::Cone(Cone {
            frame: upright(),
            radius: 0.5,
            semi_angle: FRAC_PI_6,
        });
        let n = Vec3::new(0.866, 0.2, 0.5).normalize();
        let p = Point3::new(0.0, 0.0, 0.5);
        let closed = sine_of(p, n, &cone, (Point2::new(0.0, 0.0), Point2::new(TAU, 1.0)))
            .expect("a section");
        let swept_cone = swept(p, n, &cone, whole, (0.0, 1.0));
        assert!(
            closed <= swept_cone + 1e-9 && swept_cone - closed < 1e-3,
            "cone {closed} against {swept_cone}"
        );
        // And the bent tube, cut by planes leaning every way, in the
        // quarter turn and in a face cut short.
        for (n, p, ends) in [
            (
                Vec3::new(1.0, 0.5, 0.2),
                Point3::new(-0.2, 0.05, 0.1),
                (0.0, FRAC_PI_2),
            ),
            (
                Vec3::new(1.0, -0.3, 0.4),
                Point3::new(-0.3, 0.0, 0.2),
                (0.0, FRAC_PI_2),
            ),
            (
                Vec3::new(0.3, 0.1, 1.0),
                Point3::new(-0.4, 0.0, 0.5),
                (0.0, FRAC_PI_2),
            ),
            (
                Vec3::new(1.0, 0.5, 0.2),
                Point3::new(-0.2, 0.05, 0.1),
                (0.0, 0.9),
            ),
        ] {
            let tube = bent_tube();
            let closed = sine_of(
                p,
                n,
                &tube,
                (Point2::new(ends.0, 0.0), Point2::new(ends.1, TAU)),
            );
            let swept_tube = swept(p, n, &tube, ends, (0.0, TAU));
            match closed {
                Some(closed) => assert!(
                    closed <= swept_tube + 1e-9 && swept_tube - closed < 1e-3,
                    "torus {closed} against {swept_tube} for {n:?}"
                ),
                None => assert!(swept_tube.is_infinite(), "torus none against {swept_tube}"),
            }
        }
    }

    #[test]
    fn a_tangent_point_on_a_part_the_face_lacks_is_no_pair() {
        // R = 1, r = 0.25, a face over the outer half of the tube
        // (`phi` in `[-pi/2, pi/2]`): the plane tangent to the tube at
        // `phi = 1.75`, on the part the face lacks, still cuts the face at
        // its least sine of 0.178, where a form that reads the whole tube
        // finds the tangent contact and returns 0.
        let torus = Surface::Torus(Torus {
            frame: upright(),
            major_radius: 1.0,
            minor_radius: 0.25,
        });
        let n = Vec3::new(1.75_f64.cos(), 0.0, 1.75_f64.sin());
        let p = n * (1.75_f64.cos() + 0.25);
        let half = (Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2));
        let s = sine_of(p, n, &torus, half).expect("a section");
        let brute = swept(p, n, &torus, (0.0, TAU), (-FRAC_PI_2, FRAC_PI_2));
        assert!((s - 0.178).abs() < 1e-3, "{s}");
        assert!(s <= brute + 1e-9 && brute - s < 1e-3, "{s} against {brute}");
        // The whole tube has the tangent point.
        let whole = (Point2::new(0.0, 0.0), Point2::new(TAU, TAU));
        assert!(sine_of(p, n, &torus, whole).expect("a section") < 1e-6);
        // Several faces of the tube, every plane leaning: the sweep agrees.
        for (n, p, v) in [
            (
                Vec3::new(1.0, 0.3, 0.8),
                Point3::new(0.9, 0.0, 0.1),
                (-1.0, 1.0),
            ),
            (
                Vec3::new(0.4, 0.3, 1.0),
                Point3::new(0.2, 0.1, 0.3),
                (0.5, 2.5),
            ),
            (
                Vec3::new(1.0, 0.0, 0.1),
                Point3::new(1.1, 0.0, 0.0),
                (-0.5, 0.5),
            ),
        ] {
            let closed = sine_of(p, n, &torus, (Point2::new(0.0, v.0), Point2::new(TAU, v.1)));
            let brute = swept(p, n, &torus, (0.0, TAU), v);
            match closed {
                Some(closed) => assert!(
                    closed <= brute + 1e-9 && brute - closed < 1e-3,
                    "{closed} against {brute} for {n:?} over {v:?}"
                ),
                None => assert!(brute.is_infinite(), "none against {brute}"),
            }
        }
    }

    // --- a plane that touches a face, and one that cuts through it ------

    fn curved_face(surface: Surface, (lo, hi): (Point2, Point2)) -> Face {
        Face {
            shape: Shape::Curved { surface, lo, hi },
            ball: None,
        }
    }

    fn plane_face(point: Point3, normal: Vec3) -> Face {
        Face {
            shape: Shape::Plane { point, normal },
            ball: None,
        }
    }

    /// The factor of `surface` over its box against the plane, at a
    /// tolerance of a millimetre.
    fn factor_against(surface: Surface, bounds: (Point2, Point2), p: Point3, n: Vec3) -> Scalar {
        factor(&[curved_face(surface, bounds)], &[plane_face(p, n)], 1e-3)
    }

    /// A plane through the axis of a cylinder of radius `r` and length
    /// `len`, at its middle, `degrees` to the axis.
    fn through_the_axis(r: Scalar, len: Scalar, degrees: Scalar) -> Scalar {
        let tilt = degrees.to_radians();
        factor_against(
            column(r),
            (Point2::new(0.0, 0.0), Point2::new(TAU, len)),
            Point3::new(0.0, 0.0, 0.5 * len),
            Vec3::new(tilt.cos(), 0.0, tilt.sin()),
        )
    }

    #[test]
    fn a_plane_through_a_cylinder_at_a_slant_cuts_it_deep() {
        // The depth `radius (1 - cos)` of the section's sine is that of a
        // plane parallel to the axis; a slanted plane through the axis cuts
        // the whole radius, however small that sine reads. Each of these
        // read under the millimetre by it and kept their chord.
        for (r, len, degrees, shrunk) in [
            // A 3 cm pipe cut at 20 degrees: sine 0.342, "depth" 0.9 mm,
            // real 14 mm.
            (0.015, 1.0, 20.0, 87.0 / 256.0),
            // A 0.2 m column cut lengthwise at 5 degrees.
            (0.2, 6.0, 5.0, FLOOR),
            // A 3 mm rod mitred at 45 degrees.
            (0.003, 0.1, 45.0, 181.0 / 256.0),
        ] {
            assert_eq!(through_the_axis(r, len, degrees), shrunk, "{r} {degrees}");
        }
    }

    #[test]
    fn a_plane_reaching_a_cylinder_by_the_tolerance_only_touches_it() {
        // Parallel to the axis, a hair inside the radius: the reach is the
        // hair.
        let whole = (Point2::new(0.0, 0.0), Point2::new(TAU, 3.0));
        let at = |inside: Scalar| {
            factor_against(
                column(0.1),
                whole,
                Point3::new(0.1 - inside, 0.0, 0.0),
                Vec3::X,
            )
        };
        assert_eq!(at(0.0), 1.0);
        assert_eq!(at(9e-4), 1.0);
        assert_eq!(at(2e-3), FLOOR);
        // A slanted plane that clips only the end of the face: a sliver
        // of the rim, a hair deep, is no more than a touch either.
        let tilt = 0.3_f64;
        let n = Vec3::new(tilt.cos(), 0.0, tilt.sin());
        let clip = |hair: Scalar| {
            // The rim point `(0.1, 0, 3)` lies `hair` past the plane.
            let p = Point3::new(0.1, 0.0, 3.0) - n * hair;
            factor_against(column(0.1), whole, p, n)
        };
        assert_eq!(clip(5e-4), 1.0);
        assert!(clip(5e-2) < 1.0);
    }

    #[test]
    fn a_plane_cutting_a_tube_whole_has_the_tangent_inside_its_crossing() {
        // R = 1, r = 0.25 as one face, the tube of a revolved circle: the
        // plane tangent at `phi = 1.75` also crosses the tube elsewhere,
        // at a sine of 0.178, so it is a section and its tangent point
        // reads the floor.
        let whole = (Point2::new(0.0, 0.0), Point2::new(TAU, TAU));
        let torus = Surface::Torus(Torus {
            frame: upright(),
            major_radius: 1.0,
            minor_radius: 0.25,
        });
        let n = Vec3::new(1.75_f64.cos(), 0.0, 1.75_f64.sin());
        let p = n * (1.75_f64.cos() + 0.25);
        assert_eq!(factor_against(torus, whole, p, n), FLOOR);
        // The thin tube of the issue: a plane 0.9 from the axis touches
        // the inner equator and crosses the tube at 37 degrees either side.
        let thin = Surface::Torus(Torus {
            frame: upright(),
            major_radius: 1.0,
            minor_radius: 0.1,
        });
        assert_eq!(
            factor_against(thin, whole, Point3::new(0.9, 0.0, 0.0), Vec3::X),
            FLOOR
        );
    }

    #[test]
    fn a_plane_touching_a_tube_from_outside_is_no_pair() {
        let whole = (Point2::new(0.0, 0.0), Point2::new(TAU, TAU));
        let thin = Surface::Torus(Torus {
            frame: upright(),
            major_radius: 1.0,
            minor_radius: 0.1,
        });
        // Tangent at the outer equator, and reaching in by half a
        // millimetre: the tube lies to one side of it, or all but.
        for at in [1.1, 1.1 - 5e-4] {
            assert_eq!(
                factor_against(thin.clone(), whole, Point3::new(at, 0.0, 0.0), Vec3::X),
                1.0,
                "{at}"
            );
        }
        // Reaching in by a centimetre cuts it.
        assert_eq!(
            factor_against(thin, whole, Point3::new(1.09, 0.0, 0.0), Vec3::X),
            FLOOR
        );
    }

    #[test]
    fn a_plane_touching_a_sphere_is_no_pair_and_one_cutting_it_is() {
        let whole = (Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2));
        let sphere = Surface::Sphere(Sphere {
            frame: upright(),
            radius: 0.5,
        });
        for (at, factor) in [(0.5, 1.0), (0.5 - 5e-4, 1.0), (0.49, FLOOR), (0.0, 1.0)] {
            assert_eq!(
                factor_against(sphere.clone(), whole, Point3::new(0.0, 0.0, at), Vec3::Z),
                factor,
                "{at}"
            );
        }
    }

    #[test]
    fn a_cone_clipped_square_is_a_cap_and_a_slanted_cut_is_deep() {
        let cone = Surface::Cone(Cone {
            frame: upright(),
            radius: 0.5,
            semi_angle: FRAC_PI_6,
        });
        let span = (Point2::new(0.0, 0.0), Point2::new(TAU, 1.0));
        // Square to the axis at the base: the cone lies to one side.
        assert_eq!(
            factor_against(cone.clone(), span, Point3::ZERO, Vec3::Z),
            1.0
        );
        // Through its middle, tilted by a degree: the full width on each
        // side.
        let tilt = 1.0_f64.to_radians();
        let n = Vec3::new(tilt.sin(), 0.0, tilt.cos());
        let s = sine_of(Point3::new(0.0, 0.0, 0.5), n, &cone, span).expect("a section");
        assert_eq!(
            factor_against(cone, span, Point3::new(0.0, 0.0, 0.5), n),
            quantized(s)
        );
    }

    #[test]
    fn the_reach_agrees_with_the_points_of_the_surface() {
        // The reach of every family against planes leaning every way, read
        // off a grid of the surface's own points: the closed form holds all
        // of them and is no more than the grid's step above.
        let cases: Vec<(Surface, (Point2, Point2))> = vec![
            (column(0.2), (Point2::new(0.0, 0.0), Point2::new(TAU, 1.5))),
            (
                Surface::Cone(Cone {
                    frame: upright(),
                    radius: 0.5,
                    semi_angle: FRAC_PI_6,
                }),
                (Point2::new(0.0, -0.4), Point2::new(TAU, 1.0)),
            ),
            (
                Surface::Sphere(Sphere {
                    frame: upright(),
                    radius: 0.5,
                }),
                (Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2)),
            ),
            (bent_tube(), quarter()),
            (bent_tube(), (Point2::new(0.3, 1.0), Point2::new(1.2, 4.0))),
            (
                Surface::Torus(Torus {
                    frame: upright(),
                    major_radius: 1.0,
                    minor_radius: 0.25,
                }),
                (Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2)),
            ),
            (
                Surface::Torus(Torus {
                    frame: upright(),
                    major_radius: 1.0,
                    minor_radius: 0.25,
                }),
                (Point2::new(0.0, 0.0), Point2::new(TAU, TAU)),
            ),
        ];
        let planes = [
            (Point3::new(0.1, 0.05, 0.2), Vec3::new(1.0, 0.3, 0.8)),
            (Point3::new(-0.3, 0.1, 0.4), Vec3::new(0.2, -0.3, 1.0)),
            (Point3::new(0.9, 0.0, 0.1), Vec3::new(1.0, 0.0, 0.05)),
            (Point3::new(0.0, 0.0, 0.6), Vec3::new(0.1, 0.0, 1.0)),
            (Point3::new(0.15, -0.2, 0.0), Vec3::new(0.7, 0.7, 0.1)),
            (Point3::new(-0.2, 0.0, 0.0), Vec3::new(1.0, 0.5, 0.2)),
        ];
        let steps = 400;
        for (surface, (lo, hi)) in &cases {
            for &(p, n) in &planes {
                let n = n.normalize();
                let (mut least, mut greatest) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
                for i in 0..=steps {
                    for j in 0..=steps {
                        let u = lo.x + (hi.x - lo.x) * f64::from(i) / f64::from(steps);
                        let v = lo.y + (hi.y - lo.y) * f64::from(j) / f64::from(steps);
                        let at =
                            axiolid_reference::surface::evaluate(surface, u, v).expect("a point");
                        let d = n.dot(at - p);
                        least = least.min(d);
                        greatest = greatest.max(d);
                    }
                }
                let brute = greatest.max(0.0).min((-least).max(0.0));
                let closed = shallow_reach(p, n, surface, *lo, *hi).expect("a reach");
                assert!(
                    closed >= brute - 1e-12 && closed - brute < 2e-3,
                    "{closed} against {brute} for {n:?} at {p:?}"
                );
            }
        }
    }

    #[test]
    fn a_plane_past_a_cones_apex_misses_it() {
        // The cone `r = 0.5 + v tan(30)` ends at its apex, `v = -0.866`: a
        // plane square to its axis below that has no section, though the
        // unbounded surface's box reaches it.
        let cone = Surface::Cone(Cone {
            frame: upright(),
            radius: 0.5,
            semi_angle: FRAC_PI_6,
        });
        let box_ = (Point2::new(0.0, -2.0), Point2::new(TAU, 1.0));
        assert_eq!(
            sine_of(Point3::new(0.0, 0.0, -1.5), Vec3::Z, &cone, box_),
            None
        );
        assert!(sine_of(Point3::new(0.0, 0.0, -0.5), Vec3::Z, &cone, box_).is_some());
    }

    // --- the factors of a boolean --------------------------------------

    fn options() -> ExecutionOptions {
        ExecutionOptions::new(Tolerance::new(1e-6, 1e-9).unwrap())
            .with_chord_error(1e-3)
            .unwrap()
    }

    struct Fixture(GeometryGraphBuilder);

    impl Fixture {
        fn new() -> Self {
            Self(GeometryGraphBuilder::new())
        }

        fn push(&mut self, node: GeometryNode) -> NodeId {
            self.0.push(node).expect("a valid node")
        }

        fn extrusion(&mut self, profile: Profile, depth: Scalar) -> NodeId {
            let profile = self.push(GeometryNode::Profile(profile));
            self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth,
            }))
        }

        fn rect(&mut self, x: Scalar, y: Scalar, depth: Scalar) -> NodeId {
            self.extrusion(
                Profile::Rectangle(RectangleProfile {
                    x,
                    y,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                }),
                depth,
            )
        }

        fn circle(&mut self, radius: Scalar, depth: Scalar) -> NodeId {
            self.extrusion(
                Profile::Circle(CircleProfile {
                    radius,
                    thickness: None,
                }),
                depth,
            )
        }

        fn place(&mut self, source: NodeId, transform: Transform3) -> NodeId {
            self.push(GeometryNode::Instance(Instance { source, transform }))
        }

        fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
            self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }))
        }

        fn half_space(&mut self, origin: Point3, normal: Vec3) -> NodeId {
            self.push(GeometryNode::HalfSpace(axiolid_primitive::HalfSpace {
                boundary: Plane3 { origin, normal },
                agreement: true,
            }))
        }

        fn finish(self, root: NodeId) -> GeometryGraph {
            self.0.finish(vec![root]).expect("a valid graph")
        }
    }

    fn building() -> Transform3 {
        Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
    }

    fn factors(graph: &GeometryGraph, root: NodeId) -> Factors {
        operand_factors(graph, root, &options(), &mut LeafCache::new())
    }

    /// The issue's repro (#378): a box and a quarter-turn bent tube.
    fn box_and_tube(placed: bool) -> (GeometryGraph, NodeId) {
        let mut g = Fixture::new();
        let a = g.rect(0.4, 0.5, 6.0);
        let profile = g.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius: 0.1,
            thickness: None,
        })));
        let b = g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin: Point3::new(-1.0, 0.0, 0.0),
            axis_direction: Vec3::NEG_Y,
            angle: FRAC_PI_2,
        }));
        let (a, b) = if placed {
            (g.place(a, building()), g.place(b, building()))
        } else {
            (a, b)
        };
        let root = g.boolean(a, b, BooleanOperator::Intersection);
        (g.finish(root), root)
    }

    #[test]
    fn a_box_and_a_bent_tube_shrink_the_tube_only() {
        let (graph, root) = box_and_tube(false);
        assert_eq!(
            factors(&graph, root),
            Factors {
                left: 1.0,
                right: 0.45703125
            }
        );
    }

    #[test]
    fn placing_both_operands_leaves_the_factors_unchanged() {
        let (graph, root) = box_and_tube(true);
        assert_eq!(
            factors(&graph, root),
            Factors {
                left: 1.0,
                right: 0.45703125
            }
        );
    }

    #[test]
    fn a_wall_with_a_round_window_keeps_its_chords() {
        let mut g = Fixture::new();
        let wall = g.rect(6.0, 0.3, 3.0);
        let wall = g.place(wall, building());
        let window = g.circle(0.4, 1.0);
        let across = Transform3::from_translation(Vec3::new(0.8, 0.5, 1.5))
            * Transform3::from_rotation_x(FRAC_PI_2);
        let window = g.place(window, building() * across);
        let root = g.boolean(wall, window, BooleanOperator::Difference);
        let graph = g.finish(root);
        assert_eq!(factors(&graph, root), Factors::ONE);
    }

    /// A block `1 x 1 x 1` and a column of radius `0.1` standing on its
    /// base, reaching `inside` past the block's face `x = 0.5`.
    fn block_and_column(inside: Scalar) -> (GeometryGraph, NodeId) {
        let mut g = Fixture::new();
        let block = g.rect(1.0, 1.0, 1.0);
        let column = g.circle(0.1, 3.0);
        let column = g.place(
            column,
            Transform3::from_translation(Vec3::new(0.5 - 0.1 + inside, 0.0, 0.0)),
        );
        let root = g.boolean(block, column, BooleanOperator::Difference);
        (g.finish(root), root)
    }

    #[test]
    fn a_column_touching_a_face_keeps_its_chord_and_one_grazing_it_gets_the_floor() {
        // Tangent to the face, within the rounding of its placement: no
        // section, and the mesh cut's tangent handling is left as it was.
        let (graph, root) = block_and_column(0.0);
        assert_eq!(factors(&graph, root), Factors::ONE);
        // A hundredth of the radius inside it, the column cuts a shallow
        // sliver at 8 degrees.
        let (graph, root) = block_and_column(0.001);
        assert_eq!(
            factors(&graph, root),
            Factors {
                left: 1.0,
                right: FLOOR
            }
        );
    }

    #[test]
    fn a_graze_within_the_tolerance_is_no_section() {
        // The column reaches 0.5 mm into the face: a depth of 5e-4, which a
        // tolerance of a millimetre holds and one of a micrometre does not.
        let (graph, root) = block_and_column(0.0005);
        let coarse = ExecutionOptions::new(Tolerance::new(1e-3, 1e-9).unwrap());
        assert_eq!(
            operand_factors(&graph, root, &coarse, &mut LeafCache::new()),
            Factors::ONE
        );
        assert_eq!(factors(&graph, root).right, FLOOR);
    }

    /// The beam of #243 whose round web hole touches its top flange, under
    /// `building`.
    fn beam_and_tangent_hole(building: Transform3) -> (GeometryGraph, NodeId) {
        let mut g = Fixture::new();
        let (depth, web, flange, radius) = (0.3, 0.0071, 0.0107, 0.05);
        let beam = g.extrusion(
            Profile::Section(axiolid_profile::SectionProfile::I {
                depth,
                width: 0.15,
                web_thickness: web,
                flange_thickness: flange,
                fillet_radius: None,
                flange_edge_radius: None,
                flange_slope: None,
            }),
            4.0,
        );
        let beam = g.place(beam, building);
        let hole = g.circle(radius, web);
        let across =
            Transform3::from_translation(Vec3::new(-web / 2.0, depth / 2.0 - flange - radius, 2.0))
                * Transform3::from_rotation_y(FRAC_PI_2);
        let hole = g.place(hole, building * across);
        let root = g.boolean(beam, hole, BooleanOperator::Difference);
        (g.finish(root), root)
    }

    #[test]
    fn a_tangent_hole_keeps_its_chord_at_any_placement() {
        // The rounding of a placement far from the origin turns a designed
        // tangent into a graze of a few thousandths of a radian; at a
        // millimetre tolerance it is still no section, at every scale.
        let options = ExecutionOptions::new(Tolerance::new(1e-3, 1e-9).unwrap())
            .with_chord_error(1e-3)
            .unwrap();
        for at in [
            Vec3::new(12.5, -4.0, 3.2),
            Vec3::new(250.0, -400.0, 3.2),
            Vec3::new(1000.0, 2000.0, 10.0),
            Vec3::new(3.6e5, 5.8e6, 120.0),
        ] {
            for k in 0..12 {
                let turn = Transform3::from_rotation_z(0.6 + 0.05 * f64::from(k));
                let (graph, root) = beam_and_tangent_hole(Transform3::from_translation(at) * turn);
                assert_eq!(
                    operand_factors(&graph, root, &options, &mut LeafCache::new()),
                    Factors::ONE,
                    "{at:?} turn {k}"
                );
            }
        }
    }

    #[test]
    fn a_wall_with_rectangular_doors_is_read_without_a_construction() {
        let mut g = Fixture::new();
        let wall = g.rect(6.0, 0.3, 3.0);
        let wall = g.place(wall, building());
        let door = g.rect(0.9, 1.0, 2.1);
        let doors = g.boolean(wall, door, BooleanOperator::Difference);
        let other = g.rect(0.9, 1.0, 2.1);
        let other = g.place(
            other,
            Transform3::from_translation(Vec3::new(2.0, 0.0, 0.0)),
        );
        let root = g.boolean(doors, other, BooleanOperator::Difference);
        let graph = g.finish(root);
        let GeometryNode::SolidOperation(SolidOperation::Boolean { left, right, .. }) =
            graph.get(root).unwrap()
        else {
            unreachable!()
        };
        assert!(certainly_planar(&graph, *left));
        assert!(certainly_planar(&graph, *right));
        // No leaf was compiled.
        let mut leaves = LeafCache::new();
        assert_eq!(
            operand_factors(&graph, root, &options(), &mut leaves),
            Factors::ONE
        );
        assert!(leaves.is_empty());
    }

    #[test]
    fn a_round_profile_is_not_certainly_planar() {
        let mut g = Fixture::new();
        let column = g.circle(0.2, 3.0);
        let graph = g.finish(column);
        assert!(!certainly_planar(&graph, column));
    }

    fn contour(segments: Vec<axiolid_profile::ProfileSegment>) -> Profile {
        Profile::Contour(axiolid_profile::ContourProfile {
            outer: axiolid_profile::Contour::new(segments),
            holes: Vec::new(),
        })
    }

    fn line(from: Point2, to: Point2) -> axiolid_profile::ProfileSegment {
        axiolid_profile::ProfileSegment {
            curve: Curve2::Line(axiolid_curve::Line2 {
                origin: from,
                direction: to - from,
            }),
            domain: axiolid_core::Interval::UNIT,
            same_sense: true,
        }
    }

    #[test]
    fn only_straight_profiles_are_certainly_planar() {
        let (a, b, c) = (
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        );
        let mut g = Fixture::new();
        let triangle = g.extrusion(contour(vec![line(a, b), line(b, c), line(c, a)]), 1.0);
        let rounded = g.extrusion(
            Profile::Rectangle(RectangleProfile {
                x: 1.0,
                y: 1.0,
                thickness: None,
                outer_radius: Some(0.1),
                inner_radius: None,
            }),
            1.0,
        );
        let arc = axiolid_profile::ProfileSegment {
            curve: Curve2::Circle(axiolid_curve::Circle2 {
                frame: axiolid_core::Frame2 {
                    origin: Point2::ZERO,
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: 1.0,
            }),
            domain: axiolid_core::Interval::new(0.0, 1.0),
            same_sense: true,
        };
        let bulged = g.extrusion(contour(vec![arc, line(b, a)]), 1.0);
        let graph = g.finish(bulged);
        assert!(certainly_planar(&graph, triangle));
        assert!(!certainly_planar(&graph, rounded));
        assert!(!certainly_planar(&graph, bulged));
    }

    #[test]
    fn nested_placements_compose_in_order() {
        // A tube placed by one turn, then by a move, is the tube placed by
        // the two composed: the factors agree, and are not a no-op.
        let turn = Transform3::from_rotation_z(0.3);
        let slide = Transform3::from_translation(Vec3::new(0.1, 0.0, 0.0));
        let build = |nested: bool| {
            let mut g = Fixture::new();
            let a = g.rect(0.4, 0.5, 6.0);
            let profile = g.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
                radius: 0.1,
                thickness: None,
            })));
            let tube = g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
                profile,
                axis_origin: Point3::new(-1.0, 0.0, 0.0),
                axis_direction: Vec3::NEG_Y,
                angle: FRAC_PI_2,
            }));
            let tube = if nested {
                let first = g.place(tube, turn);
                g.place(first, slide)
            } else {
                g.place(tube, slide * turn)
            };
            let root = g.boolean(a, tube, BooleanOperator::Intersection);
            let graph = g.finish(root);
            factors(&graph, root)
        };
        let flat = build(false);
        assert_eq!(build(true), flat);
        assert!(flat.right < 1.0, "{flat:?}");
    }

    #[test]
    fn a_column_clipped_obliquely_shrinks_the_subject() {
        let tilt = FRAC_PI_6 * 2.0;
        let mut g = Fixture::new();
        let column = g.circle(0.2, 3.0);
        let roof = g.half_space(
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(tilt.sin(), 0.0, tilt.cos()),
        );
        let root = g.boolean(column, roof, BooleanOperator::Difference);
        let graph = g.finish(root);
        assert_eq!(
            factors(&graph, root),
            Factors {
                left: 0.5,
                right: 1.0
            }
        );
        // An intersection commutes: the region on the left reads the same.
        let mut g = Fixture::new();
        let column = g.circle(0.2, 3.0);
        let roof = g.half_space(
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(tilt.sin(), 0.0, tilt.cos()),
        );
        let root = g.boolean(roof, column, BooleanOperator::Intersection);
        let graph = g.finish(root);
        assert_eq!(
            factors(&graph, root),
            Factors {
                left: 1.0,
                right: 0.5
            }
        );
    }

    #[test]
    fn a_half_space_under_a_stretch_keeps_its_plane() {
        // The plane `x + z = 0` stretched twice along `z` is `x + z / 2 = 0`:
        // a normal moves by the inverse transpose, not by the matrix.
        let stretch = Transform3::from_scale(Vec3::new(1.0, 1.0, 2.0));
        let half_space = axiolid_primitive::HalfSpace {
            boundary: Plane3 {
                origin: Point3::ZERO,
                normal: Vec3::new(1.0, 0.0, 1.0),
            },
            agreement: true,
        };
        let face = half_space_plane(&half_space, &stretch).expect("a plane");
        let Shape::Plane { point, normal } = face.shape else {
            unreachable!()
        };
        assert!((normal.length() - 1.0).abs() < 1e-12);
        for on in [
            Point3::new(1.0, 0.0, -1.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(2.0, 5.0, -2.0),
        ] {
            let off = normal.dot(stretch.transform_point3(on) - point);
            assert!(off.abs() < 1e-12, "{off}");
        }
        // No inverse, no plane.
        let flat = Transform3::from_scale(Vec3::new(1.0, 1.0, 0.0));
        assert!(half_space_plane(&half_space, &flat).is_none());
    }

    #[test]
    fn a_pair_whose_balls_do_not_meet_is_not_read() {
        // The tube far from the box: the box's plane `x = -0.2` meets the
        // tube's torus unbounded, but not within the box's face.
        let mut g = Fixture::new();
        let a = g.rect(0.4, 0.5, 6.0);
        let a = g.place(a, Transform3::from_translation(Vec3::new(0.0, 0.0, 40.0)));
        let profile = g.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius: 0.1,
            thickness: None,
        })));
        let b = g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin: Point3::new(-1.0, 0.0, 0.0),
            axis_direction: Vec3::NEG_Y,
            angle: FRAC_PI_2,
        }));
        let root = g.boolean(a, b, BooleanOperator::Union);
        let graph = g.finish(root);
        assert_eq!(factors(&graph, root), Factors::ONE);
    }
}
